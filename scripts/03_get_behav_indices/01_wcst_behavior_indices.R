# stan_data <- readRDS(here::here("data", "processed", "wcst_stan_list.RDS"))
#
# saveRDS(
#   stan_data,
#   here::here(
#     "data",
#     "processed",
#     "interim",
#     "wcst_stan_data.RDS"
#   )
# )

# ============================================================
# WCST — Indici comportamentali da lista RDS (input Stan)
# Autore: ChatGPT
# Data: generato automaticamente
# ------------------------------------------------------------
# Cosa calcola (per soggetto):
# - trials                : numero di prove
# - total_errors          : errori totali
# - perseverative_errors  : errori perseverativi (PE)
# - non_persev_errors     : errori non perseverativi (NPE)
# - perseverative_resp    : risposte perseverative (PR, incl. corrette)
# - categories_achieved   : categorie completate (se rilevabile)
# - failure_maintain_set  : fallimenti nel mantenere il set (FMS, se rilevabile)
# - Proporzioni su prove: *_prop_trials  (es. PE / trials)
# - Proporzioni su errori: *_prop_errors (es. PE / total_errors)
#
# Requisiti del .RDS:
# - Il file RDS deve contenere dati per 1+ soggetti in una di queste forme:
#   A) LONG DATA FRAME per prove con colonne (nomi flessibili, vedi mapping):
#       user_id, group, trial, choice, correct, rule, a_color, a_form, a_number
#      (rule = 1:colore, 2:forma, 3:num; correct ∈ {0,1})
#      a_color/a_form/a_number = indice (1..4) della chiave corretta per ciascuna regola in quel trial.
#   B) LISTA in stile Stan, con matrici o liste per soggetto, ad es.:
#       S (n.soggetti), T (n.trial), y (SxT scelte 1..4), r (SxT 0/1),
#       rule (SxT 1..3, opz.), a_color/a_form/a_number (SxT, opz.),
#       user_id (lunghezza S), group (lunghezza S)
#
# Se le colonne a_* non sono disponibili, lo script fornirà indici "generali"
# (errori totali, win-stay/lose-shift, ecc.), ma non PE/PR "classici" Heaton.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})

# ------------------------------
# 1) Parametri di input/output
# ------------------------------

# Percorso del file RDS (modifica se necessario)
path_rds <- here::here(
  "data",
  "processed",
  "interim",
  "wcst_stan_data.RDS"
)

# File di output con gli indici per-soggetto
out_csv <- here::here(
  "data",
  "processed",
  "wcst_behavioral_indices.csv"
)

# ------------------------------
# 2) Utility: stampa struttura
# ------------------------------
peek_obj <- function(x) {
  cat("CLASS:", paste(class(x), collapse = " "), "\n")
  if (is.list(x)) {
    cat("LUNGHEZZA LISTA:", length(x), "\n")
    nms <- names(x)
    if (!is.null(nms)) {
      cat("NOMI LISTA:", paste(head(nms, 30), collapse = ", "), "\n")
    }
  }
  if (is.data.frame(x)) {
    cat("COLONNE DF:", paste(names(x), collapse = ", "), "\n")
  }
  invisible(NULL)
}

# ------------------------------
# 3) Normalizzazione a long DF
# ------------------------------

# Heuristics per mappare nomi flessibili -> standard
standardize_names <- function(df) {
  nm <- names(df)
  map_to <- c(
    user_id = nm[str_detect(
      nm,
      regex("^user(_|)id$|^subj(_|)id$|^id$", ignore_case = TRUE)
    )][1],
    group = nm[str_detect(
      nm,
      regex("^group$|^grp$|^condition$|^arm$", ignore_case = TRUE)
    )][1],
    trial = nm[str_detect(
      nm,
      regex("^trial$|^t$|^time$|^trial_index$", ignore_case = TRUE)
    )][1],
    choice = nm[str_detect(
      nm,
      regex("^choice$|^resp$|^y$|^action$", ignore_case = TRUE)
    )][1],
    correct = nm[str_detect(
      nm,
      regex("^correct$|^reward$|^r$|^feedback$", ignore_case = TRUE)
    )][1],
    rule = nm[str_detect(
      nm,
      regex("^rule$|^cat$|^category$|^dimension$", ignore_case = TRUE)
    )][1],
    a_color = nm[str_detect(
      nm,
      regex("a(_|)color|ans(_|)color|opt(_|)color", ignore_case = TRUE)
    )][1],
    a_form = nm[str_detect(
      nm,
      regex("a(_|)form|ans(_|)form|opt(_|)form|shape", ignore_case = TRUE)
    )][1],
    a_number = nm[str_detect(
      nm,
      regex("a(_|)number|ans(_|)number|opt(_|)number|count", ignore_case = TRUE)
    )][1]
  )
  # Rinominare solo quelli trovati
  map_to <- map_to[!is.na(map_to)]
  dplyr::rename(df, !!!setNames(names(map_to), map_to))
}

# Da formato lista "stan-like" a long DF
list_to_long <- function(L) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # Deduzione S e T
  S <- L$S %||% L$N %||% L$n_subj %||% length(L$user_id)
  Tt <- L$T %||%
    L$Tmax %||%
    L$n_trial %||%
    (if (!is.null(L$resp_choice)) dim(L$resp_choice)[2] else NULL)
  if (is.null(S) || is.null(Tt)) {
    stop("Impossibile dedurre S (n soggetti) e T (n trial) dalla lista.")
  }

  # ID e gruppo
  user_id <- L$user_id %||% seq_len(S)
  group <- L$group %||% rep(NA, S)

  # === MAPPATURA SPECIFICA PER IL TUO FILE ===
  # scelte e correttezza
  y <- L$resp_choice %||% L$choice %||% L$y
  r <- L$rew %||% L$reward %||% L$correct %||% L$r

  if (is.null(y) || is.null(r)) {
    stop("Non trovo le matrici delle scelte/correttezza (resp_choice/rew).")
  }
  if (!all(dim(y) == c(S, Tt))) stop("Dimensioni inattese per resp_choice.")
  if (!all(dim(r) == c(S, Tt))) stop("Dimensioni inattese per rew.")

  # regola e target-by-rule (per PR/PE)
  rule <- L$rule_choice %||% L$rule %||% L$category %||% L$cat
  a_color <- L$resp_color
  a_form <- L$resp_shape
  a_number <- L$resp_number

  base <- tibble::tibble(
    user_id = rep(user_id, each = Tt),
    group = rep(group, each = Tt),
    trial = rep(seq_len(Tt), times = S),
    choice = as.integer(as.vector(t(y))),
    correct = as.integer(as.vector(t(r)))
  )

  add_if <- function(df, mat, name) {
    if (!is.null(mat)) {
      if (!all(dim(mat) == c(S, Tt)))
        stop(paste("Dimensioni inattese per", name))
      df[[name]] <- as.integer(as.vector(t(mat)))
    }
    df
  }

  base <- add_if(base, rule, "rule")
  base <- add_if(base, a_color, "a_color")
  base <- add_if(base, a_form, "a_form")
  base <- add_if(base, a_number, "a_number")

  base
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

to_long_df <- function(obj) {
  if (is.data.frame(obj)) {
    df <- standardize_names(obj)
    # Controlli minimi
    need <- c("user_id", "trial", "choice", "correct")
    if (!all(need %in% names(df)))
      stop(
        "Data frame privo di colonne minime: ",
        paste(setdiff(need, names(df)), collapse = ", ")
      )
    # Assicurare group se manca
    if (!"group" %in% names(df)) df$group <- NA
    df <- df %>% arrange(user_id, trial)
    return(df)
  }
  if (is.list(obj)) {
    return(list_to_long(obj))
  }
  stop("Oggetto RDS non riconosciuto (né data.frame né lista).")
}

# ------------------------------
# 4) Scoring WCST
# ------------------------------

# Identifica i trial di cambio regola (shift)
flag_shifts <- function(df) {
  if (!"rule" %in% names(df)) {
    df <- df %>%
      group_by(user_id) %>%
      mutate(
        shift = if_else(
          lag(correct, default = 1L) == 1L & correct == 0L,
          0L,
          0L
        )
      ) %>%
      ungroup()
    df$shift <- NA_integer_
    return(df)
  }
  df %>%
    group_by(user_id) %>%
    arrange(trial, .by_group = TRUE) %>%
    mutate(
      shift = if_else(is.na(lag(rule)), 1L, as.integer(rule != lag(rule)))
    ) %>%
    replace_na(list(shift = 0L)) %>%
    ungroup()
}

# PE/PR basate su conoscenza del "target per regola"
# Definizione:
# - PR (perseverative response): una risposta che segue la regola precedente (prev_rule),
#   ossia choice == a_{prev_rule} al trial corrente.
# - PE (perseverative error): PR & incorrect == 1.
# - NPE (non-perseverative error): incorrect == 1 & non PR.
classify_perseveration <- function(df, rule_col = "rule_rec") {
  have_targets <- all(c("a_color", "a_form", "a_number") %in% names(df))
  if (!have_targets || !(rule_col %in% names(df))) {
    df$PR <- NA_integer_
    df$PE <- NA_integer_
    df$NPE <- NA_integer_
    return(df)
  }
  df %>%
    group_by(user_id) %>%
    arrange(trial, .by_group = TRUE) %>%
    mutate(
      prev_rule = dplyr::lag(.data[[rule_col]]),
      # Identifica quando c'è uno shift di regola
      rule_switched = !is.na(prev_rule) &
        !is.na(.data[[rule_col]]) &
        (.data[[rule_col]] != prev_rule),
      # Target della regola PRECEDENTE (prima dello shift)
      prev_target = case_when(
        prev_rule == 1L ~ a_color,
        prev_rule == 2L ~ a_form,
        prev_rule == 3L ~ a_number,
        TRUE ~ NA_integer_
      ),
      # PR: risposta che matcha il target della regola precedente SOLO dopo uno shift
      PR = as.integer(
        rule_switched & !is.na(prev_target) & choice == prev_target
      ),
      # PE: errore perseverativo (PR che è anche errore)
      PE = as.integer(PR == 1L & correct == 0L),
      # NPE: errore non perseverativo (errore che non è PR)
      NPE = as.integer(correct == 0L & PR == 0L)
    ) %>%
    select(-prev_rule, -rule_switched, -prev_target) %>%
    ungroup()
}


# Failure to Maintain Set (FMS) — euristico:
# errore isolato all'interno di una sequenza di ≥5 corrette senza cambio regola
flag_fms <- function(df, run_len = 5L) {
  if (!"rule" %in% names(df)) {
    df$FMS <- NA_integer_
    return(df)
  }
  df %>%
    group_by(user_id, rule) %>%
    arrange(trial, .by_group = TRUE) %>%
    mutate(
      correct_run = if_else(correct == 1L, 1L, 0L),
      correct_run = ave(correct_run, user_id, rule, FUN = function(x) {
        # calcolo lunghezze delle run di corretti (base R)
        r <- integer(length(x))
        run_id <- cumsum(c(1L, diff(x != 1L)))
        tapply(seq_along(x), run_id, function(idx) {
          if (x[idx[1]] == 1L) r[idx] <<- seq_along(idx) else r[idx] <<- 0L
        })
        r
      })
    ) %>%
    group_by(user_id) %>%
    mutate(
      FMS = as.integer(
        correct == 0L &
          lag(correct_run, default = 0L) >= run_len &
          (is.na(lag(rule)) | rule == lag(rule))
      )
    ) %>%
    ungroup() %>%
    select(-correct_run)
}

# Win–Stay / Lose–Shift (utili come indici generali se non si possono calcolare PE/PR)
compute_wsls <- function(df) {
  df %>%
    group_by(user_id) %>%
    arrange(trial, .by_group = TRUE) %>%
    mutate(
      stay = as.integer(choice == lag(choice)),
      win_stay = as.integer(lag(correct) == 1L & stay == 1L),
      lose_shift = as.integer(lag(correct) == 0L & stay == 0L)
    ) %>%
    summarise(
      n_ws = sum(win_stay, na.rm = TRUE),
      n_ws_denom = sum(lag(correct) == 1L, na.rm = TRUE),
      p_win_stay = ifelse(n_ws_denom > 0, n_ws / n_ws_denom, NA_real_),

      n_ls = sum(lose_shift, na.rm = TRUE),
      n_ls_denom = sum(lag(correct) == 0L, na.rm = TRUE),
      p_lose_shift = ifelse(n_ls_denom > 0, n_ls / n_ls_denom, NA_real_)
    ) %>%
    ungroup()
}

# Categorie completate (euristico): conta i blocchi di regola con ≥10 corretti consecutivi
categories_achieved_df <- function(df, criterion = 10L) {
  if (!"rule" %in% names(df)) {
    return(
      df %>%
        group_by(user_id) %>%
        summarise(categories_achieved = NA_integer_) %>%
        ungroup()
    )
  }
  df %>%
    group_by(user_id, rule) %>%
    arrange(trial, .by_group = TRUE) %>%
    summarise(
      max_streak = {
        x <- correct
        best <- 0L
        cur <- 0L
        for (v in x) {
          if (isTRUE(v == 1L)) {
            cur <- cur + 1L
            if (cur > best) best <- cur
          } else cur <- 0L
        }
        best
      },
      .groups = "drop_last"
    ) %>%
    summarise(
      categories_achieved = sum(max_streak >= criterion, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    ungroup()
}

# Wrapper: calcola tutti gli indici richiesti
score_wcst <- function(long_df) {
  df <- long_df %>%
    mutate(
      correct = as.integer(correct),
      choice = as.integer(choice),
      trial = as.integer(trial)
    )

  df <- flag_shifts(df)
  df <- classify_perseveration(df)
  df <- flag_fms(df)

  # Aggregazione per soggetto
  agg_basic <- df %>%
    group_by(user_id, group) %>%
    summarise(
      trials = n(),
      total_errors = sum(correct == 0L, na.rm = TRUE),
      perseverative_errors = sum(PE == 1L, na.rm = TRUE),
      non_persev_errors = sum(NPE == 1L, na.rm = TRUE),
      perseverative_resp = sum(PR == 1L, na.rm = TRUE),
      fms_count = sum(FMS == 1L, na.rm = TRUE),
      .groups = "drop"
    )

  # Proporzioni
  agg_basic <- agg_basic %>%
    mutate(
      pe_prop_trials = ifelse(
        trials > 0,
        perseverative_errors / trials,
        NA_real_
      ),
      npe_prop_trials = ifelse(
        trials > 0,
        non_persev_errors / trials,
        NA_real_
      ),
      pr_prop_trials = ifelse(
        trials > 0,
        perseverative_resp / trials,
        NA_real_
      ),
      te_prop_trials = ifelse(trials > 0, total_errors / trials, NA_real_),

      pe_prop_errors = ifelse(
        total_errors > 0,
        perseverative_errors / total_errors,
        NA_real_
      ),
      npe_prop_errors = ifelse(
        total_errors > 0,
        non_persev_errors / total_errors,
        NA_real_
      )
    )

  # WSLS
  wsls <- compute_wsls(df)

  # Categorie completate
  cats <- categories_achieved_df(df)

  # Merge
  out <- agg_basic %>%
    left_join(wsls, by = "user_id") %>%
    left_join(cats, by = "user_id") %>%
    relocate(group, .after = user_id)

  # Ordinamento
  out %>% arrange(group, user_id)
}

# ------------------------------
# 5) Esecuzione
# ------------------------------

rds_obj <- readRDS(path_rds)
peek_obj(rds_obj)
long_df <- to_long_df(rds_obj)

indices <- score_wcst(long_df)

# Salvataggio CSV
readr::write_csv(indices, out_csv)
message("Salvato: ", out_csv)

# Stampa anteprima
print(dplyr::glimpse(indices))
print(head(indices, 10))


###########

old <- rio::import(
  here::here("data", "processed", "wcst_behav_indices.csv")
) |>
  dplyr::rename(
    user_id = subj_name,
    gr = group,
    old_prop_pers_err = prop_pers_err,
    old_prop_non_pers_err = prop_non_pers_err,
    old_prop_pers_resp = prop_pers_resp
  )

tot <- full_join(indices, old, by = "user_id")

plot(tot$prop_non_pers_err, tot$old_prop_non_pers_err)
