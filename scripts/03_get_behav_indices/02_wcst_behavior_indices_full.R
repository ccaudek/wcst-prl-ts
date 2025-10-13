# ============================================================
# WCST — Indici comportamentali da lista RDS (input Stan)
# Versione: patch "rule recovery" + PR solo dopo shift
# Autore: ChatGPT
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  # library(here) # scommenta se vuoi usare here::here(...)
})

# ------------------------------
# 1) Parametri di input/output
# ------------------------------

# Percorso del file RDS
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
    if (!is.null(nms))
      cat("NOMI LISTA:", paste(head(nms, 30), collapse = ", "), "\n")
  }
  if (is.data.frame(x))
    cat("COLONNE DF:", paste(names(x), collapse = ", "), "\n")
  invisible(NULL)
}

# ------------------------------
# 3) Normalizzazione a long DF
# ------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Da formato lista "stan-like" a long DF — adattato al tuo RDS
list_to_long <- function(L) {
  S <- L$S %||% L$N %||% L$n_subj %||% length(L$user_id)
  Tt <- L$T %||%
    L$Tmax %||%
    L$n_trial %||%
    (if (!is.null(L$resp_choice)) dim(L$resp_choice)[2] else NULL)
  if (is.null(S) || is.null(Tt))
    stop("Impossibile dedurre S (n soggetti) e T (n trial) dalla lista.")

  user_id <- L$user_id %||% seq_len(S)
  group <- L$group %||% rep(NA, S)

  # scelte e correttezza
  y <- L$resp_choice %||% L$choice %||% L$y
  r <- L$rew %||% L$reward %||% L$correct %||% L$r
  if (is.null(y) || is.null(r))
    stop("Non trovo le matrici delle scelte/correttezza (resp_choice/rew).")
  if (!all(dim(y) == c(S, Tt))) stop("Dimensioni inattese per resp_choice.")
  if (!all(dim(r) == c(S, Tt))) stop("Dimensioni inattese per rew.")

  # regola dichiarata e target-by-rule
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

to_long_df <- function(obj) {
  if (is.data.frame(obj)) {
    df <- obj
    need <- c("user_id", "trial", "choice", "correct")
    if (!all(need %in% names(df)))
      stop(
        "Data frame privo di colonne minime: ",
        paste(setdiff(need, names(df)), collapse = ", ")
      )
    if (!"group" %in% names(df)) df$group <- NA
    df %>% arrange(user_id, trial)
  } else if (is.list(obj)) {
    list_to_long(obj)
  } else {
    stop("Oggetto RDS non riconosciuto (né data.frame né lista).")
  }
}

# ------------------------------
# 4) Ricostruzione regola + scoring WCST
# ------------------------------

# Ricostruisce la regola attiva dai trial corretti e propaga sui trial errati
recover_rule <- function(df) {
  df %>%
    group_by(user_id) %>%
    arrange(trial, .by_group = TRUE) %>%
    mutate(
      cand_color = as.integer(correct == 1L & choice == a_color),
      cand_shape = as.integer(correct == 1L & choice == a_form),
      cand_number = as.integer(correct == 1L & choice == a_number),
      cand_sum = cand_color + cand_shape + cand_number,
      rule_obs = case_when(
        cand_sum == 1L & cand_color == 1L ~ 1L,
        cand_sum == 1L & cand_shape == 1L ~ 2L,
        cand_sum == 1L & cand_number == 1L ~ 3L,
        TRUE ~ NA_integer_
      )
    ) %>%
    tidyr::fill(rule_obs, .direction = "down") %>%
    rename(rule_rec = rule_obs) %>%
    select(-cand_color, -cand_shape, -cand_number, -cand_sum) %>%
    ungroup()
}

.flag_shifts <- function(df, rule_col = "rule_rec") {
  df %>%
    group_by(user_id) %>%
    arrange(trial, .by_group = TRUE) %>%
    mutate(
      prev_rule = dplyr::lag(.data[[rule_col]]),
      shift = as.integer(
        !is.na(prev_rule) &
          !is.na(.data[[rule_col]]) &
          .data[[rule_col]] != prev_rule
      )
    ) %>%
    ungroup()
}

.classify_perseveration <- function(df, rule_col = "rule_rec") {
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
      # Identifica quando c'è stato uno shift di regola
      rule_switched = !is.na(prev_rule) &
        !is.na(.data[[rule_col]]) &
        (.data[[rule_col]] != prev_rule),
      # Crea un "periodo post-shift" che persiste finché non c'è un nuovo shift
      shift_period = cumsum(coalesce(rule_switched, FALSE)),
      # Target della regola precedente
      prev_target = case_when(
        prev_rule == 1L ~ a_color,
        prev_rule == 2L ~ a_form,
        prev_rule == 3L ~ a_number,
        TRUE ~ NA_integer_
      ),
      # PR: risposta perseverativa (segue target vecchio dopo shift)
      PR = as.integer(
        shift_period > 0 & !is.na(prev_target) & choice == prev_target
      ),
      # PE: errore perseverativo (PR che è anche errore)
      PE = as.integer(PR == 1L & correct == 0L),
      # NPE: errore non perseverativo (errore che non è PR)
      NPE = as.integer(correct == 0L & PR == 0L)
    ) %>%
    select(-rule_switched, -shift_period, -prev_rule, -prev_target) %>%
    ungroup()
}

.flag_fms <- function(df, rule_col = "rule_rec", run_len = 5L) {
  if (!(rule_col %in% names(df))) {
    df$FMS <- NA_integer_
    return(df)
  }
  df %>%
    group_by(user_id, .data[[rule_col]]) %>%
    arrange(trial, .by_group = TRUE) %>%
    mutate(
      correct_run = if_else(correct == 1L, 1L, 0L),
      correct_run = ave(
        correct_run,
        user_id,
        .data[[rule_col]],
        FUN = function(x) {
          r <- integer(length(x))
          run_id <- cumsum(c(1L, diff(x != 1L)))
          tapply(seq_along(x), run_id, function(idx) {
            if (x[idx[1]] == 1L) r[idx] <<- seq_along(idx) else r[idx] <<- 0L
          })
          r
        }
      )
    ) %>%
    group_by(user_id) %>%
    mutate(
      FMS = as.integer(
        correct == 0L &
          dplyr::lag(correct_run, default = 0L) >= run_len &
          (.data[[rule_col]] == dplyr::lag(.data[[rule_col]]) |
            is.na(dplyr::lag(.data[[rule_col]])))
      )
    ) %>%
    ungroup() %>%
    select(-correct_run)
}

.categories_achieved_df <- function(
  df,
  rule_col = "rule_rec",
  criterion = 10L
) {
  if (!(rule_col %in% names(df))) {
    return(
      df %>%
        group_by(user_id) %>%
        summarise(categories_achieved = NA_integer_) %>%
        ungroup()
    )
  }
  df %>%
    group_by(user_id, .data[[rule_col]]) %>%
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

# Wrapper principale
score_wcst <- function(long_df) {
  df0 <- long_df %>%
    mutate(
      correct = as.integer(correct == 1L),
      choice = as.integer(choice),
      trial = as.integer(trial)
    )

  # Regola ricostruita dai trial corretti
  df1 <- recover_rule(df0)

  # Usa la regola ricostruita per tutti i derivati
  df2 <- df1 %>%
    .flag_shifts(rule_col = "rule_rec") %>%
    .classify_perseveration(rule_col = "rule_rec") %>%
    .flag_fms(rule_col = "rule_rec")

  # Aggregazione per soggetto
  agg_basic <- df2 %>%
    group_by(user_id, group) %>%
    summarise(
      trials = n(),
      total_errors = sum(correct == 0L, na.rm = TRUE),
      perseverative_errors = sum(PE == 1L, na.rm = TRUE),
      non_persev_errors = sum(NPE == 1L, na.rm = TRUE),
      perseverative_resp = sum(PR == 1L, na.rm = TRUE),
      fms_count = sum(FMS == 1L, na.rm = TRUE),
      .groups = "drop"
    ) %>%
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
      ),

      # Alias per compatibilità con vecchi nomi (proporzioni su 60 prove)
      prop_pers_err = pe_prop_trials,
      prop_non_pers_err = npe_prop_trials,
      prop_pers_resp = pr_prop_trials
    ) %>%
    # Mappa group 1=Patients, 2=Controls
    mutate(
      group_label = dplyr::case_when(
        group == 1 ~ "Patients",
        group == 2 ~ "Controls",
        TRUE ~ "Unknown"
      )
    ) %>%
    relocate(group_label, .after = group)

  # WSLS
  wsls <- df2 %>%
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

  # Categorie completate (criterio 10 corrette consecutive)
  cats <- .categories_achieved_df(df2, rule_col = "rule_rec", criterion = 10L)

  out <- agg_basic %>%
    left_join(wsls, by = "user_id") %>%
    left_join(cats, by = "user_id") %>%
    arrange(group, user_id)

  out
}

# ------------------------------
# 5) Esecuzione
# ------------------------------

rds_obj <- readRDS(path_rds)
peek_obj(rds_obj)
long_df <- to_long_df(rds_obj)
indices <- score_wcst(long_df)

# Salva
readr::write_csv(indices, out_csv)
message("Salvato: ", out_csv)

# Anteprima
print(dplyr::glimpse(indices))
print(utils::head(indices, 10))
