# ============================================================
# WCST — Regola fissa ogni 10 trial; PR/PE su tutto il blocco post-shift
# Tre varianti PR/PE: classic-prevTarget, StaticKey-last, StaticKey-mode
# Include confronto opzionale con "old_*" in un data frame 'tot' (se presente)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble); library(stringr)
})

# ---------- PATH (modifica se serve)
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

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------- Loader: lista Stan -> long df
list_to_long <- function(L) {
  S  <- L$S %||% L$N %||% length(L$user_id)
  Tt <- L$T %||% (if (!is.null(L$resp_choice)) dim(L$resp_choice)[2] else NULL)
  if (is.null(S) || is.null(Tt)) stop("Non deduco S/T.")
  
  user_id <- L$user_id %||% seq_len(S)
  group   <- L$group   %||% rep(NA, S)
  
  y <- L$resp_choice %||% L$choice %||% L$y
  r <- L$rew         %||% L$reward %||% L$correct %||% L$r
  if (is.null(y) || is.null(r)) stop("Mancano resp_choice/rew (scelte/correttezza).")
  stopifnot(all(dim(y) == c(S, Tt)), all(dim(r) == c(S, Tt)))
  
  rule     <- L$rule_choice %||% L$rule %||% L$category %||% L$cat
  a_color  <- L$resp_color
  a_shape  <- L$resp_shape
  a_number <- L$resp_number
  
  tibble(
    user_id = rep(user_id, each = Tt),
    group   = rep(group,   each = Tt),
    trial   = rep(seq_len(Tt), times = S),
    choice  = as.integer(as.vector(t(y))),
    correct = as.integer(as.vector(t(r))),
    rule    = if (!is.null(rule)) as.integer(as.vector(t(rule))) else NA_integer_,
    a_color  = if (!is.null(a_color))  as.integer(as.vector(t(a_color)))  else NA_integer_,
    a_shape  = if (!is.null(a_shape))  as.integer(as.vector(t(a_shape)))  else NA_integer_,
    a_number = if (!is.null(a_number)) as.integer(as.vector(t(a_number))) else NA_integer_
  )
}

to_long_df <- function(obj) if (is.list(obj)) list_to_long(obj) else obj

# ---------- Utility
max_streak_of_ones <- function(x) {
  best <- 0L; cur <- 0L
  for (v in x) { if (isTRUE(v==1L)) { cur <- cur+1L; if (cur>best) best <- cur } else cur <- 0L }
  best
}
ModeInt <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_integer_)
  as.integer(names(sort(table(x), decreasing = TRUE)[1]))
}

# ---------- Scoring blocchi da 10 con 3 varianti PR/PE
score_wcst_blocks10_multi <- function(long_df, block_len = 10L, fms_run = 5L) {
  
  df <- long_df %>%
    mutate(
      correct = as.integer(correct == 1L),
      choice  = as.integer(choice),
      trial   = as.integer(trial),
      rule    = as.integer(rule),
      seg_id  = ((trial - 1L) %/% block_len) + 1L  # 1..6 per 60 trial
    )
  
  # Table dei segmenti con regole del blocco corrente e precedente
  seg_tbl <- df %>%
    group_by(user_id, seg_id) %>%
    summarise(
      seg_rule_first = dplyr::first(rule),
      seg_rule_mode  = ModeInt(rule),
      seg_rule = dplyr::coalesce(seg_rule_first, seg_rule_mode),  # regola del blocco corrente
      .groups = "drop"
    ) %>%
    arrange(user_id, seg_id) %>%
    group_by(user_id) %>%
    mutate(prev_rule_seg = dplyr::lag(seg_rule)) %>%  # NA nel 1° blocco
    ungroup()
  
  # Per StaticKey: calcolo 2 chiavi “perseverative” derivate dal BLOCCO PRECEDENTE
  # (a) last: a_{prev_rule} all'ULTIMO trial del blocco precedente
  # (b) mode: MODA di a_{prev_rule} nel blocco precedente
  build_prev_keys <- function(d) {
    d <- d %>% arrange(trial)
    d %>%
      group_by(seg_id) %>%
      summarise(
        seg_rule = dplyr::first(rule),
        last_t   = dplyr::last(trial),
        .groups = "drop_last"
      ) %>% 
      ungroup() %>%
      mutate(prev_seg_id = seg_id - 1L) %>%
      left_join(
        d %>%
          mutate(
            a_prev_color  = a_color,
            a_prev_shape  = a_shape,
            a_prev_number = a_number
          ) %>%
          group_by(seg_id) %>%
          summarise(
            # chiave perseverativa "last": a_{seg_rule} all'ultimo trial del BLOCCO
            last_key = {
              t_last <- dplyr::last(trial)
              r_this <- dplyr::first(rule)
              if (is.na(r_this)) NA_integer_ else
                switch(as.character(r_this),
                       "1" = a_color[trial == t_last][1],
                       "2" = a_shape[trial == t_last][1],
                       "3" = a_number[trial == t_last][1],
                       NA_integer_)
            },
            # chiave perseverativa "mode": moda di a_{seg_rule} nel blocco
            mode_key = {
              r_this <- dplyr::first(rule)
              vec <- switch(as.character(r_this),
                            "1" = a_color,
                            "2" = a_shape,
                            "3" = a_number,
                            NA_integer_)
              ModeInt(vec)
            },
            .groups = "drop"
          ),
        by = "seg_id"
      ) %>%
      transmute(seg_id = seg_id + 1L,  # chiave da usare nel blocco successivo
                prev_key_last = last_key,
                prev_key_mode = mode_key)
  }
  
  prev_keys <- df %>% group_by(user_id) %>% group_modify(~build_prev_keys(.x)) %>% ungroup()
  
  df <- df %>%
    left_join(seg_tbl %>% select(user_id, seg_id, seg_rule, prev_rule_seg),
              by = c("user_id","seg_id")) %>%
    left_join(prev_keys, by = c("user_id","seg_id"))
  
  # --------- PR/PE: tre varianti sul BLOCCO k>=2
  # 1) classic-prevTarget (trial-wise)
  prev_target_trial <- dplyr::case_when(
    df$prev_rule_seg == 1L ~ df$a_color,
    df$prev_rule_seg == 2L ~ df$a_shape,
    df$prev_rule_seg == 3L ~ df$a_number,
    TRUE                   ~ NA_integer_
  )
  df <- df %>%
    mutate(
      classic_PR  = as.integer(seg_id >= 2L & !is.na(prev_target_trial) & choice == prev_target_trial),
      classic_PE  = as.integer(classic_PR == 1L & correct == 0L),
      classic_NPE = as.integer(correct == 0L & classic_PR == 0L),
      
      # 2) StaticKey-last
      stat_last_PR  = as.integer(seg_id >= 2L & !is.na(prev_key_last) & choice == prev_key_last),
      stat_last_PE  = as.integer(stat_last_PR == 1L & correct == 0L),
      stat_last_NPE = as.integer(correct == 0L & stat_last_PR == 0L),
      
      # 3) StaticKey-mode
      stat_mode_PR  = as.integer(seg_id >= 2L & !is.na(prev_key_mode) & choice == prev_key_mode),
      stat_mode_PE  = as.integer(stat_mode_PR == 1L & correct == 0L),
      stat_mode_NPE = as.integer(correct == 0L & stat_mode_PR == 0L)
    )
  
  # --------- FMS (per blocco) e WS/LS (globali)
  df <- df %>%
    group_by(user_id, seg_id) %>%
    arrange(trial, .by_group = TRUE) %>%
    mutate(
      corr_flag = as.integer(correct == 1L),
      corr_run  = {
        x <- corr_flag; r <- integer(length(x))
        rid <- cumsum(c(1L, diff(x != 1L)))
        tapply(seq_along(x), rid, function(idx) {
          if (x[idx[1]] == 1L) r[idx] <<- seq_along(idx) else r[idx] <<- 0L
        }); r
      },
      FMS = as.integer(correct == 0L & dplyr::lag(corr_run, default = 0L) >= fms_run)
    ) %>%
    ungroup() %>%
    select(-corr_flag, -corr_run)
  
  wsls <- df %>%
    group_by(user_id) %>%
    arrange(trial, .by_group = TRUE) %>%
    mutate(
      stay      = as.integer(choice == dplyr::lag(choice)),
      win_stay  = as.integer(dplyr::lag(correct) == 1L & stay == 1L),
      lose_shift= as.integer(dplyr::lag(correct) == 0L & stay == 0L)
    ) %>%
    summarise(
      n_ws        = sum(win_stay, na.rm = TRUE),
      n_ws_denom  = sum(dplyr::lag(correct) == 1L, na.rm = TRUE),
      p_win_stay  = ifelse(n_ws_denom > 0, n_ws / n_ws_denom, NA_real_),
      n_ls        = sum(lose_shift, na.rm = TRUE),
      n_ls_denom  = sum(dplyr::lag(correct) == 0L, na.rm = TRUE),
      p_lose_shift= ifelse(n_ls_denom > 0, n_ls / n_ls_denom, NA_real_)
    ) %>% ungroup()
  
  cats <- df %>%
    group_by(user_id, seg_id) %>%
    summarise(max_streak = max_streak_of_ones(correct), .groups = "drop_last") %>%
    summarise(categories_achieved = sum(max_streak >= 10L, na.rm = TRUE), .groups = "drop")
  
  # --------- Aggregazione per soggetto — tre varianti
  agg_one <- function(prefix_PR, prefix_PE, prefix_NPE) {
    df %>%
      group_by(user_id, group) %>%
      summarise(
        trials               = n(),
        total_errors         = sum(correct == 0L, na.rm = TRUE),
        perseverative_errors = sum(.data[[paste0(prefix_PE)] ] == 1L, na.rm = TRUE),
        non_persev_errors    = sum(.data[[paste0(prefix_NPE)]] == 1L, na.rm = TRUE),
        perseverative_resp   = sum(.data[[paste0(prefix_PR)] ] == 1L, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        pe_prop_trials  = ifelse(trials > 0, perseverative_errors / trials, NA_real_),
        npe_prop_trials = ifelse(trials > 0, non_persev_errors    / trials, NA_real_),
        pr_prop_trials  = ifelse(trials > 0, perseverative_resp   / trials, NA_real_),
        te_prop_trials  = ifelse(trials > 0, total_errors         / trials, NA_real_),
        pe_prop_errors  = ifelse(total_errors > 0, perseverative_errors / total_errors, NA_real_),
        npe_prop_errors = ifelse(total_errors > 0, non_persev_errors    / total_errors, NA_real_)
      )
  }
  
  out_classic   <- agg_one("classic_PR","classic_PE","classic_NPE") %>%
    rename_with(~paste0("classic_", .), -c(user_id, group))
  out_stat_last <- agg_one("stat_last_PR","stat_last_PE","stat_last_NPE") %>%
    rename_with(~paste0("stat_last_", .), -c(user_id, group))
  out_stat_mode <- agg_one("stat_mode_PR","stat_mode_PE","stat_mode_NPE") %>%
    rename_with(~paste0("stat_mode_", .), -c(user_id, group))
  
  out <- out_classic %>%
    inner_join(out_stat_last, by = c("user_id","group")) %>%
    inner_join(out_stat_mode, by = c("user_id","group")) %>%
    mutate(
      group_label = dplyr::case_when(group == 1 ~ "Patients",
                                     group == 2 ~ "Controls",
                                     TRUE ~ "Unknown"),
      gr = dplyr::case_when(group == 1 ~ "an",
                            group == 2 ~ "hc",
                            TRUE ~ "na")
    ) %>%
    relocate(group_label, gr, .after = group) %>%
    left_join(wsls, by = "user_id") %>%
    left_join(cats, by = "user_id") %>%
    arrange(group, user_id)
  
  out
}

# ---------- Run
rds_obj <- readRDS(path_rds)
long_df <- to_long_df(rds_obj)
indices_all <- score_wcst_blocks10_multi(long_df, block_len = 10L, fms_run = 5L)

readr::write_csv(indices_all, out_csv)
message("Salvato: ", out_csv)
print(dplyr::glimpse(indices_all))
print(utils::head(indices_all, 5))


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

tot <- full_join(indices_all, old, by = "user_id")

plot(tot$classic_pe_prop_errors, tot$old_prop_pers_err)




# ---------- Confronto opzionale con 'tot' se presente in ambiente:
if (exists("tot")) {
  # si assume che 'tot' contenga: user_id, old_prop_pers_err, old_prop_non_pers_err, old_prop_pers_resp
  cmp <- indices_all %>%
    select(user_id,
           classic_pe_prop_trials, classic_npe_prop_trials, classic_pr_prop_trials,
           stat_last_pe_prop_trials, stat_last_npe_prop_trials, stat_last_pr_prop_trials,
           stat_mode_pe_prop_trials, stat_mode_npe_prop_trials, stat_mode_pr_prop_trials) %>%
    inner_join(tot %>% select(user_id,
                              old_prop_pers_err   = old_prop_pers_err,
                              old_prop_non_pers_err = old_prop_non_pers_err,
                              old_prop_pers_resp  = old_prop_pers_resp),
               by = "user_id") %>%
    mutate(
      d_classic_pe = classic_pe_prop_trials - old_prop_pers_err,
      d_classic_npe= classic_npe_prop_trials - old_prop_non_pers_err,
      d_classic_pr = classic_pr_prop_trials - old_prop_pers_resp,
      
      d_last_pe = stat_last_pe_prop_trials - old_prop_pers_err,
      d_last_npe= stat_last_npe_prop_trials - old_prop_non_pers_err,
      d_last_pr = stat_last_pr_prop_trials - old_prop_pers_resp,
      
      d_mode_pe = stat_mode_pe_prop_trials - old_prop_pers_err,
      d_mode_npe= stat_mode_npe_prop_trials - old_prop_non_pers_err,
      d_mode_pr = stat_mode_pr_prop_trials - old_prop_pers_resp
    )
  
  summ <- function(x) c(N = sum(!is.na(x)),
                        med = stats::median(x, na.rm = TRUE),
                        mean = mean(x, na.rm = TRUE),
                        mad = stats::mad(x, na.rm = TRUE),
                        max_abs = max(abs(x), na.rm = TRUE))
  
  message("\n--- Confronto con OLD ---")
  message("classic  ΔPE/ΔNPE/ΔPR:"); print(rbind(PE=summ(cmp$d_classic_pe), NPE=summ(cmp$d_classic_npe), PR=summ(cmp$d_classic_pr)))
  message("stat_last ΔPE/ΔNPE/ΔPR:"); print(rbind(PE=summ(cmp$d_last_pe),    NPE=summ(cmp$d_last_npe),    PR=summ(cmp$d_last_pr)))
  message("stat_mode ΔPE/ΔNPE/ΔPR:"); print(rbind(PE=summ(cmp$d_mode_pe),    NPE=summ(cmp$d_mode_npe),    PR=summ(cmp$d_mode_pr)))
  
  # suggerimento: quale variante è più vicina?
  scores <- tibble(
    variant = c("classic","stat_last","stat_mode"),
    pe_med_abs  = c(median(abs(cmp$d_classic_pe), na.rm = TRUE),
                    median(abs(cmp$d_last_pe),    na.rm = TRUE),
                    median(abs(cmp$d_mode_pe),    na.rm = TRUE)),
    pr_med_abs  = c(median(abs(cmp$d_classic_pr), na.rm = TRUE),
                    median(abs(cmp$d_last_pr),    na.rm = TRUE),
                    median(abs(cmp$d_mode_pr),    na.rm = TRUE))
  ) %>% mutate(total = pe_med_abs + pr_med_abs) %>% arrange(total)
  message("\nVariante più vicina (mediana ΔPE+ΔPR più bassa):")
  print(scores)
}






old2 <- readRDS(
  here::here(
    "data", "processed", "interim", "wcst_behavioral_stats.RDS"
  )
) |> 
  dplyr::filter(group != "ri") |> 
  dplyr::rename(user_id = subj_name)


dd <- left_join(old, old2, by = "user_id")

glimpse(dd)



