# ============================================================
# WCST — Model comparison (PSIS-LOO), parameter extraction (best model),
# and classification AUC (multivariate & univariate)
# ============================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  here,
  tidyverse,
  cmdstanr,
  posterior,
  loo,
  purrr,
  matrixStats,
  pROC,
  glmnet
)

set.seed(20251002)

# ------------------------------------------------------------
# 0) Paths & data
# ------------------------------------------------------------
stan_dir <- here::here(
  "scripts",
  "02_wcst_gen_params_hmm",
  "stan"
)
stan_data <- readRDS(here::here("data", "processed", "wcst_stan_list.RDS"))

req_core <- c("rw_dual_sticky_lapse.stan", "rw_single_alpha_dic.stan")
req_bishara <- c(
  "bishara_single_subject_with_lapse.stan",
  "bishara_single_subject_patched_modern.stan",
  "bishara_single_subject.stan"
)
bishara_file <- req_bishara[which(file.exists(file.path(
  stan_dir,
  req_bishara
)))[1]]
if (is.na(bishara_file))
  stop("File Stan per Bishara non trovato nella cartella 'stan'.")

req_extra <- c("rw_feature_based.stan", "wcst_rule_inference_hmm.stan")
all_req <- c(req_core, bishara_file, req_extra)
if (!all(file.exists(file.path(stan_dir, all_req)))) {
  missing <- all_req[!file.exists(file.path(stan_dir, all_req))]
  stop("Mancano file .stan: ", paste(missing, collapse = ", "))
}

# ------------------------------------------------------------
# 1) Data prep helpers (per modello)
# ------------------------------------------------------------
# RW classici: resp_* indicatori 0/1
prep_subject_data_rw <- function(i, use_kappa = 1L, use_lapse = 1L) {
  Ti <- stan_data$Tsubj[i]
  list(
    T = as.integer(Ti),
    choice = as.integer(stan_data$resp_choice[i, 1:Ti]),
    rew = as.integer(stan_data$rew[i, 1:Ti] == 1L),
    los = as.integer(stan_data$los[i, 1:Ti] == -1L),
    resp_color = as.integer(stan_data$resp_color[i, 1:Ti] == 1L),
    resp_shape = as.integer(stan_data$resp_shape[i, 1:Ti] == 1L),
    resp_number = as.integer(stan_data$resp_number[i, 1:Ti] == 1L),
    use_kappa = as.integer(use_kappa),
    use_lapse = as.integer(use_lapse)
  )
}

# Bishara: resp_* = INDICI 1..4; rew=1 corretto, altrimenti !=1
prep_subject_data_bishara <- function(i) {
  Ti <- stan_data$Tsubj[i]
  list(
    T = as.integer(Ti),
    resp_choice = as.integer(stan_data$resp_choice[i, 1:Ti]),
    resp_color = as.integer(stan_data$resp_color[i, 1:Ti]),
    resp_shape = as.integer(stan_data$resp_shape[i, 1:Ti]),
    resp_number = as.integer(stan_data$resp_number[i, 1:Ti]),
    rew = as.integer(stan_data$rew[i, 1:Ti])
  )
}

# RW feature-based: resp_* = INDICI 1..4; rew/los binari; sticky/lapse opzionali
prep_subject_data_feature <- function(i, use_kappa = 1L, use_lapse = 1L) {
  Ti <- stan_data$Tsubj[i]
  list(
    T = as.integer(Ti),
    choice = as.integer(stan_data$resp_choice[i, 1:Ti]),
    rew = as.integer(stan_data$rew[i, 1:Ti] == 1L),
    los = as.integer(stan_data$los[i, 1:Ti] == -1L),
    resp_color = as.integer(stan_data$resp_color[i, 1:Ti]),
    resp_shape = as.integer(stan_data$resp_shape[i, 1:Ti]),
    resp_number = as.integer(stan_data$resp_number[i, 1:Ti]),
    use_kappa = as.integer(use_kappa),
    use_lapse = as.integer(use_lapse)
  )
}

# HMM rule-inference: resp_* = INDICI 1..4; correct=1/0
prep_subject_data_hmm <- function(i) {
  Ti <- stan_data$Tsubj[i]
  list(
    T = as.integer(Ti),
    resp_choice = as.integer(stan_data$resp_choice[i, 1:Ti]),
    resp_color = as.integer(stan_data$resp_color[i, 1:Ti]),
    resp_shape = as.integer(stan_data$resp_shape[i, 1:Ti]),
    resp_number = as.integer(stan_data$resp_number[i, 1:Ti]),
    correct = as.integer(stan_data$rew[i, 1:Ti] == 1L)
  )
}

# ------------------------------------------------------------
# 2) Registry modelli + compilazione
# ------------------------------------------------------------
registry <- tibble::tribble(
  ~model_id,
  ~stan_file,
  ~prep_fun,
  ~flags,
  "full",
  "rw_dual_sticky_lapse.stan",
  prep_subject_data_rw,
  list(use_kappa = 1L, use_lapse = 1L),
  "no_kappa",
  "rw_dual_sticky_lapse.stan",
  prep_subject_data_rw,
  list(use_kappa = 0L, use_lapse = 1L),
  "no_lapse",
  "rw_dual_sticky_lapse.stan",
  prep_subject_data_rw,
  list(use_kappa = 1L, use_lapse = 0L),
  "basic_rw",
  "rw_dual_sticky_lapse.stan",
  prep_subject_data_rw,
  list(use_kappa = 0L, use_lapse = 0L),
  "rw_single_alpha_dic",
  "rw_single_alpha_dic.stan",
  prep_subject_data_rw,
  list(),
  "bishara",
  bishara_file,
  prep_subject_data_bishara,
  list(),
  "rw_feature",
  "rw_feature_based.stan",
  prep_subject_data_feature,
  list(use_kappa = 1L, use_lapse = 1L),
  "hmm_rule",
  "wcst_rule_inference_hmm.stan",
  prep_subject_data_hmm,
  list()
) %>%
  mutate(
    prep_fun = as.list(prep_fun),
    stan_path = file.path(stan_dir, stan_file),
    mod_obj = purrr::map(
      stan_path,
      ~ cmdstan_model(.x, force_recompile = FALSE)
    )
  )

# ------------------------------------------------------------
# 3) Helpers: estrazione log_lik pointwise & fit per soggetto
# ------------------------------------------------------------
extract_pointwise <- function(fit) {
  # preferisci log_lik_t; fallback a log_lik[T]
  mat <- try(fit$draws("log_lik_t", format = "draws_matrix"), silent = TRUE)
  if (!inherits(mat, "try-error") && ncol(mat) >= 1) return(mat)
  mat <- try(fit$draws("log_lik", format = "draws_matrix"), silent = TRUE)
  if (!inherits(mat, "try-error") && ncol(mat) >= 1) return(mat)
  dm <- fit$draws(format = "draws_matrix")
  all_names <- colnames(dm)
  cand <- unique(sub(
    "\\[.*$",
    "",
    all_names[grepl("^log_lik(_t)?\\[", all_names)]
  ))
  if (length(cand) >= 1) {
    mat <- try(fit$draws(cand[1], format = "draws_matrix"), silent = TRUE)
    if (!inherits(mat, "try-error") && ncol(mat) >= 1) return(mat)
  }
  message("Variabili nel fit:\n", paste(all_names, collapse = ", "))
  stop(
    "Il modello non emette log_lik per-trial (serve log_lik_t[T] o log_lik[T])."
  )
}

fit_subject_one <- function(
  row,
  i,
  chains = 2,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.9,
  max_treedepth = 12,
  keep_fit = FALSE
) {
  data_i <- do.call(row$prep_fun[[1]], c(list(i = i), row$flags[[1]]))
  fit <- row$mod_obj[[1]]$sample(
    data = data_i,
    seed = 20251002 + i,
    chains = chains,
    parallel_chains = min(chains, max(1, parallel::detectCores() - 1)),
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = 0,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
  ll_mat <- extract_pointwise(fit)
  out <- list(ll_mat = ll_mat, T_i = ncol(ll_mat))
  if (keep_fit) out$fit <- fit
  out
}

# ------------------------------------------------------------
# 4) Runner per modello: LOO per-soggetto + totale
#    (opzionalmente conserva i fit; lo useremo per hmm_rule)
# ------------------------------------------------------------
fit_model_all <- function(
  row,
  keep_fits = FALSE,
  chains = 2,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.9,
  max_treedepth = 12
) {
  N <- stan_data$N
  message("Model ", row$model_id)
  safe <- purrr::safely(fit_subject_one, otherwise = NULL, quiet = FALSE)
  res <- purrr::map(
    seq_len(N),
    ~ safe(
      row = row,
      i = .x,
      chains = chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      keep_fit = keep_fits
    )
  )
  fail <- which(purrr::map_lgl(res, ~ !is.null(.x$error)))
  if (length(fail)) {
    message("Falliti soggetti: ", paste(fail, collapse = ", "))
    message("Primo errore: ", res[[fail[1]]]$error$message)
  }
  ok <- purrr::compact(purrr::map(res, "result"))
  stopifnot(length(ok) > 0)

  # allinea numero di draws
  n_draws <- unique(purrr::map_int(ok, ~ nrow(.x$ll_mat)))
  if (length(n_draws) > 1L) {
    minD <- min(n_draws)
    ok <- purrr::map(
      ok,
      ~ {
        .x$ll_mat <- .x$ll_mat[seq_len(minD), , drop = FALSE]
        .x
      }
    )
  }

  # LOO per soggetto (ELPD additivo)
  loo_list <- purrr::map(ok, ~ loo::loo(.x$ll_mat))
  elpd_sum <- sum(purrr::map_dbl(
    loo_list,
    ~ .x$estimates["elpd_loo", "Estimate"]
  ))
  se_sum <- sqrt(sum(purrr::map_dbl(
    loo_list,
    ~ (.x$estimates["elpd_loo", "SE"])^2
  )))

  # LOO totale concatenando i trials (per loo_compare)
  big_ll <- do.call(cbind, purrr::map(ok, "ll_mat"))
  loo_total <- loo::loo(big_ll)

  # Pareto-k summary
  pareto_k_all <- unlist(purrr::map(loo_list, ~ .x$diagnostics$pareto_k))
  pk_summary <- list(
    frac_k_gt_0.7 = mean(pareto_k_all > 0.7),
    frac_k_gt_1 = mean(pareto_k_all > 1)
  )

  out <- list(
    model_id = row$model_id,
    loo = loo_total,
    loo_list = loo_list,
    elpd_loo_total = elpd_sum,
    se_elpd_total = se_sum,
    pareto_k_summary = pk_summary
  )
  if (keep_fits) out$fits <- purrr::map(ok, "fit")
  out
}

# ------------------------------------------------------------
# 5) Esecuzione: tutti i modelli (conserva i fit solo per hmm_rule)
# ------------------------------------------------------------
fits_loo <- purrr::map(seq_len(nrow(registry)), function(ix) {
  row <- dplyr::slice(registry, ix)
  keep <- (row$model_id == "hmm_rule")
  fit_model_all(
    row,
    keep_fits = keep,
    chains = 2,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.9,
    max_treedepth = 12
  )
})

# Tabella LOO e ranking
loo_tbl <- tibble(
  model_id = purrr::map_chr(fits_loo, "model_id"),
  elpd_loo = purrr::map_dbl(fits_loo, "elpd_loo_total"),
  se_elpd = purrr::map_dbl(fits_loo, "se_elpd_total"),
  frac_k_gt_0.7 = purrr::map_dbl(fits_loo, ~ .x$pareto_k_summary$frac_k_gt_0.7),
  frac_k_gt_1 = purrr::map_dbl(fits_loo, ~ .x$pareto_k_summary$frac_k_gt_1)
) %>%
  arrange(desc(elpd_loo))
print(loo_tbl, n = Inf)
# model_id            elpd_loo se_elpd frac_k_gt_0.7 frac_k_gt_1
# <chr>                  <dbl>   <dbl>         <dbl>       <dbl>
# 1 hmm_rule              -3869.   49.5         0           0
# 2 bishara               -4511.   69.2         0.0428      0.0311
# 3 rw_feature            -5961.   45.0         0.0170      0.0167
# 4 no_kappa              -7419.    9.04        0.0254      0.0212
# 5 basic_rw              -7430.    9.76        0.0267      0.0205
# 6 rw_single_alpha_dic   -7434.    9.57        0.0167      0.0167
# 7 full                  -7441.   12.8         0.0167      0.0167
# 8 no_lapse              -7453.   13.8         0.0167      0.0167

# Confronto ΔELPD ± SE
loo_list <- setNames(
  purrr::map(fits_loo, "loo"),
  purrr::map_chr(fits_loo, "model_id")
)
print(loo::loo_compare(loo_list))
#                     elpd_diff se_diff
# hmm_rule                0.0       0.0
# bishara              -641.7      33.7
# rw_feature          -2092.4      33.6
# no_kappa            -3549.8      50.0
# basic_rw            -3560.8      50.1
# rw_single_alpha_dic -3565.0      50.3
# full                -3572.3      51.7
# no_lapse            -3584.3      52.0

# Best model
best_model_id <- loo_tbl$model_id[1]
message("Miglior modello (ELPD): ", best_model_id)
# Miglior modello (ELPD): hmm_rule

# ------------------------------------------------------------
# 6) Estrazione parametri per soggetto (best model = hmm_rule atteso)
#    Usa i fit conservati; se mancassero, rifitta solo hmm_rule.
# ------------------------------------------------------------

get_row <- function(mid) dplyr::slice(registry, which(registry$model_id == mid))

extract_params_hmm <- function(fit) {
  # Variabili chiave: h, d, lapse, eta
  summ <- fit$summary(variables = c("h", "d", "lapse", "eta"))
  summ %>%
    dplyr::select(variable, mean, median, q5, q95, rhat, ess_bulk, ess_tail)
}

if (best_model_id != "hmm_rule") {
  warning(
    "Il modello migliore non è 'hmm_rule' (è ",
    best_model_id,
    "). Estrarrò comunque i parametri di hmm_rule rifittando."
  )
}

# Recupera fit list (se disponibile) o rifitta
hmm_ix <- which(purrr::map_chr(fits_loo, "model_id") == "hmm_rule")
if (length(hmm_ix) == 1 && !is.null(fits_loo[[hmm_ix]]$fits)) {
  fits_hmm <- fits_loo[[hmm_ix]]$fits
} else {
  row_hmm <- get_row("hmm_rule")
  fits_hmm <- purrr::map(seq_len(stan_data$N), function(i) {
    fit_subject_one(
      row_hmm,
      i,
      chains = 4,
      iter_warmup = 1000,
      iter_sampling = 1000,
      adapt_delta = 0.9,
      max_treedepth = 12,
      keep_fit = TRUE
    )$fit
  })
}

pars_hmm <- purrr::imap_dfr(fits_hmm, function(ft, sid) {
  extract_params_hmm(ft) %>%
    tidyr::pivot_wider(
      names_from = variable,
      values_from = c(mean, median, q5, q95, rhat, ess_bulk, ess_tail),
      names_sep = "_"
    ) %>%
    dplyr::mutate(subject = sid)
})

# ---------------------------
# Patch etichette gruppi: 42 -> Patients, 46 -> Controls
# ---------------------------

# livelli numerici effettivi presenti in stan_data$group (es. 1 e 2, oppure 0 e 1)
lvl_found <- sort(unique(stan_data$group))
stopifnot(length(lvl_found) == 2L)

# conteggi per ciascun livello nell'ordine di lvl_found
grp_counts <- as.integer(table(factor(stan_data$group, levels = lvl_found)))

# Verifica che i conteggi siano proprio 42 e 46 (in qualunque ordine)
if (!setequal(grp_counts, c(42L, 46L))) {
  stop(
    "I conteggi dei gruppi non sono 42 e 46. Ho trovato: ",
    paste(sprintf("%s=%d", lvl_found, grp_counts), collapse = ", "),
    ". Controlla stan_data$group."
  )
}

# Mappa etichette: 42 -> "Patients", 46 -> "Controls"
labs <- if (grp_counts[1] == 42L) c("Patients", "Controls") else
  c("Controls", "Patients")

# Aggiungi colonna group (codici grezzi) e group_label (etichette coerenti ai conteggi)
pars_hmm <- pars_hmm %>%
  dplyr::mutate(
    group = stan_data$group[subject],
    group_label = factor(group, levels = lvl_found, labels = labs)
  ) %>%
  # Ordine livelli coerente con pipeline ROC/CV: "Controls" prima, poi "Patients"
  dplyr::mutate(
    group_label = forcats::fct_relevel(group_label, "Controls", "Patients")
  )

# Controlli diagnostici MCMC (facoltativi)
max_rhat <- pars_hmm %>%
  dplyr::select(dplyr::matches("rhat_")) %>%
  as.matrix() %>%
  max(na.rm = TRUE)
message("Max R-hat (hmm_rule): ", round(max_rhat, 3))

# Sanity checks etichette/conteggi
print(table(stan_data$group)) # conteggi per i codici grezzi
print(table(pars_hmm$group_label)) # atteso: Controls 46, Patients 42

# ------------------------------------------------------------
# 6-bis) Trasformazioni su scala interpretabile + check colonne
#       (DA INSERIRE SUBITO DOPO la costruzione di `pars_hmm`)
# ------------------------------------------------------------
eps <- 1e-6

# Costruisci il feature set trasformato
pars_hmm_x <- pars_hmm %>%
  dplyr::transmute(
    subject,
    group,
    group_label, # già correttamente etichettato e riordinato
    # trasformazioni come nel classificatore
    logit_h = qlogis(pmin(pmax(mean_h, eps), 1 - eps)),
    log_d = log(pmax(mean_d, eps)),
    logit_lapse = qlogis(pmin(pmax(mean_lapse, eps), 1 - eps)),
    logit_eta = qlogis(pmin(pmax(mean_eta, eps), 1 - eps))
  )

# Check livelli fattoriali e conteggi attesi
stopifnot(identical(levels(pars_hmm_x$group_label), c("Controls", "Patients")))
print(table(pars_hmm_x$group_label)) # atteso: Controls 46, Patients 42

# Controllo: le colonne attese esistono davvero in pars_hmm_x
stopifnot(all(
  c("group_label", "logit_h", "log_d", "logit_lapse", "logit_eta") %in%
    names(pars_hmm_x)
))

# Save pars_hmm_x for later use (Bayesian assurance analysis)
saveRDS(
  pars_hmm_x,
  here::here(
    "data",
    "processed",
    "params",
    "pars_hmm_x.RDS"
  )
)

# --- 2) Long format + etichette leggibili
plot_df <- pars_hmm_x %>%
  pivot_longer(
    cols = c(logit_h, log_d, logit_lapse, logit_eta),
    names_to = "param",
    values_to = "value"
  ) %>%
  mutate(
    param = forcats::fct_recode(
      param,
      "h (hazard) — logit" = "logit_h",
      "d (decision) — log" = "log_d",
      "lapse — logit" = "logit_lapse",
      "eta (fb noise) — logit" = "logit_eta"
    )
  )

# --- 3) Grafico: boxplot + jitter (punti)
p <- ggplot(plot_df, aes(x = group_label, y = value, fill = group_label)) +
  geom_boxplot(width = 0.55, alpha = 0.55, outlier.shape = NA) +
  geom_jitter(
    aes(color = group_label),
    width = 0.12,
    height = 0,
    size = 1.8,
    alpha = 0.70
  ) +
  facet_wrap(~param, scales = "free_y", ncol = 2) +
  labs(
    x = NULL,
    y = "Parametro (scala trasformata)",
    title = "Parametri HMM-rule (WCST) per gruppo",
    subtitle = "Boxplot + singoli partecipanti"
  ) +
  guides(fill = "none", color = "none") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

print(p)

# --- 4) Salvataggio
out_dir <- here::here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
ggsave(
  file.path(out_dir, "hmm_rule_params_boxpoints.png"),
  p,
  width = 9,
  height = 7,
  dpi = 300
)
ggsave(
  file.path(out_dir, "hmm_rule_params_boxpoints.pdf"),
  p,
  width = 9,
  height = 7
)


# --- (Opzionale) Variante con punti "impaccati" (beeswarm) ---
# pacman::p_load(ggbeeswarm)
# p2 <- ggplot(plot_df, aes(x = group_label, y = value, fill = group_label)) +
#   geom_boxplot(width = 0.55, alpha = 0.55, outlier.shape = NA) +
#   ggbeeswarm::geom_quasirandom(aes(color = group_label),
#                                width = 0.25, size = 1.6, alpha = 0.75) +
#   facet_wrap(~ param, scales = "free_y", ncol = 2) +
#   labs(x = NULL, y = "Parametro (scala trasformata)",
#        title = "Parametri HMM-rule (WCST) per gruppo — Beeswarm") +
#   guides(fill = "none", color = "none") +
#   theme_minimal(base_size = 12) +
#   theme(panel.grid.minor = element_blank(),
#         strip.text = element_text(face = "bold"),
#         plot.title = element_text(face = "bold"))
# ggsave(file.path(out_dir, "hmm_rule_params_boxpoints_beeswarm.png"), p2, width = 9, height = 7, dpi = 300)

# ------------------------------------------------------------
# 7) Feature set trasformato per classificazione
#    (orientamento coerente: Patients = classe positiva)
# ------------------------------------------------------------
eps <- 1e-6
pars_hmm_x <- pars_hmm %>%
  transmute(
    subject,
    group,
    group_label,
    logit_h = qlogis(pmin(pmax(mean_h, eps), 1 - eps)),
    log_d = log(pmax(mean_d, eps)),
    logit_lapse = qlogis(pmin(pmax(mean_lapse, eps), 1 - eps)),
    logit_eta = qlogis(pmin(pmax(mean_eta, eps), 1 - eps))
  )

# ------------------------------------------------------------
# 8) AUC multivariata con CV stratificata (ridge/enet) + CI DeLong
# ------------------------------------------------------------
make_stratified_folds <- function(y, K = 10, seed = 1L) {
  set.seed(seed)
  y <- factor(y)
  folds <- rep(NA_integer_, length(y))
  for (lvl in levels(y)) {
    idx <- which(y == lvl)
    idx <- sample(idx)
    parts <- split(idx, rep(1:K, length.out = length(idx)))
    for (k in seq_along(parts)) folds[parts[[k]]] <- k
  }
  folds
}

cv_auc_glmnet <- function(X, y_fac, alpha = 0, K = 10, seed = 1L) {
  stopifnot(
    is.factor(y_fac) && identical(levels(y_fac), c("Controls", "Patients"))
  )
  folds <- make_stratified_folds(y_fac, K = K, seed = seed)
  oof_pred <- rep(NA_real_, nrow(X))
  for (k in 1:K) {
    tr <- which(folds != k)
    te <- which(folds == k)
    Xtr <- X[tr, , drop = FALSE]
    ytr <- (y_fac[tr] == "Patients")
    Xte <- X[te, , drop = FALSE]
    cvfit <- glmnet::cv.glmnet(
      x = Xtr,
      y = ytr,
      family = "binomial",
      alpha = alpha,
      type.measure = "auc",
      nfolds = 5,
      standardize = TRUE
    )
    oof_pred[te] <- as.numeric(predict(
      cvfit,
      newx = Xte,
      s = "lambda.min",
      type = "response"
    ))
  }
  roc_obj <- pROC::roc(
    response = y_fac,
    predictor = oof_pred,
    levels = c("Controls", "Patients"),
    direction = "auto",
    quiet = TRUE
  )
  list(
    auc = as.numeric(pROC::auc(roc_obj)),
    ci_auc = as.numeric(pROC::ci.auc(roc_obj, method = "delong")),
    oof_pred = oof_pred,
    folds = folds,
    roc = roc_obj
  )
}

X_all <- as.matrix(
  pars_hmm_x %>% dplyr::select(logit_h, log_d, logit_lapse, logit_eta)
)
y_fac <- factor(pars_hmm_x$group_label, levels = c("Controls", "Patients"))

ridge_res <- cv_auc_glmnet(X_all, y_fac, alpha = 0, K = 10, seed = 20251002)
enet_res <- cv_auc_glmnet(X_all, y_fac, alpha = 0.5, K = 10, seed = 20251002)

cat(sprintf(
  "Ridge AUC (10-fold CV):  %.3f  [%.3f, %.3f]\n",
  ridge_res$auc,
  ridge_res$ci_auc[1],
  ridge_res$ci_auc[3]
))
# Ridge AUC (10-fold CV):  0.661  [0.547, 0.776]

cat(sprintf(
  "ENet  AUC (10-fold CV):  %.3f  [%.3f, %.3f]\n",
  enet_res$auc,
  enet_res$ci_auc[1],
  enet_res$ci_auc[3]
))
# ENet  AUC (10-fold CV):  0.617  [0.498, 0.737]

# CV ripetuta + grid alpha (robustezza)
alpha_grid <- c(0, 0.25, 0.5, 0.75, 1)
R <- 50
set.seed(20251002)
res_rep <- purrr::map(alpha_grid, function(a) {
  aucs <- replicate(
    R,
    cv_auc_glmnet(
      X_all,
      y_fac,
      alpha = a,
      K = 10,
      seed = sample.int(1e6, 1)
    )$auc
  )
  tibble(
    alpha = a,
    auc_mean = mean(aucs),
    auc_sd = sd(aucs),
    auc_ci_low = quantile(aucs, .025),
    auc_ci_high = quantile(aucs, .975)
  )
}) %>%
  bind_rows() %>%
  arrange(desc(auc_mean))
print(res_rep)
# alpha auc_mean auc_sd auc_ci_low auc_ci_high
# <dbl>    <dbl>  <dbl>      <dbl>       <dbl>
# 1  0.25    0.655 0.0264      0.591       0.697
# 2  0.5     0.655 0.0263      0.612       0.707
# 3  1       0.649 0.0258      0.590       0.688
# 4  0.75    0.648 0.0307      0.594       0.695
# 5  0       0.636 0.0314      0.571       0.680
# ------------------------------------------------------------
# 9) AUC univariata per parametro + Hedges g + rank-biserial
# ------------------------------------------------------------
hedges_g <- function(x_pat, x_ctrl) {
  n1 <- length(x_pat)
  n2 <- length(x_ctrl)
  m1 <- mean(x_pat)
  m2 <- mean(x_ctrl)
  s1 <- stats::var(x_pat)
  s2 <- stats::var(x_ctrl)
  s_p <- sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
  g_unbiased <- (m1 - m2) / s_p
  J <- 1 - 3 / (4 * (n1 + n2) - 9)
  g_unbiased * J
}

params <- c("logit_h", "log_d", "logit_lapse", "logit_eta")
uni_tbl <- purrr::map_dfr(params, function(par) {
  df <- pars_hmm_x %>%
    dplyr::select(group_label, !!sym(par)) %>%
    dplyr::rename(value = !!sym(par))
  roc_obj <- pROC::roc(
    response = factor(df$group_label, levels = c("Controls", "Patients")),
    predictor = df$value,
    levels = c("Controls", "Patients"),
    direction = "auto",
    quiet = TRUE
  )
  auc_val <- as.numeric(pROC::auc(roc_obj))
  ci_val <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  x_pat <- df$value[df$group_label == "Patients"]
  x_ctrl <- df$value[df$group_label == "Controls"]
  g_val <- hedges_g(x_pat, x_ctrl)
  rb <- 2 * auc_val - 1
  tibble(
    parameter = par,
    auc = auc_val,
    auc_ci_low = ci_val[1],
    auc_ci_high = ci_val[3],
    hedges_g = g_val,
    rank_biserial = rb
  )
})
print(uni_tbl)
# parameter     auc auc_ci_low auc_ci_high hedges_g rank_biserial
# <chr>       <dbl>      <dbl>       <dbl>    <dbl>         <dbl>
# 1 logit_h     0.669      0.554       0.784  -0.550        0.339
# 2 log_d       0.498      0.375       0.622  -0.0429      -0.00311
# 3 logit_lapse 0.503      0.379       0.627   0.0332       0.00621
# 4 logit_eta   0.695      0.582       0.807  -0.695        0.389

# ------------------------------------------------------------
# 10) Salvataggi
# ------------------------------------------------------------
out_dir <- here::here("src", "wcst", "documentation", "outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

readr::write_csv(
  loo_tbl,
  file.path(out_dir, "model_comparison_loo_summary.csv")
)
saveRDS(
  list(
    registry = registry,
    fits_loo = fits_loo,
    loo_compare = loo::loo_compare(loo_list)
  ),
  file.path(out_dir, "model_comparison_loo.rds")
)

readr::write_csv(pars_hmm, file.path(out_dir, "hmm_rule_params_raw.csv"))
readr::write_csv(pars_hmm_x, file.path(out_dir, "hmm_rule_params_features.csv"))

readr::write_csv(res_rep, file.path(out_dir, "auc_cv_repeated_alpha_grid.csv"))
readr::write_csv(uni_tbl, file.path(out_dir, "auc_univariate_params.csv"))

# curva ROC (ridge)
png(
  file.path(out_dir, "roc_ridge_oof.png"),
  width = 900,
  height = 700,
  res = 120
)
plot(
  ridge_res$roc,
  print.auc = TRUE,
  main = "ROC — hmm_rule params (Ridge, out-of-fold)"
)
dev.off()

#' Implicazioni per i deficit nelle pazienti
#' Basso h nelle pazienti è un marcatore diretto di minor flessibilità/set-shifting:
#' presuppongono che la regola non cambi e, quindi, aggiornano più lentamente
#' dopo gli switch — un profilo classico di perseverazione.
#' Basso eta nelle pazienti indica che percepiscono il feedback come molto
#' affidabile. Nel WCST il feedback è effettivamente deterministico, quindi in
#' sé questo non è “patologico”. Tuttavia, combinato con h basso significa che:
#' finché la regola sembra stabile, un singolo errore potrebbe essere
#' inizialmente “assorbito” dall’aspettativa di stabilità (basso h), ma più
#' errori consecutivi forniscono evidenza molto forte (basso eta ⇒ feedback “pesa”
#' tanto), portando a switch tardivi ma talvolta drastici.
#' In pratica: resistenza iniziale al cambiamento, poi aggiornamenti bruschi
#' quando l’evidenza diventa schiacciante. Questo pattern è compatibile con un
#' deficit di flessibilità più che con “rumore” decisionale o lapsi casuali (che
#' qui non differiscono).
#'
#' Risultato chiave: i parametri inferenziali legati allo switching di regola (h)
#' e all’affidabilità del feedback (eta) sono quelli che discriminano i gruppi,
#' mentre determinismo decisionale (d) e lapse no.
#' Interpretazione: il deficit non è “rumore” o incoerenza della risposta, ma una
#' credenza di stabilità (basso h) che rende lente le pazienti nell’adattarsi ai
#' cambi di regola; l’affidabilità percepita del feedback (basso eta) accentua
#' l’aggiornamento solo quando l’evidenza contraria accumulata supera la forte
#' aspettativa di stabilità.
#' Forza dell’evidenza: AUC univariate moderate per eta e h; AUC multivariata
#' ~0.67 con CV conferma che i parametri computazionali hanno potere
#' discriminante superiore agli indici comportamentali “grezzi” attesi (da
#' verificare nello stesso schema di CV per un confronto fair).

# Save params in .csv format ---------------------------------------------------

temp <- pars_hmm_x |>
  mutate(is_patient = ifelse(group_label == "Patients", 1, 0)) |>
  dplyr::select(-c("group", "group_label"))

rio::export(
  temp,
  here::here(
    "src",
    "stan_auc",
    "models_params",
    "wcst_hmm_params.csv"
  )
)

# Additional checks ------------------------------------------------------------

library(pROC)
for (par in c("logit_h", "logit_eta", "log_d", "logit_lapse")) {
  roc_obj <- roc(
    response = factor(
      pars_hmm_x$group_label,
      levels = c("Controls", "Patients")
    ),
    predictor = pars_hmm_x[[par]],
    quiet = TRUE
  )
  cat(par, "AUC =", as.numeric(auc(roc_obj)), "\n")
}


set.seed(1)
boot_diff <- function(v, g, B = 5000) {
  x <- v[g == "Patients"]
  y <- v[g == "Controls"]
  n1 <- length(x)
  n0 <- length(y)
  diffs <- replicate(B, {
    mean(sample(x, n1, TRUE)) - mean(sample(y, n0, TRUE))
  })
  c(
    diff = mean(x) - mean(y),
    lcl95 = quantile(diffs, .025),
    ucl95 = quantile(diffs, .975),
    p_dir = mean(diffs > 0)
  )
}

res_h <- boot_diff(pars_hmm_x$logit_h, pars_hmm_x$group_label)
res_eta <- boot_diff(pars_hmm_x$logit_eta, pars_hmm_x$group_label)
print(res_h)
print(res_eta)


invlogit <- function(x) 1 / (1 + exp(-x))

set.seed(1)
boot_diff_prob <- function(v_logit, g, B = 5000) {
  x <- v_logit[g == "Patients"]
  y <- v_logit[g == "Controls"]
  n1 <- length(x)
  n0 <- length(y)
  diffs <- replicate(B, {
    invlogit(mean(sample(x, n1, TRUE))) - invlogit(mean(sample(y, n0, TRUE)))
  })
  c(
    diff = invlogit(mean(x)) - invlogit(mean(y)),
    lcl95 = quantile(diffs, .025),
    ucl95 = quantile(diffs, .975),
    p_dir = mean(diffs > 0)
  )
}

res_h_p <- boot_diff_prob(pars_hmm_x$logit_h, pars_hmm_x$group_label)
res_eta_p <- boot_diff_prob(pars_hmm_x$logit_eta, pars_hmm_x$group_label)
res_h_p
res_eta_p
