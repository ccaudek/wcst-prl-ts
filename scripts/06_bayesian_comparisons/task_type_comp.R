if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  dplyr,
  tidyr,
  stringr,
  purrr,
  tibble,
  ggplot2,
  brms,
  posterior,
  pROC
)


wcst_params <- rio::import(
  here::here(
    "data",
    "processed",
    "wcst_hmm_params.csv"
  )
)
wcst_params$task <- "wcst"
wcst_params$type <- "params"

wcst_params$is_patient <- ifelse(wcst_params$group_label == "Patients", 1, 0)
wcst_params$group_label <- NULL

wcst_bi <- rio::import(here::here(
  "data",
  "processed",
  "wcst_behav_indices.csv"
))
wcst_bi$task <- "wcst"
wcst_bi$type <- "behav"

wcst_bi$is_patient <- ifelse(wcst_bi$group == "an", 1, 0)
wcst_bi$group <- NULL

prl_params <- rio::import(here::here(
  "data",
  "processed",
  "prl_params.csv"
))
prl_params$task <- "prl"
prl_params$type <- "params"

prl_params$is_patient <- ifelse(prl_params$group == "patients", 1, 0)
prl_params$group <- NULL


prl_bi <- rio::import(here::here(
  "data",
  "processed",
  "prl_behav_indices.csv"
))
prl_bi$task <- "prl"
prl_bi$type <- "behav"

prl_bi$is_patient <- ifelse(prl_bi$group == "AN", 1, 0)
prl_bi$group <- NULL


ts_params <- rio::import(here::here(
  "data",
  "processed",
  "task_switching_params.csv"
))
ts_params$task <- "ts"
ts_params$type <- "params"

ts_params$is_patient <- ifelse(ts_params$group == "AN", 1, 0)
ts_params$group <- NULL


ts_bi <- rio::import(here::here(
  "data",
  "processed",
  "task_switching_behav_indices.csv"
))
ts_bi$task <- "ts"
ts_bi$type <- "behav"


wcst_params |> head()
wcst_bi |> head()

ts_params |> head()
ts_bi |> head()

prl_params |> head()
prl_bi |> head()


# ============================================================
# brms: modello gerarchico per discriminare AN vs HC
# con compito (WCST/PRL/TS), tipo (behav/params) e feature multivariate
# ============================================================

# ---------------------------
# 0) Funzioni di utilità
# ---------------------------
to_long <- function(df, id_col, task, type, feature_cols) {
  stopifnot(all(c(id_col, "is_patient") %in% names(df)))
  df %>%
    select(all_of(c(id_col, "is_patient", feature_cols))) %>%
    rename(user_id = !!id_col) %>%
    pivot_longer(
      cols = all_of(feature_cols),
      names_to = "feature",
      values_to = "value"
    ) %>%
    mutate(task = task, type = type, .before = feature)
}

z_within_feature <- function(dat) {
  dat %>%
    group_by(feature) %>%
    mutate(z_value = as.numeric(scale(value))) %>%
    ungroup()
}

# ---------------------------
# 1) Armonizza i dataset
#    (usa solo quelli che esistono in ambiente)
# ---------------------------
have <- ls()

pieces <- list()

# --- WCST ---
if (all(c("wcst_params", "wcst_bi") %in% have)) {
  w_params_features <- c("logit_h", "log_d", "logit_lapse", "logit_eta")
  w_bi_features <- c("prop_pers_err", "prop_non_pers_err", "prop_pers_resp")

  stopifnot(all(w_params_features %in% names(wcst_params)))
  stopifnot(all(w_bi_features %in% names(wcst_bi)))

  pieces <- c(
    pieces,
    list(to_long(
      wcst_params,
      id_col = "user_id",
      task = "wcst",
      type = "params",
      feature_cols = w_params_features
    )),
    list(to_long(
      wcst_bi %>% rename(user_id = subj_name),
      id_col = "user_id",
      task = "wcst",
      type = "behav",
      feature_cols = w_bi_features
    ))
  )
}

# --- Task switching (ts_*) ---
if (all(c("ts_params", "ts_bi") %in% have)) {
  # Nota: nel tuo esempio ts_* ha task="prl" -> forzo "ts"
  ts_params2 <- ts_params %>% mutate(task = "ts")
  ts_bi2 <- ts_bi %>% mutate(task = "ts")

  # Parametri: adattali ai tuoi nomi effettivi
  ts_params_features <- setdiff(
    names(ts_params2),
    c("user_id", "is_patient", "task", "type")
  )
  ts_params_features <- ts_params_features[
    !str_detect(ts_params_features, "^is_patient$")
  ]

  # Indici comportamentali (es. repetition/switch RT): adatta ai tuoi nomi
  ts_bi_features <- setdiff(
    names(ts_bi2),
    c("user_id", "is_patient", "task", "type")
  )
  ts_bi_features <- ts_bi_features[!str_detect(ts_bi_features, "^is_patient$")]

  pieces <- c(
    pieces,
    list(to_long(
      ts_params2 %>% mutate(type = "params"),
      id_col = "user_id",
      task = "ts",
      type = "params",
      feature_cols = ts_params_features
    )),
    list(to_long(
      ts_bi2 %>% mutate(type = "behav"),
      id_col = "user_id",
      task = "ts",
      type = "behav",
      feature_cols = ts_bi_features
    ))
  )
}

# --- PRL ---
if (all(c("prl_params", "prl_bi") %in% have)) {
  p_params_features <- c("a", "alpha", "pos_alpha", "t", "v")
  # Se i tuoi nomi includono prefissi (es. a, alpha, pos_alpha, t, v), lascia così;
  # altrimenti adegua ai nomi reali delle colonne.
  p_params_features <- p_params_features[
    p_params_features %in% names(prl_params)
  ]

  p_bi_features <- c("win_stay", "lose_shift", "e_pers")
  p_bi_features <- p_bi_features[p_bi_features %in% names(prl_bi)]

  pieces <- c(
    pieces,
    list(to_long(
      prl_params,
      id_col = "user_id",
      task = "prl",
      type = "params",
      feature_cols = p_params_features
    )),
    list(to_long(
      prl_bi,
      id_col = "user_id",
      task = "prl",
      type = "behav",
      feature_cols = p_bi_features
    ))
  )
}

stopifnot(length(pieces) > 0)
dat_long <- bind_rows(pieces)

# Coerenza outcome
dat_long <- dat_long %>%
  mutate(
    is_patient = as.integer(is_patient),
    task = factor(task, levels = c("wcst", "prl", "ts")),
    type = factor(type, levels = c("behav", "params"))
  )

# Rimuovi righe non informative (NA, inf)
dat_long <- dat_long %>%
  filter(is.finite(value)) %>%
  drop_na(is_patient, value, task, type, feature, user_id)

# Z-score per feature (maggior sensibilità e comparabilità)
dat_long <- z_within_feature(dat_long)

# Controllo rapido
message("Righe totali: ", nrow(dat_long))
message("Soggetti unici: ", length(unique(dat_long$user_id)))
message("Feature uniche: ", length(unique(dat_long$feature)))
print(table(dat_long$task, dat_long$type))

# ---------------------------
# 2) Modello brms
# ---------------------------
# Formula:
# - effetti fissi per task, type e interazione
# - un singolo 'z_value' come predittore generico (perché ogni riga = 1 feature)
# - effetti casuali per feature (intercetto e slope di z_value)
# - intercetto casuale per soggetto (corregge la dipendenza intra-soggetto)
bf_form <- bf(
  is_patient ~
    0 +
      task +
      type +
      task:type +
      z_value +
      (1 + z_value | feature) +
      (1 | user_id),
  family = bernoulli()
)

priors <- c(
  prior(normal(0, 1), class = "b"), # fissati (task/type)
  prior(normal(0, 1), class = "b", coef = "z_value"),
  prior(student_t(3, 0, 2.5), class = "sd"), # deviazioni random
  prior(lkj(2), class = "cor") # correlazioni random
)

# Opzioni MCMC sensate; aumenta iter per stabilità se necessario
fit <- brm(
  formula = bf_form,
  data = dat_long,
  prior = priors,
  cores = max(2, parallel::detectCores() - 1),
  chains = 4,
  iter = 4000,
  warmup = 1500,
  seed = 20251013,
  control = list(adapt_delta = 0.97, max_treedepth = 12)
)

print(fit)
loo_fit <- loo(fit)
print(loo_fit)

# ---------------------------
# 3) Contrasti tra compiti (effetti fissi)
# ---------------------------
# Quanto (in media) un compito discrimina rispetto a un altro (a parità di tipo)?
# Esempi:
hypothesis(fit, "taskprl - taskwcst > 0") # PRL vs WCST
hypothesis(fit, "taskts - taskwcst > 0") # TS vs WCST
hypothesis(fit, "taskts - taskprl > 0") # TS vs PRL

# Interazione tipo:params vs behav (effetto medio dei parametri vs indici)
hypothesis(fit, "typeparams - typebehav > 0")

# Interazioni task:type (es. params meglio di behav nel PRL?)
# Nota: i nomi dei coefficienti seguono la codifica 0+task+type+task:type
fixef_names <- names(fixef(fit)[, "Estimate"])
print(fixef_names)

# ---------------------------
# 4) Ranking delle feature (slopes specifici)
# ---------------------------
# Estrai gli slope di z_value per ciascuna feature (quanto è utile una feature)
draws_ranef <- ranef(fit)$feature
# colonne tipiche: "(Intercept)", "z_value", e relative SD/CI
feature_slopes <- tibble(
  feature = rownames(draws_ranef[,, "z_value"]),
  slope_est = draws_ranef[,, "z_value"][, "Estimate"],
  slope_l95 = draws_ranef[,, "z_value"][, "Q2.5"],
  slope_u95 = draws_ranef[,, "z_value"][, "Q97.5"]
) %>%
  arrange(desc(slope_est))

print(head(feature_slopes, 10))

# ---------------------------
# 5) AUC posteriore per task (check predittivo)
# ---------------------------
# Calcola p_hat a livello di riga, poi sintetizza a livello soggetto per task
# (media delle predizioni sulle feature di quel soggetto e task)
epred <- posterior_epred(fit, newdata = dat_long, ndraws = 800)
# epred: draws x N_obs
p_hat <- apply(epred, 2, mean) # media sui draws

dat_pred <- dat_long %>%
  mutate(p_hat = p_hat) %>%
  group_by(user_id, is_patient, task) %>%
  summarise(p_hat_task = mean(p_hat), .groups = "drop")

# AUC per ciascun task
auc_by_task <- dat_pred %>%
  group_by(task) %>%
  summarise(
    AUC = tryCatch(
      {
        roc_obj <- pROC::roc(
          response = factor(is_patient, levels = c(0, 1)),
          predictor = p_hat_task,
          quiet = TRUE
        )
        as.numeric(pROC::auc(roc_obj))
      },
      error = function(e) NA_real_
    )
  )
print(auc_by_task)

# ---------------------------
# 6) (Opzionale) confronto "stacking" per task:
#    rifitta 3 modelli filtrando per task e confronta LOO/stacking
# ---------------------------
refit_by_task <- function(task_name) {
  dat_t <- dat_long %>% filter(task == task_name)
  brm(
    is_patient ~ 0 + type + z_value + (1 + z_value | feature) + (1 | user_id),
    data = dat_t,
    prior = c(
      prior(normal(0, 1), class = "b"),
      prior(student_t(3, 0, 2.5), class = "sd"),
      prior(lkj(2), class = "cor")
    ),
    cores = max(2, parallel::detectCores() - 1),
    chains = 4,
    iter = 4000,
    warmup = 1500,
    seed = 20251013,
    control = list(adapt_delta = 0.97, max_treedepth = 12)
  )
}

fit_wcst <- refit_by_task("wcst")
fit_prl <- refit_by_task("prl")
fit_ts <- refit_by_task("ts")

loo_wcst <- loo(fit_wcst)
loo_prl <- loo(fit_prl)
loo_ts <- loo(fit_ts)
print(loo_compare(loo_wcst, loo_prl, loo_ts))

# Pesi di stacking (quale compito contribuisce di più alla predizione)
stacking_weights(list(loo_wcst, loo_prl, loo_ts))


library(brms)

bf_form2 <- bf(
  is_patient ~
    0 +
      task +
      type +
      task:type +
      z_value:task +
      z_value:type +
      (1 + z_value | feature) + # parziale pooling per feature
      (1 | user_id), # intercetto soggetto (necessario, ma non dominante)
  family = bernoulli()
)

priors2 <- c(
  prior(normal(0, 1), class = "b"),
  prior(student_t(3, 0, 2.5), class = "sd"),
  prior(lkj(2), class = "cor")
)

fit2 <- brm(
  formula = bf_form2,
  data = dat_long, # lo stesso dat_long che hai già
  prior = priors2,
  chains = 4,
  iter = 4000,
  warmup = 1500,
  seed = 20251013,
  cores = max(2, parallel::detectCores() - 1),
  control = list(adapt_delta = 0.97, max_treedepth = 12),
  backend = "cmdstanr"
)

loo_fit2 <- loo(fit2)
print(loo_fit2)

# Contrasti fra task nei coefficienti su z_value (informatività differenziale per task):
hypothesis(fit2, "z_value:taskprl - z_value:taskwcst > 0")
hypothesis(fit2, "z_value:taskts  - z_value:taskwcst > 0")
hypothesis(fit2, "z_value:taskts  - z_value:taskprl  > 0")


#######
str(dat_long)

rosenberg <- quest |>
  dplyr::select(ros_tot, subj_code) |>
  dplyr::rename(
    user_id = subj_code
  )

dat_long2 <- left_join(dat_long, rosenberg, by = "user_id") %>%
  mutate(z_ros = scale(ros_tot))

bf_mod <- bf(
  is_patient ~
    0 +
      task +
      type +
      task:type +
      z_value * task +
      z_value:z_ros +
      z_value:task:z_ros +
      (1 + z_value | feature) +
      (1 | user_id),
  family = bernoulli()
)

fit_mod <- brm(
  bf_mod,
  data = dat_long2,
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(student_t(3, 0, 2.5), class = "sd"),
    prior(lkj(2), class = "cor")
  ),
  backend = "cmdstanr",
  algorithm = "meanfield"
  # chains=4, iter=4000, warmup=1500, seed=2,
  # control=list(adapt_delta=0.97, max_treedepth=12)
)

# Test: la moderazione del task PRL da SelfEsteem
hypothesis(fit_mod, "taskprl:z_value:z_ros> 0")
