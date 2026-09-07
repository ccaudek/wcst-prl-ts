# ============================================================
# AUC con CV stratificata (ridge/enet) su wcst_hmm_params.csv
# - Out-of-fold predictions
# - AUC + IC DeLong
# - CV ripetuta (robustezza) e grid su alpha
# - Salvataggio risultati e OOF su disco
# ============================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readr, dplyr, tibble, stringr, glmnet, pROC, purrr)

# ----------------------------
# 0) Impostazioni e path
# ----------------------------
set.seed(20251003)
# file_in <- here::here(
#   "data",
#   "processed",
#   "wcst_hmm_params.csv"
# )
# file_in <- here::here(
#   "data",
#   "processed",
#   "wcst_behav_indices.csv"
# )

# file_in <- here::here(
#   "data",
#   "processed",
#   "prl_params.csv"
# )
# file_in <- here::here(
#   "data",
#   "processed",
#   "prl_behav_indices.csv"
# )

# file_in <- here::here(
#   "data",
#   "processed",
#   "task_switching_params.csv"
# )
file_in <- here::here(
  "data",
  "processed",
  "task_switching_behav_indices.csv"
)

out_dir <- "auc_results"
dir.create(out_dir, showWarnings = FALSE)

# ------------------------------------------------------------
# 1) Lettura dati e preparazione X,y
#    Assunzioni:
#    - esiste una colonna con le etichette di gruppo:
#         "group" o "group_label" (Controls/Patients oppure 0/1)
#    - le colonne predittori sono numeriche (p.es. logit_h, log_d, logit_lapse, logit_eta)
# ------------------------------------------------------------
raw <- readr::read_csv(file_in, show_col_types = FALSE)

# Trova colonna gruppo
group_col <- dplyr::coalesce(
  if ("group" %in% names(raw)) "group" else NA_character_,
  if ("group_label" %in% names(raw)) "group_label" else NA_character_,
  if ("is_patient" %in% names(raw)) "is_patient" else NA_character_
)
if (is.na(group_col))
  stop(
    "Non ho trovato una colonna 'group', 'group_label' o 'is_patient' nel CSV."
  )

# Rende la colonna y come fattore con livelli (Controls, Patients)
y_raw <- raw[[group_col]]

# Ricodifica robusta: gestisce 0/1, patients/controls, etc.
if (is.numeric(y_raw)) {
  # aspettati 0/1
  stopifnot(all(y_raw %in% c(0, 1)))
  y_fac <- factor(
    ifelse(y_raw == 1, "Patients", "Controls"),
    levels = c("Controls", "Patients")
  )
} else {
  y_chr <- as.character(y_raw)
  # normalizza testo
  y_chr <- stringr::str_to_lower(y_chr)
  y_chr <- dplyr::case_when(
    y_chr %in% c("patient", "patients", "paziente", "pazienti", "1", "an") ~
      "Patients",
    y_chr %in% c("control", "controls", "controllo", "controlli", "0", "hc") ~
      "Controls",
    TRUE ~ y_raw
  )
  y_fac <- factor(y_chr, levels = c("Controls", "Patients"))
  # se ancora non combacia, prova a invertirle per contare i casi
  if (any(is.na(y_fac)))
    stop(
      "Impossibile ricodificare la colonna di gruppo in (Controls, Patients)."
    )
}

# Seleziona predittori numerici, escludendo id e gruppo
drop_cols <- c(
  group_col,
  "id",
  "ID",
  "is_patient",
  "subject",
  "subject_id",
  "user_id"
)
num_cols <- names(raw)[sapply(raw, is.numeric)]
pred_cols <- setdiff(num_cols, intersect(num_cols, drop_cols))

if (length(pred_cols) == 0)
  stop("Nessun predittore numerico trovato (esclusi id/gruppo).")
message("Userò questi predittori: ", paste(pred_cols, collapse = ", "))

dat <- raw %>%
  dplyr::select(all_of(pred_cols)) %>%
  mutate(.group = y_fac) %>%
  # rimuovi righe con NA
  tidyr::drop_na()

X_all <- as.matrix(dat %>% dplyr::select(all_of(pred_cols)))
y_fac <- dat$.group

# Controllo bilanciamento
tbl_y <- table(y_fac)
message(
  "Distribuzione gruppi: ",
  paste(names(tbl_y), as.integer(tbl_y), collapse = " | ")
)

# ------------------------------------------------------------
# 2) Utilità: fold stratificati e CV
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

cv_auc_glmnet <- function(
  X,
  y_fac,
  alpha = 0,
  K = 10,
  seed = 1L,
  inner_nfolds = 5,
  standardize = TRUE
) {
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
      nfolds = inner_nfolds,
      standardize = standardize
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

# ------------------------------------------------------------
# 3) Un run base: Ridge (alpha=0) ed Elastic Net (alpha=0.5)
# ------------------------------------------------------------
ridge_res <- cv_auc_glmnet(X_all, y_fac, alpha = 0, K = 10, seed = 20251003)
enet_res <- cv_auc_glmnet(X_all, y_fac, alpha = 0.5, K = 10, seed = 20251003)

cat(sprintf(
  "\nRidge  (10-fold OOF) AUC = %.3f  [%.3f, %.3f] (DeLong)\n",
  ridge_res$auc,
  ridge_res$ci_auc[1],
  ridge_res$ci_auc[3]
))
cat(sprintf(
  "ENet α=0.5 (10-fold OOF) AUC = %.3f  [%.3f, %.3f] (DeLong)\n\n",
  enet_res$auc,
  enet_res$ci_auc[1],
  enet_res$ci_auc[3]
))

# Salva OOF predizioni (utile per confronti diretti tra modelli sullo stesso sample)
readr::write_csv(
  tibble(
    model = "ridge_alpha0",
    y = as.character(y_fac),
    oof_pred = ridge_res$oof_pred
  ),
  # file.path(out_dir, "wcst_params_oof_ridge.csv")
  # file.path(out_dir, "wcst_behav_idices_oof_ridge.csv")
  # file.path(out_dir, "prl_params_oof_ridge.csv")
  file.path(out_dir, "prl_behav_indices_oof_ridge.csv")
)
readr::write_csv(
  tibble(
    model = "enet_alpha0.5",
    y = as.character(y_fac),
    oof_pred = enet_res$oof_pred
  ),
  # file.path(out_dir, "wcst_params_oof_enet.csv")
  # file.path(out_dir, "wcst_behav_idices_oof_enet.csv")
  # file.path(out_dir, "prl_params_oof_enet.csv")
  file.path(out_dir, "prl_behav_indices_oof_enet.csv")
)

# ------------------------------------------------------------
# 4) CV ripetuta + grid alpha (robustezza)
# ------------------------------------------------------------
alpha_grid <- c(0, 0.25, 0.5, 0.75, 1)
R <- 10 # aumenta se vuoi maggiore precisione

set.seed(20251003)
res_rep <- purrr::map_dfr(alpha_grid, function(a) {
  aucs <- replicate(
    R,
    cv_auc_glmnet(
      X_all,
      y_fac,
      alpha = a,
      K = 10,
      seed = sample.int(1e7, 1)
    )$auc
  )
  tibble(
    alpha = a,
    auc_mean = mean(aucs),
    auc_sd = sd(aucs),
    auc_q025 = quantile(aucs, .025),
    auc_q975 = quantile(aucs, .975),
    R = R
  )
}) %>%
  arrange(desc(auc_mean))

print(res_rep)

# readr::write_csv(res_rep, file.path(out_dir, "wcst_params_auc_repeated_cv.csv"))
# readr::write_csv(
#   res_rep,
#   file.path(out_dir, "wcst_behav_indices_auc_repeated_cv.csv")
# )
# readr::write_csv(res_rep, file.path(out_dir, "prl_params_auc_repeated_cv.csv"))
readr::write_csv(
  res_rep,
  file.path(out_dir, "prl_behav_indices_auc_repeated_cv.csv")
)


# ------------------------------------------------------------
# 5) (Opzionale) Confronto tra due modelli sullo stesso campione:
#    differenza AUC via bootstrap su predizioni OOF accoppiate
# ------------------------------------------------------------
diff_auc_boot <- function(y_fac, pred1, pred2, B = 5000, seed = 1L) {
  # bootstrap sui casi (paired), stima distribuzione AUC2 - AUC1
  set.seed(seed)
  n <- length(y_fac)
  stopifnot(length(pred1) == n, length(pred2) == n)
  # AUC osservate
  roc1 <- pROC::roc(
    y_fac,
    pred1,
    levels = c("Controls", "Patients"),
    quiet = TRUE
  )
  roc2 <- pROC::roc(
    y_fac,
    pred2,
    levels = c("Controls", "Patients"),
    quiet = TRUE
  )
  auc1 <- as.numeric(pROC::auc(roc1))
  auc2 <- as.numeric(pROC::auc(roc2))
  # bootstrap
  diffs <- numeric(B)
  for (b in 1:B) {
    idx <- sample.int(n, n, replace = TRUE)
    rb1 <- pROC::roc(
      y_fac[idx],
      pred1[idx],
      levels = c("Controls", "Patients"),
      quiet = TRUE
    )
    rb2 <- pROC::roc(
      y_fac[idx],
      pred2[idx],
      levels = c("Controls", "Patients"),
      quiet = TRUE
    )
    diffs[b] <- as.numeric(pROC::auc(rb2)) - as.numeric(pROC::auc(rb1))
  }
  tibble(
    auc1 = auc1,
    auc2 = auc2,
    diff_mean = mean(diffs),
    diff_q025 = quantile(diffs, .025),
    diff_q975 = quantile(diffs, .975)
  )
}

# Esempio: ENet (alpha=.5) vs Ridge
cmp_tbl <- diff_auc_boot(
  y_fac,
  ridge_res$oof_pred,
  enet_res$oof_pred,
  B = 5000,
  seed = 123
)
print(cmp_tbl)
readr::write_csv(
  cmp_tbl,
  file.path(out_dir, "wcst_params_diff_auc_enet_vs_ridge.csv")
)

# ------------------------------------------------------------
# 6) Consigli per replicare su altri dataset:
#    - Per "PRL params": punta file_in al CSV corrispondente; pred_cols auto-detect come sopra.
#    - Per "indici comportamentali": stesso codice (basta che le colonne siano numeriche).
#    - Mantieni nomenclatura file in output in modo coerente (es. prl_params_*.csv, wcst_indices_*.csv)
# ------------------------------------------------------------
