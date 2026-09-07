# ============================================================
# AUC con CV stratificata (ridge/enet) su EAT-26 (dieting, bulimia, oral_control)
# - Out-of-fold predictions
# - AUC + IC DeLong
# - CV ripetuta (robustezza) e grid su alpha
# - Salvataggio risultati e OOF su disco
# ============================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tibble, stringr, glmnet, pROC, purrr, readr, tidyr)

set.seed(20251013)


quest <- readRDS(
  here::here(
    "data", "quest", "quest_data_tot.rds"
  )
)


good_ids <- c(
  "cr_gi_1994_10_14_378_f", "gr_bo_1996_07_31_547_f", "an_de_1998_11_10_289_f",
  "al_zu_1997_04_02_880_f", "am_gu_1999_02_11_937_f", "ir_pi_2002_01_22_765_f",
  "el_ma_1986_06_14_839_f", "ir_bo_1981_03_29_325_f", "da_de_1998_08_15_141_m",
  "em_al_1989_07_27_200_f", "bi_an_2001_09_16_735_f", "em_gr_2002_08_25_628_f",
  "em_or_2003_01_02_101_f", "ch_ca_2000_09_26_406_f", "lu_mu_1997_03_18_059_f",
  "al_ro_1989_04_25_160_f", "an_am_1996_05_12_176_f", "he_ha_2006_04_21_874_f",
  "fe_ma_1998_06_29_257_f", "sa_ta_2003_11_14_150_f", "fe_sa_2002_05_09_08_f" ,
  "ch_pi_2004_02_25_126_f", "so_be_2008_12_15_399_f", "fr_la_2004_05_17_363_f",
  "ca_so_2001_01_09_118_f", "fr_bo_1993_09_09_170_f", "ch_ri_1993_05_05_564_f",
  "ch_pi_2001_10_08_418_f", "em_bi_2007_12_28_766_f", "fr_ro_1982_08_15_048_f",
  "ch_ma_2001_10_27_331_f", "gr_de_2002_09_21_426_f", "gi_va_1992_04_14_174_f",
  "fr_au_1987_12_16_221_f", "il_fu_2002_12_30_306_f", "ma_za_2002_02_28_051_f",
  "es_bo_2004_07_23_474_f", "ir_to_2007_08_01_838_f", "ar_ce_2005_04_20_937_f",
  "gi_ba_2008_01_31_376_f", "cl_pu_2007_05_24_423_f", "bi_di_2006_04_20_725_f",
  "gi_fi_1996_03_09_339_f", "ch_ca_1998_01_31_179_f", "mi_pr_1997_07_03_575_f",
  "fr_tr_1997_09_19_223_f", "el_gi_1994_06_25_443_f", "el_ch_1994_09_27_944_f",
  "fr_pl_2002_02_14_755_f", "si_sc_1992_12_23_119_f", "su_or_1993_01_13_765_f",
  "al_an_1996_06_03_205_f", "na_ge_1996_05_29_070_f", "ma_ce_2002_04_17_755_f",
  "gi_ma_1998_10_27_642_f", "ga_gi_1992_02_10_570_f", "gi_gi_1993_06_23_188_f",
  "el_ma_2001_07_17_978_f", "ca_sa_2001_06_01_608_f", "el_pa_2000_09_28_331_f",
  "el_ta_2000_01_02_570_f", "re_ve_2001_03_28_201_f", "an_ma_2001_06_27_832_f",
  "al_lo_2001_02_10_286_f", "el_nu_2001_02_09_373_f", "be_ba_1995_04_27_656_f",
  "al_me_2001_06_24_456_f", "la_sa_2001_11_15_307_f", "gi_to_1997_07_30_762_f",
  "fr_ba_1997_10_29_663_f", "so_go_1997_05_16_135_f", "al_la_2001_07_14_104_f",
  "ar_cr_1996_08_05_738_f", "ma_sp_2000_08_01_464_f", "el_la_1996_02_10_228_f",
  "sa_pa_2001_05_14_311_f", "gi_lu_2001_07_05_641_f", "ca_va_2001_08_28_797_f",
  "so_mo_2001_10_27_943_f", "ma_pe_2000_09_06_022_f", "ma_ca_1995_10_01_691_f",
  "gi_sa_2000_07_05_104_f", "an_am_1993_05_20_789_f", "la_or_2001_09_18_400_f",
  "fe_te_2000_05_17_086_f", "ar_lo_2001_11_02_365_f", "ma_pr_2000_09_18_430_f",
  "fe_ro_1999_10_25_558_f"
)

quest <- quest[quest$subj_code %in% good_ids, ]



# ----------------------------
# 0) Impostazioni I/O
# ----------------------------
out_dir <- "auc_results_eat26"
dir.create(out_dir, showWarnings = FALSE)

# ----------------------------
# 1) Selezione colonne da `quest`
# ----------------------------
stopifnot(exists("quest"))

needed <- c("dieting", "bulimia", "oral_control", "is_patient")
missing <- setdiff(needed, names(quest))
if (length(missing) > 0) stop("Mancano queste colonne in `quest`: ", paste(missing, collapse=", "))

raw <- quest %>%
  select(all_of(needed)) %>%
  # assicurati che i predittori siano numerici
  mutate(
    dieting = as.numeric(dieting),
    bulimia = as.numeric(bulimia),
    oral_control = as.numeric(oral_control)
  ) %>%
  drop_na(dieting, bulimia, oral_control, is_patient)

# ----------------------------
# 2) y come fattore (Controls, Patients) + controlli
# ----------------------------
y_raw <- raw$is_patient
if (is.numeric(y_raw)) {
  if (!all(y_raw %in% c(0,1))) stop("`is_patient` numerico ma non binario {0,1}.")
  y_fac <- factor(ifelse(y_raw == 1, "Patients", "Controls"),
                  levels = c("Controls","Patients"))
} else {
  y_chr <- tolower(as.character(y_raw))
  y_chr <- dplyr::case_when(
    y_chr %in% c("1","patient","patients","paziente","pazienti","an") ~ "Patients",
    y_chr %in% c("0","control","controls","controllo","controlli","hc") ~ "Controls",
    TRUE ~ NA_character_
  )
  if (anyNA(y_chr)) stop("Impossibile ricodificare `is_patient` in (Controls, Patients).")
  y_fac <- factor(y_chr, levels=c("Controls","Patients"))
}

tbl_y <- table(y_fac)
message("Distribuzione gruppi: ", paste(names(tbl_y), as.integer(tbl_y), collapse=" | "))
if (length(tbl_y) < 2) stop("C'è una sola classe in `is_patient` dopo pulizia: serve almeno una osservazione per classe.")

# ----------------------------
# 3) Matrice X (solo EAT-26)
# ----------------------------
pred_cols <- c("dieting","bulimia","oral_control")
X_all <- as.matrix(raw %>% select(all_of(pred_cols)))

# ----------------------------
# 4) Utility: fold stratificati + CV esterna
# ----------------------------
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
    X, y_fac,
    alpha = 0,
    K = 10,
    seed = 1L,
    inner_nfolds = 5,
    standardize = TRUE
) {
  stopifnot(is.factor(y_fac) && identical(levels(y_fac), c("Controls","Patients")))
  folds <- make_stratified_folds(y_fac, K = K, seed = seed)
  oof_pred <- rep(NA_real_, nrow(X))
  for (k in seq_len(K)) {
    tr <- which(folds != k)
    te <- which(folds == k)
    Xtr <- X[tr, , drop = FALSE]
    ytr <- (y_fac[tr] == "Patients")
    Xte <- X[te, , drop = FALSE]
    
    cvfit <- glmnet::cv.glmnet(
      x = Xtr, y = ytr,
      family = "binomial",
      alpha = alpha,
      type.measure = "auc",
      nfolds = inner_nfolds,
      standardize = standardize
    )
    oof_pred[te] <- as.numeric(predict(cvfit, newx = Xte, s = "lambda.min", type = "response"))
  }
  
  roc_obj <- pROC::roc(response = y_fac, predictor = oof_pred,
                       levels = c("Controls","Patients"), direction = "auto", quiet = TRUE)
  
  list(
    auc = as.numeric(pROC::auc(roc_obj)),
    ci_auc = as.numeric(pROC::ci.auc(roc_obj, method = "delong")),
    oof_pred = oof_pred,
    folds = folds,
    roc = roc_obj
  )
}

# ----------------------------
# 5) Un run base: Ridge (α=0) ed Elastic Net (α=0.5)
# ----------------------------
ridge_res <- cv_auc_glmnet(X_all, y_fac, alpha = 0,   K = 10, seed = 20251013)
enet_res  <- cv_auc_glmnet(X_all, y_fac, alpha = 0.5, K = 10, seed = 20251013)

cat(sprintf(
  "\nEAT-26 | Ridge  (10-fold OOF) AUC = %.3f  [%.3f, %.3f] (DeLong)\n",
  ridge_res$auc, ridge_res$ci_auc[1], ridge_res$ci_auc[3]
))
cat(sprintf(
  "EAT-26 | ENet α=0.5 (10-fold OOF) AUC = %.3f  [%.3f, %.3f] (DeLong)\n\n",
  enet_res$auc, enet_res$ci_auc[1], enet_res$ci_auc[3]
))

# Salva OOF (per confronti paired su stesso campione)
readr::write_csv(
  tibble(model = "eat26_ridge_alpha0", y = as.character(y_fac), oof_pred = ridge_res$oof_pred),
  file.path(out_dir, "eat26_oof_ridge.csv")
)
readr::write_csv(
  tibble(model = "eat26_enet_alpha0.5", y = as.character(y_fac), oof_pred = enet_res$oof_pred),
  file.path(out_dir, "eat26_oof_enet.csv")
)

# ----------------------------
# 6) CV ripetuta + grid α (robustezza)
# ----------------------------
alpha_grid <- c(0, 0.25, 0.5, 0.75, 1)
R <- 20  # aumenta se vuoi intervalli più stretti

set.seed(20251013)
res_rep <- purrr::map_dfr(alpha_grid, function(a) {
  aucs <- replicate(
    R,
    cv_auc_glmnet(X_all, y_fac, alpha = a, K = 10, seed = sample.int(1e7, 1))$auc
  )
  tibble(
    alpha = a,
    auc_mean = mean(aucs),
    auc_sd = sd(aucs),
    auc_q025 = quantile(aucs, .025),
    auc_q975 = quantile(aucs, .975),
    R = R
  )
}) %>% arrange(desc(auc_mean))

print(res_rep)
readr::write_csv(res_rep, file.path(out_dir, "eat26_auc_repeated_cv.csv"))

# ----------------------------
# 7) Confronto paired tra ENet e Ridge (bootstrap su OOF)
# ----------------------------
diff_auc_boot <- function(y_fac, pred1, pred2, B = 5000, seed = 1L) {
  set.seed(seed)
  n <- length(y_fac)
  stopifnot(length(pred1) == n, length(pred2) == n)
  diffs <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    r1 <- pROC::roc(y_fac[idx], pred1[idx], levels = c("Controls","Patients"), quiet = TRUE)
    r2 <- pROC::roc(y_fac[idx], pred2[idx], levels = c("Controls","Patients"), quiet = TRUE)
    diffs[b] <- as.numeric(pROC::auc(r2)) - as.numeric(pROC::auc(r1))
  }
  tibble(
    diff_mean = mean(diffs),
    diff_q025 = quantile(diffs, .025),
    diff_q975 = quantile(diffs, .975)
  )
}

cmp_tbl <- diff_auc_boot(
  y_fac,
  ridge_res$oof_pred,
  enet_res$oof_pred,
  B = 5000, seed = 123
)
print(cmp_tbl)
readr::write_csv(cmp_tbl, file.path(out_dir, "eat26_diff_auc_enet_vs_ridge.csv"))



rosenberg <- quest |> 
  dplyr::select(ros_tot, subj_code) |> 
  dplyr::rename(
    user_id = subj_code
  )








