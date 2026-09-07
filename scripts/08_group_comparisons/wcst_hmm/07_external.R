# 07_external.R ------------------------------------------------------------
# Validita' esterna dei parametri e confronto tra gruppi.
# (1) effetto di gruppo stimato dentro il modello (b_h, b_d): nessun test su
#     stime puntuali, l'incertezza e' quella a posteriori;
# (2) validita' convergente: correlazione con gli indici comportamentali
#     classici, calcolata draw per draw;
# (3) correlazione con i parametri di PRL e task-switching (stime puntuali:
#     l'incertezza propagata e' solo quella dei parametri WCST).

suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })
source("wcst_rl/funs_simulate.R")

sdh <- readRDS("wcst_rl/wcst_stan_data_hier.RDS")
fit <- readRDS("wcst_rl/fits/hmm_hier.RDS")
N <- sdh$N
ids <- sdh$user_id

# --- (1) effetto di gruppo -------------------------------------------------
cat("--- effetto di gruppo (pazienti - controlli), scala non vincolata ---\n")
grp_tab <- do.call(rbind, lapply(c("b_h", "b_d"), function(v) {
  x <- as.vector(fit$draws(v, format = "draws_matrix"))
  data.frame(par = v, media = mean(x), q2.5 = quantile(x, .025),
             q97.5 = quantile(x, .975), p_dir = mean(x > 0))
}))
print(grp_tab, digits = 3, row.names = FALSE)

# dimensione dell'effetto standardizzata sulla SD tra soggetti
sig <- fit$draws("sigma", format = "draws_matrix")
cat("\neffetto standardizzato (b / sigma):\n")
for (k in 1:2) {
  v <- c("b_h", "b_d")[k]
  x <- as.vector(fit$draws(v, format = "draws_matrix")) / as.vector(sig[, k])
  cat(sprintf("  %s: %.2f [%.2f, %.2f]\n", v, mean(x), quantile(x, .025), quantile(x, .975)))
}

# --- (2) validita' convergente --------------------------------------------
obs <- wcst_signatures(sdh, sdh$choice, sdh$rew)
bw <- read.csv("wcst_behav_indices.csv")           # indici salvati dalla pipeline originale
bw <- bw[match(ids, bw$subj_name), ]

dh <- fit$draws("logit_h", format = "draws_matrix")
dd <- fit$draws("log_d",   format = "draws_matrix")
cor_draws <- function(dr, y) {
  ok <- !is.na(y)
  apply(dr[, ok, drop = FALSE], 1, function(x) cor(x, y[ok], method = "spearman"))
}
targets <- list(
  `accuratezza`            = obs$acc,
  `err. perseverativi`     = obs$prop_pers_err,
  `err. non persev.`       = obs$prop_non_pers_err,
  `ripetizione dimensione` = obs$prop_rep_dim,
  `lose-shift`             = obs$lose_shift,
  `recupero post-switch`   = obs$acc_pos1_3,
  `err. persev. (pipeline)` = bw$prop_pers_err)
conv <- do.call(rbind, lapply(names(targets), function(nm) {
  a <- cor_draws(dh, targets[[nm]]); b <- cor_draws(dd, targets[[nm]])
  data.frame(indice = nm,
             r_h = mean(a), h_q2.5 = quantile(a, .025), h_q97.5 = quantile(a, .975),
             r_d = mean(b), d_q2.5 = quantile(b, .025), d_q97.5 = quantile(b, .975))
}))
cat("\n--- validita' convergente (rho di Spearman, media a posteriori) ---\n")
print(conv, digits = 2, row.names = FALSE)

# --- (3) correlazioni tra compiti -----------------------------------------
prl <- read.csv("prl_params.csv"); ts <- read.csv("task_switching_params.csv")
prl <- prl[match(ids, prl$user_id), ]; ts <- ts[match(ids, ts$user_id), ]
cat("\nsoggetti WCST con dati PRL:", sum(!is.na(prl$user_id)),
    "| con task-switching:", sum(!is.na(ts$user_id)), "\n")
cross_targets <- c(setNames(as.list(prl[, c("alpha", "pos_alpha", "t", "v")]),
                            paste0("PRL: ", c("alpha", "pos_alpha", "t", "v"))),
                   setNames(as.list(ts[, c("a_1", "v_1", "t_1")]),
                            paste0("TS: ", c("a_1", "v_1", "t_1"))))
cross <- do.call(rbind, lapply(names(cross_targets), function(nm) {
  y <- as.numeric(cross_targets[[nm]])
  a <- cor_draws(dh, y); b <- cor_draws(dd, y)
  data.frame(parametro = nm, n = sum(!is.na(y)),
             r_h = mean(a), h_q2.5 = quantile(a, .025), h_q97.5 = quantile(a, .975),
             r_d = mean(b), d_q2.5 = quantile(b, .025), d_q97.5 = quantile(b, .975))
}))
cat("\n--- correlazioni tra compiti ---\n"); print(cross, digits = 2, row.names = FALSE)

# --- esportazione dei parametri per soggetto ------------------------------
out <- data.frame(
  user_id = ids, group = sdh$group_label,
  h_mean = colMeans(fit$draws("h", format = "draws_matrix")),
  h_sd   = apply(fit$draws("h", format = "draws_matrix"), 2, sd),
  d_mean = colMeans(fit$draws("d", format = "draws_matrix")),
  d_sd   = apply(fit$draws("d", format = "draws_matrix"), 2, sd),
  logit_h_mean = colMeans(dh), logit_h_sd = apply(dh, 2, sd),
  log_d_mean   = colMeans(dd), log_d_sd   = apply(dd, 2, sd))
out <- cbind(out, obs[, c("acc", "prop_pers_err", "prop_non_pers_err",
                          "prop_rep_dim", "win_stay", "lose_shift",
                          "acc_pos1_3", "acc_pos8_10")])
write.csv(out, "wcst_rl/wcst_hmm_params_hier.csv", row.names = FALSE)
write.csv(conv, "wcst_rl/convergent_validity.csv", row.names = FALSE)
write.csv(cross, "wcst_rl/cross_task.csv", row.names = FALSE)
write.csv(grp_tab, "wcst_rl/group_effects.csv", row.names = FALSE)
message("\nsalvati wcst_hmm_params_hier.csv, convergent_validity.csv, cross_task.csv, group_effects.csv")
