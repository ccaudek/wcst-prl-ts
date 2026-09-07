# 09_external_sticky.R -----------------------------------------------------
# Validita' esterna, effetti di gruppo e correlazioni tra compiti sul modello
# finale (hmm_sticky). Esporta anche tutte le tabelle che servono alle figure.
#
# Ogni correlazione e' calcolata draw per draw: la stima per soggetto porta
# incertezza, e ignorarla gonfia la precisione apparente.

suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })
sdh <- readRDS("wcst_rl/wcst_stan_data_hier.RDS")
fit <- readRDS("wcst_rl/fits/hmm_sticky.RDS")
N <- sdh$N; g <- sdh$grp
source("wcst_rl/funs_simulate.R")

dh <- qlogis(fit$draws("h", format = "draws_matrix"))
dd <- log(fit$draws("d", format = "draws_matrix"))
dk <- fit$draws("kappa", format = "draws_matrix")
PAR <- list("logit h" = dh, "log d" = dd, "kappa" = dk)

# --- (1) effetti di gruppo -------------------------------------------------
bg <- fit$draws("bgrp", format = "draws_matrix")
sg <- fit$draws("sigma", format = "draws_matrix")
grp_tab <- do.call(rbind, lapply(1:3, function(k) {
  x <- as.vector(bg[, k]); s <- as.vector(x / sg[, k])
  data.frame(par = names(PAR)[k], b = mean(x), q05 = quantile(x, .05),
             q95 = quantile(x, .95), p_direzione = max(mean(x > 0), mean(x < 0)),
             b_su_sigma = median(s))
}))
cat("--- effetto di gruppo (pazienti - controlli) ---\n")
print(grp_tab, digits = 3, row.names = FALSE)

# --- (2) validita' convergente contro gli indici comportamentali ----------
obs <- wcst_signatures(sdh, sdh$choice, sdh$rew)
bw <- read.csv("wcst_behav_indices.csv")   # la chiave qui e' subj_name, non user_id
bw <- bw[match(sdh$user_id, bw$subj_name), ]
targets <- list(
  "accuratezza" = obs$acc,
  "errori perseverativi" = obs$prop_pers_err,
  "errori non perseverativi" = obs$prop_non_pers_err,
  "ripetizione dimensione" = obs$prop_rep_dim,
  "lose-shift" = obs$lose_shift,
  "recupero post-switch" = obs$acc_pos1_3,
  "err. persev. (pipeline)" = bw$prop_pers_err,
  "err. non persev. (pipeline)" = bw$prop_non_pers_err,
  "risp. persev. (pipeline)" = bw$prop_pers_resp)
cor_draws <- function(D, y) {
  keep <- !is.na(y)
  apply(D[, keep, drop = FALSE], 1, function(x) cor(x, y[keep], method = "spearman"))
}
conv <- do.call(rbind, lapply(names(targets), function(nm) {
  r <- lapply(PAR, cor_draws, y = targets[[nm]])
  data.frame(indice = nm,
             rho_h = mean(r[["logit h"]]), h_q05 = quantile(r[["logit h"]], .05),
             h_q95 = quantile(r[["logit h"]], .95),
             rho_d = mean(r[["log d"]]), rho_kappa = mean(r[["kappa"]]))
}))
cat("\n--- validita' convergente (Spearman, media a posteriori) ---\n")
print(conv, digits = 2, row.names = FALSE)

# --- (3) correlazioni tra compiti -----------------------------------------
prl <- read.csv("prl_params.csv"); ts <- read.csv("task_switching_params.csv")
prl <- prl[match(sdh$user_id, prl$user_id), ]; ts <- ts[match(sdh$user_id, ts$user_id), ]
cat("\nsoggetti WCST con dati PRL:", sum(!is.na(prl$user_id)),
    "| con dati task-switching:", sum(!is.na(ts$user_id)), "\n")
cross_src <- list()
for (v in setdiff(names(prl), c("user_id", "group", "group_label")))
  if (is.numeric(prl[[v]])) cross_src[[paste0("PRL: ", v)]] <- prl[[v]]
for (v in setdiff(names(ts), c("user_id", "group", "group_label")))
  if (is.numeric(ts[[v]])) cross_src[[paste0("TS: ", v)]] <- ts[[v]]
cross <- do.call(rbind, lapply(names(cross_src), function(nm) {
  y <- cross_src[[nm]]
  if (sum(!is.na(y)) < 30) return(NULL)
  r <- cor_draws(dh, y)
  data.frame(variabile = nm, n = sum(!is.na(y)), rho_h = mean(r),
             q05 = quantile(r, .05), q95 = quantile(r, .95),
             rho_d = mean(cor_draws(dd, y)))
}))
cross <- cross[order(-abs(cross$rho_h)), ]
cat("\n--- correlazioni tra compiti (ordinate per |rho| con h) ---\n")
print(head(cross, 10), digits = 2, row.names = FALSE)

# --- esportazioni ----------------------------------------------------------
pars <- data.frame(
  user_id = sdh$user_id, group = sdh$group_label,
  h = colMeans(fit$draws("h", format = "draws_matrix")),
  h_sd = apply(fit$draws("h", format = "draws_matrix"), 2, sd),
  d = colMeans(fit$draws("d", format = "draws_matrix")),
  d_sd = apply(fit$draws("d", format = "draws_matrix"), 2, sd),
  kappa = colMeans(dk), kappa_sd = apply(dk, 2, sd))
pars <- cbind(pars, obs[, c("acc", "prop_pers_err", "prop_non_pers_err",
                            "prop_rep_dim", "win_stay", "lose_shift",
                            "acc_pos1_3", "acc_pos8_10")])
write.csv(pars, "wcst_rl/wcst_params_final.csv", row.names = FALSE)
write.csv(grp_tab, "wcst_rl/group_effects_final.csv", row.names = FALSE)
write.csv(conv, "wcst_rl/convergent_validity_final.csv", row.names = FALSE)
write.csv(cross, "wcst_rl/cross_task_final.csv", row.names = FALSE)

# dati per le figure
write.csv(as.data.frame(bg), "wcst_rl/fig_bgrp_draws.csv", row.names = FALSE)
ss <- readRDS("wcst_rl/fits/sticky_summary.RDS")
write.csv(ss$ppc, "wcst_rl/fig_ppc_sticky.csv", row.names = FALSE)
ppc_ind <- data.frame(user_id = sdh$user_id, group = sdh$group_label)
for (k in seq_along(ss$SIG)) {
  ppc_ind[[paste0("obs_", ss$SIG[k])]] <- ss$obs[[ss$SIG[k]]]
  ppc_ind[[paste0("pred_", ss$SIG[k])]] <- colMeans(ss$rep_subj[, , k], na.rm = TRUE)
}
write.csv(ppc_ind, "wcst_rl/fig_ppc_individual.csv", row.names = FALSE)
rc <- read.csv("wcst_rl/recovery_sticky.csv")
write.csv(rc, "wcst_rl/fig_recovery.csv", row.names = FALSE)
message("\nesportati wcst_params_final.csv e le tabelle per le figure")
