# 08_sticky.R --------------------------------------------------------------
# Stima di hmm_sticky, confronto LOO con hmm_hier, PPC, affidabilita' e
# recovery dei tre parametri per soggetto. Tutto quello che serve per decidere
# se la perseverazione va nel modello finale.

suppressPackageStartupMessages({ library(cmdstanr); library(posterior); library(loo) })
set_cmdstan_path(file.path(Sys.getenv("CONDA_PREFIX"), "bin", "cmdstan"))
source("wcst_rl/funs_simulate.R")

sdh <- readRDS("wcst_rl/wcst_stan_data_hier.RDS")
N <- sdh$N; g <- sdh$grp
stan_input <- list(N = N, T = sdh$T, choice = sdh$choice, kc = sdh$k_color,
                   ks = sdh$k_shape, kn = sdh$k_number, rew = sdh$rew,
                   valid = sdh$valid, grp = as.numeric(g))
init3 <- function() list(mu = c(-1.5, 1.8, 0.5), bgrp = rep(0, 3),
                         sigma = rep(0.3, 3), z = matrix(0, 3, N),
                         Lcorr = diag(3), mu_lapse = -4.5, mu_eta = -2.8)

mod <- cmdstan_model("wcst_rl/stan/hmm_sticky.stan")
fit <- mod$sample(data = stan_input, chains = 4, parallel_chains = 4,
                  iter_warmup = 1000, iter_sampling = 1000, seed = 20260905,
                  adapt_delta = 0.95, max_treedepth = 12, refresh = 500, init = init3)
fit$save_object("wcst_rl/fits/hmm_sticky.RDS")
dg <- fit$diagnostic_summary()
message("divergenze = ", sum(dg$num_divergent), " | R-hat max = ",
        round(max(fit$summary()$rhat, na.rm = TRUE), 4))
print(as.data.frame(fit$summary(c("mu", "bgrp", "sigma", "mu_lapse", "mu_eta"))[
  , c("variable", "mean", "q5", "q95", "rhat", "ess_bulk")]), digits = 3)

# --- LOO contro hmm_hier ---------------------------------------------------
ll_new <- fit$draws("log_lik", format = "draws_matrix")
cid <- rep(1:4, each = nrow(ll_new) / 4)
loo_new <- loo(ll_new, r_eff = relative_eff(exp(ll_new), chain_id = cid))
old <- readRDS("wcst_rl/fits/hmm_hier.RDS")
ll_old <- old$draws("log_lik", format = "draws_matrix")
loo_old <- loo(ll_old, r_eff = relative_eff(exp(ll_old), chain_id = cid))
rm(old, ll_old); gc()
cat("\n--- LOO: hmm_sticky vs hmm_hier ---\n")
print(loo_compare(list(hmm_sticky = loo_new, hmm_hier = loo_old)))
cat("elpd per trial, hmm_sticky:",
    round(loo_new$estimates["elpd_loo", "Estimate"] / (N * sdh$T), 4),
    "| Pareto k > 0.7:", sum(loo_new$diagnostics$pareto_k > 0.7), "\n")

# --- PPC -------------------------------------------------------------------
obs <- wcst_signatures(sdh, sdh$choice, sdh$rew)
SIG <- c("acc", "prop_pers_err", "prop_non_pers_err", "prop_rep_dim",
         "win_stay", "lose_shift", "acc_pos1_3", "acc_pos8_10")
dh <- fit$draws("h", format = "draws_matrix"); dd <- fit$draws("d", format = "draws_matrix")
dk <- fit$draws("kappa", format = "draws_matrix")
dl <- as.vector(fit$draws("lapse", format = "draws_matrix"))
de <- as.vector(fit$draws("eta", format = "draws_matrix"))
NREP <- 200; set.seed(11); idx <- sample(nrow(dh), NREP)
rep_group <- matrix(NA_real_, NREP, length(SIG), dimnames = list(NULL, SIG))
rep_subj <- array(NA_real_, c(NREP, N, length(SIG)))
for (k in seq_along(idx)) {
  j <- idx[k]
  s <- sim_hmm_sticky(sdh, dh[j, ], dd[j, ], dk[j, ], rep(dl[j], N), rep(de[j], N),
                      seed = 2000 + k)
  sg <- wcst_signatures(sdh, s$choice, s$rew)
  rep_group[k, ] <- colMeans(sg[, SIG], na.rm = TRUE)
  rep_subj[k, , ] <- as.matrix(sg[, SIG])
  if (k %% 100 == 0) message("  replicazione ", k, "/", NREP)
}
ppc <- data.frame(indice = SIG,
                  osservato = sapply(SIG, function(v) mean(obs[[v]], na.rm = TRUE)),
                  predetto = colMeans(rep_group, na.rm = TRUE),
                  q05 = apply(rep_group, 2, quantile, .05, na.rm = TRUE),
                  q95 = apply(rep_group, 2, quantile, .95, na.rm = TRUE))
ppc$dentro_90 <- ppc$osservato >= ppc$q05 & ppc$osservato <= ppc$q95
ppc$r_individuale <- sapply(seq_along(SIG), function(k)
  suppressWarnings(cor(obs[[SIG[k]]], colMeans(rep_subj[, , k], na.rm = TRUE),
                       use = "complete.obs")))
cat("\n--- PPC hmm_sticky ---\n"); print(ppc, digits = 3, row.names = FALSE)

# --- affidabilita' modello-based e peso del gruppo -------------------------
PN <- c("logit_h" = 1, "log_d" = 2, "kappa" = 3)
sig <- fit$draws("sigma", format = "draws_matrix")
rel <- do.call(rbind, lapply(names(PN), function(v) {
  dr <- if (v == "logit_h") qlogis(dh) else if (v == "log_d") log(dd) else dk
  tau2 <- mean(sig[, PN[[v]]]^2); vw <- mean(apply(dr, 2, var))
  pm <- colMeans(dr)
  data.frame(par = v, sd_between = sqrt(tau2), sd_within = sqrt(vw),
             reliability = tau2 / (tau2 + vw),
             var_da_gruppo = summary(lm(pm ~ g))$r.squared)
}))
cat("\n--- affidabilita' e quota di varianza spiegata dal gruppo ---\n")
print(rel, digits = 3, row.names = FALSE)

# --- recovery dei tre parametri -------------------------------------------
tr_h <- colMeans(dh); tr_d <- colMeans(dd); tr_k <- colMeans(dk)
sim <- sim_hmm_sticky(sdh, tr_h, tr_d, tr_k, rep(mean(dl), N), rep(mean(de), N), seed = 777)
f2 <- mod$sample(data = c(stan_input[c("N", "T", "kc", "ks", "kn", "grp")],
                          list(choice = sim$choice, rew = sim$rew,
                               valid = matrix(1L, N, sdh$T))),
                 chains = 4, parallel_chains = 4, iter_warmup = 1000,
                 iter_sampling = 1000, seed = 20260905, adapt_delta = 0.95,
                 refresh = 0, init = init3)
message("recovery, R-hat max = ", round(max(f2$summary()$rhat, na.rm = TRUE), 3))
rec <- data.frame(
  par = c("h", "log d", "kappa"),
  r_grezza = c(cor(tr_h, colMeans(f2$draws("h", format = "draws_matrix"))),
               cor(log(tr_d), log(colMeans(f2$draws("d", format = "draws_matrix")))),
               cor(tr_k, colMeans(f2$draws("kappa", format = "draws_matrix")))))
rec$r_netto_gruppo <- c(
  cor(resid(lm(tr_h ~ g)), resid(lm(colMeans(f2$draws("h", format = "draws_matrix")) ~ g))),
  cor(resid(lm(log(tr_d) ~ g)), resid(lm(log(colMeans(f2$draws("d", format = "draws_matrix"))) ~ g))),
  cor(resid(lm(tr_k ~ g)), resid(lm(colMeans(f2$draws("kappa", format = "draws_matrix")) ~ g))))
cat("\n--- parameter recovery hmm_sticky ---\n"); print(rec, digits = 3, row.names = FALSE)

write.csv(ppc, "wcst_rl/ppc_sticky.csv", row.names = FALSE)
write.csv(rel, "wcst_rl/reliability_sticky.csv", row.names = FALSE)
write.csv(rec, "wcst_rl/recovery_sticky.csv", row.names = FALSE)
saveRDS(list(ppc = ppc, rel = rel, rec = rec, obs = obs, rep_group = rep_group,
             rep_subj = rep_subj, SIG = SIG, loo = loo_new),
        "wcst_rl/fits/sticky_summary.RDS")
message("\nsalvati ppc_sticky.csv, reliability_sticky.csv, recovery_sticky.csv")
