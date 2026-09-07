# 10_nogroup.R -------------------------------------------------------------
# Le stime per soggetto del modello gerarchico sono restringite verso la media
# del PROPRIO gruppo: usarle per classificare i gruppi e' circolare. Qui lo
# stesso modello viene stimato con una sola popolazione (grp = 0 per tutti),
# quindi le stime per soggetto non contengono l'etichetta. Sono queste che
# vanno usate per tutte le analisi a livello individuale; la stima dell'effetto
# di gruppo resta quella del modello con il predittore.
suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })
set_cmdstan_path(file.path(Sys.getenv("CONDA_PREFIX"), "bin", "cmdstan"))

sdh <- readRDS("wcst_rl/wcst_stan_data_hier.RDS")
N <- sdh$N
stan_input <- list(N = N, T = sdh$T, choice = sdh$choice, kc = sdh$k_color,
                   ks = sdh$k_shape, kn = sdh$k_number, rew = sdh$rew,
                   valid = sdh$valid, grp = rep(0, N))
init3 <- function() list(mu = c(-1.5, 1.8, 0.5), bgrp = rep(0, 3),
                         sigma = rep(0.3, 3), z = matrix(0, 3, N),
                         Lcorr = diag(3), mu_lapse = -4.5, mu_eta = -2.8)

mod <- cmdstan_model("wcst_rl/stan/hmm_sticky.stan")
fit <- mod$sample(data = stan_input, chains = 4, parallel_chains = 4,
                  iter_warmup = 1000, iter_sampling = 1000, seed = 20260905,
                  adapt_delta = 0.95, max_treedepth = 12, refresh = 0, init = init3)
fit$save_object("wcst_rl/fits/hmm_sticky_nogroup.RDS")
dg <- fit$diagnostic_summary()
ss <- fit$summary(c("mu", "sigma", "mu_lapse", "mu_eta"))
cat("divergenze:", sum(dg$num_divergent), "| R-hat max:", round(max(ss$rhat, na.rm = TRUE), 4),
    "| ESS min:", round(min(ss$ess_bulk, na.rm = TRUE)), "\n")

dh <- fit$draws("h", format = "draws_matrix")
dd <- fit$draws("d", format = "draws_matrix")
dk <- fit$draws("kappa", format = "draws_matrix")
lg <- function(x) log(x / (1 - x))
out <- data.frame(
  user_id = sdh$user_id,
  wcst_logit_h    = colMeans(lg(dh)), wcst_logit_h_sd = apply(lg(dh), 2, sd),
  wcst_log_d      = colMeans(log(dd)), wcst_log_d_sd  = apply(log(dd), 2, sd),
  wcst_kappa      = colMeans(dk),      wcst_kappa_sd  = apply(dk, 2, sd))
write.csv(out, "wcst_rl/wcst_params_nogroup.csv", row.names = FALSE)

# quanta varianza tra soggetti resta spiegata dal gruppo con stime non informate
g <- as.numeric(sdh$grp)
r2 <- sapply(c("wcst_logit_h", "wcst_log_d", "wcst_kappa"),
             function(v) summary(lm(out[[v]] ~ g))$r.squared)
cat("\nR2 del gruppo sulle stime senza predittore:\n"); print(round(r2, 3))
cat("\ncorrelazione con le stime del modello con gruppo:\n")
og <- read.csv("wcst_rl/wcst_params_final.csv")
og <- og[match(out$user_id, og$user_id), ]
print(round(sapply(c("wcst_logit_h", "wcst_log_d", "wcst_kappa"),
                   function(v) cor(out[[v]], og[[v]])), 3))
cat("\nsd a posteriori media (senza gruppo):\n")
print(round(colMeans(out[, c("wcst_logit_h_sd", "wcst_log_d_sd", "wcst_kappa_sd")]), 3))
message("salvato wcst_params_nogroup.csv")
