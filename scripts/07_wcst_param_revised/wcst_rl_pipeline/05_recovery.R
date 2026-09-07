# 05_recovery.R ------------------------------------------------------------
# Parameter recovery: si simulano 88 soggetti con i parametri stimati sui dati
# reali (quindi con la stessa dispersione individuale osservata), si rifitta il
# modello sui dati simulati e si correlano valori veri e recuperati.
# E' il test decisivo: se il recupero e' scarso, le differenze individuali
# stimate non sono interpretabili anche se il modello vince il LOO.

suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })
set_cmdstan_path(file.path(Sys.getenv("CONDA_PREFIX"), "bin", "cmdstan"))
source("wcst_rl/funs_simulate.R")

sdh <- readRDS("wcst_rl/wcst_stan_data_hier.RDS")
fit <- readRDS("wcst_rl/fits/hmm_hier.RDS")
N <- sdh$N

true_h <- colMeans(fit$draws("h", format = "draws_matrix"))
true_d <- colMeans(fit$draws("d", format = "draws_matrix"))
lapse  <- rep(mean(fit$draws("lapse", format = "draws_matrix")), N)
eta    <- rep(mean(fit$draws("eta",   format = "draws_matrix")), N)

sim <- sim_hmm(sdh, true_h, true_d, lapse, eta, seed = 4242)
message("accuratezza simulata: ", round(mean(sim$rew), 3),
        " | osservata: ", round(mean(sdh$rew), 3))

si <- list(N = N, T = sdh$T, choice = sim$choice, kc = sdh$k_color,
           ks = sdh$k_shape, kn = sdh$k_number, rew = sim$rew,
           valid = matrix(1L, N, sdh$T), grp = as.numeric(sdh$grp))
mod <- cmdstan_model("wcst_rl/stan/hmm_hier.stan")
f2 <- mod$sample(data = si, chains = 4, parallel_chains = 4,
                 iter_warmup = 1000, iter_sampling = 1000, seed = 20260905,
                 adapt_delta = 0.95, refresh = 500,
                 init = function() list(mu_h = -1.5, mu_d = 2.0, b_h = 0, b_d = 0,
                                        mu_lapse = -4.5, mu_eta = -2.8,
                                        sigma = c(.3, .3), z = matrix(0, 2, N),
                                        Lcorr = diag(2)))
message("R-hat max = ", round(max(f2$summary()$rhat, na.rm = TRUE), 3))

rec_h <- colMeans(f2$draws("h", format = "draws_matrix"))
rec_d <- colMeans(f2$draws("d", format = "draws_matrix"))
rec <- data.frame(
  par = c("h", "d"),
  r_pearson  = c(cor(true_h, rec_h), cor(log(true_d), log(rec_d))),
  r_spearman = c(cor(true_h, rec_h, method = "spearman"),
                 cor(true_d, rec_d, method = "spearman")))
cat("--- parameter recovery ---\n"); print(rec, digits = 3, row.names = FALSE)
# recupero anche dell'effetto di gruppo (che nei dati simulati e' quello stimato)
cat("\neffetto di gruppo, vero vs recuperato:\n")
for (v in c("b_h", "b_d"))
  cat(sprintf("  %s: vero %.3f | recuperato %.3f [%.3f, %.3f]\n", v,
      mean(fit$draws(v, format = "draws_matrix")),
      mean(f2$draws(v, format = "draws_matrix")),
      quantile(f2$draws(v, format = "draws_matrix"), .05),
      quantile(f2$draws(v, format = "draws_matrix"), .95)))

saveRDS(list(true_h = true_h, true_d = true_d, rec_h = rec_h, rec_d = rec_d, tab = rec),
        "wcst_rl/fits/recovery.RDS")
write.csv(rec, "wcst_rl/recovery_table.csv", row.names = FALSE)
message("\nsalvato wcst_rl/recovery_table.csv")
