# 04_reliability.R ---------------------------------------------------------
# Quanta della differenziazione tra soggetti e' reale?
# (1) affidabilita' modello-based: tau^2 / (tau^2 + varianza a posteriori media
#     entro soggetto), calcolata draw per draw cosi' che l'incertezza si
#     propaghi;
# (2) split-half: rifit sui blocchi 1-3-5 e 2-4-6, correlazione tra le stime
#     per soggetto, con correzione di Spearman-Brown.

suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })
set_cmdstan_path(file.path(Sys.getenv("CONDA_PREFIX"), "bin", "cmdstan"))

sdh <- readRDS("wcst_rl/wcst_stan_data_hier.RDS")
fit <- readRDS("wcst_rl/fits/hmm_hier.RDS")
N <- sdh$N

# --- (1) affidabilita' modello-based --------------------------------------
# parametri per soggetto sulla scala non vincolata, dove tau e' definito
draws_lh <- fit$draws("logit_h", format = "draws_matrix")
draws_ld <- fit$draws("log_d",   format = "draws_matrix")
sig <- fit$draws("sigma", format = "draws_matrix")

rel_table <- function(dr, tau) {
  var_within <- mean(apply(dr, 2, var))          # varianza a posteriori media entro soggetto
  tau2 <- mean(tau^2)
  c(sd_between = sqrt(tau2),
    sd_within  = sqrt(var_within),
    reliability = tau2 / (tau2 + var_within),
    sd_post_means = sd(colMeans(dr)))
}
rel <- rbind(`logit h` = rel_table(draws_lh, sig[, 1]),
             `log d`   = rel_table(draws_ld, sig[, 2]))
cat("--- affidabilita' modello-based ---\n"); print(round(rel, 3))

# --- (2) split-half --------------------------------------------------------
halves <- list(A = c(1, 3, 5), B = c(2, 4, 6))
mod <- cmdstan_model("wcst_rl/stan/hmm_hier.stan")
fits_h <- list()
for (nm in names(halves)) {
  keep <- which(sdh$block[1, ] %in% halves[[nm]])
  si <- list(N = N, T = length(keep),
             choice = sdh$choice[, keep], kc = sdh$k_color[, keep],
             ks = sdh$k_shape[, keep], kn = sdh$k_number[, keep],
             rew = sdh$rew[, keep], valid = sdh$valid[, keep],
             grp = as.numeric(sdh$grp))
  message("--- meta' ", nm, " (", length(keep), " trial) ---")
  f <- mod$sample(data = si, chains = 4, parallel_chains = 4,
                  iter_warmup = 1000, iter_sampling = 1000, seed = 20260905,
                  adapt_delta = 0.95, refresh = 0,
                  init = function() list(mu_h = -1.5, mu_d = 2.0, b_h = 0, b_d = 0,
                                         mu_lapse = -4.5, mu_eta = -2.8,
                                         sigma = c(.3, .3), z = matrix(0, 2, N),
                                         Lcorr = diag(2)))
  message("  R-hat max = ", round(max(f$summary()$rhat, na.rm = TRUE), 3))
  fits_h[[nm]] <- f
}
sh <- data.frame(
  par = c("logit h", "log d"),
  r_pearson = NA_real_, r_spearman = NA_real_, sb_corrected = NA_real_)
for (k in 1:2) {
  vn <- c("logit_h", "log_d")[k]
  a <- colMeans(fits_h$A$draws(vn, format = "draws_matrix"))
  b <- colMeans(fits_h$B$draws(vn, format = "draws_matrix"))
  r <- cor(a, b); rs <- cor(a, b, method = "spearman")
  sh[k, 2:4] <- c(r, rs, 2 * r / (1 + r))
}
cat("\n--- split-half (blocchi 1-3-5 vs 2-4-6, 30 trial per meta') ---\n")
print(sh, digits = 3, row.names = FALSE)

saveRDS(list(rel = rel, split_half = sh, fits_h = lapply(fits_h, function(f)
  list(logit_h = colMeans(f$draws("logit_h", format = "draws_matrix")),
       log_d   = colMeans(f$draws("log_d",   format = "draws_matrix"))))),
  "wcst_rl/fits/reliability.RDS")
write.csv(cbind(as.data.frame(rel), par = rownames(rel)), "wcst_rl/reliability_table.csv", row.names = FALSE)
message("\nsalvato wcst_rl/reliability_table.csv")
