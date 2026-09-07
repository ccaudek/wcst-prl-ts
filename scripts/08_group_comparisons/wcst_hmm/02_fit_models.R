# 02_fit_models.R ----------------------------------------------------------
# Compila e stima i tre modelli gerarchici; salva fit e log_lik per il LOO.

suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })
set_cmdstan_path(file.path(Sys.getenv("CONDA_PREFIX"), "bin", "cmdstan"))

sdh <- readRDS("wcst_rl/wcst_stan_data_hier.RDS")
stan_input <- list(
  N = sdh$N, T = sdh$T, choice = sdh$choice,
  kc = sdh$k_color, ks = sdh$k_shape, kn = sdh$k_number,
  rew = sdh$rew, valid = sdh$valid, grp = as.numeric(sdh$grp))

MODELS <- c("hmm_hier", "rw_dim_hier", "hmm_hier4")
dir.create("wcst_rl/fits", showWarnings = FALSE)

for (mm in MODELS) {
  message("\n=== ", mm, " ===")
  mod <- cmdstan_model(file.path("wcst_rl/stan", paste0(mm, ".stan")))
  # inizializzazioni ragionevoli: accorciano l'adattamento e tengono le catene
  # lontane dai bordi
  npar <- if (mm == "hmm_hier") 2 else 4
  init_f <- function() {
    base <- list(sigma = rep(0.3, npar), z = matrix(0, npar, sdh$N),
                 Lcorr = diag(npar))
    if (mm == "hmm_hier") {
      c(base, list(mu_h = -1.5, mu_d = 2.0, b_h = 0, b_d = 0,
                   mu_lapse = -4.5, mu_eta = -2.8))
    } else if (mm == "hmm_hier4") {
      c(base, list(mu = c(-1.5, 2.0, -4.5, -2.8), bgrp = rep(0, 4)))
    } else {
      c(base, list(mu = c(0, 0, 1.0, 0), bgrp = rep(0, 4)))
    }
  }
  fit <- mod$sample(data = stan_input, chains = 4, parallel_chains = 4,
                    iter_warmup = 1000, iter_sampling = 1000, seed = 20260905,
                    adapt_delta = 0.95, max_treedepth = 12, refresh = 500,
                    init = init_f)
  fit$save_object(file.path("wcst_rl/fits", paste0(mm, ".RDS")))
  dg <- fit$diagnostic_summary()
  message(mm, ": divergenze = ", sum(dg$num_divergent),
          " | treedepth = ", sum(dg$max_treedepth),
          " | E-BFMI min = ", round(min(dg$ebfmi), 3))
  pars <- if (mm == "rw_dim_hier") c("mu", "bgrp", "sigma") else
          if (mm == "hmm_hier")   c("mu_h","mu_d","b_h","b_d","sigma","mu_lapse","mu_eta") else
                                   c("mu","bgrp","sigma")
  s <- fit$summary(pars)
  print(as.data.frame(s[, c("variable","mean","q5","q95","rhat","ess_bulk")]), digits = 3)
  message("R-hat max (tutti i parametri): ",
          round(max(fit$summary()$rhat, na.rm = TRUE), 4))
}
message("\nfit salvati in wcst_rl/fits/")
