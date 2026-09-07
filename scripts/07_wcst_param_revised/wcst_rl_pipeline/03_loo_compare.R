# 03_loo_compare.R ---------------------------------------------------------
# Confronto tra modelli su base PSIS-LOO, a livello di trial e di soggetto.
#
# Nota: il LOO per trial in un modello sequenziale e' approssimato (i trial di
# uno stesso soggetto non sono scambiabili). Lo riportiamo perche' e' lo
# standard in letteratura, ma la conclusione va appoggiata anche al LOO per
# soggetto, che rispetta la struttura di scambiabilita' del modello gerarchico.

suppressPackageStartupMessages({ library(cmdstanr); library(loo); library(posterior) })

sdh <- readRDS("wcst_rl/wcst_stan_data_hier.RDS")
MODELS <- c("hmm_hier", "hmm_hier4", "rw_dim_hier")
NT <- sdh$N * sdh$T

loos_trial <- list(); loos_subj <- list(); tab <- list()
for (mm in MODELS) {
  fit <- readRDS(file.path("wcst_rl/fits", paste0(mm, ".RDS")))
  ll_all <- fit$draws("log_lik", format = "draws_matrix")  # draws x (N*T)
  rm(fit); gc()
  keep <- apply(ll_all, 2, function(x) any(x != 0))        # scarta i trial non validi
  ll <- ll_all[, keep, drop = FALSE]
  r_eff <- relative_eff(exp(ll), chain_id = rep(1:4, each = nrow(ll) / 4))
  loos_trial[[mm]] <- loo(ll, r_eff = r_eff)

  # log_lik per soggetto: somma sui trial (leave-one-subject-out)
  idx <- matrix(seq_len(NT), sdh$N, sdh$T)                 # log_lik[i,t] in ordine colonna
  lls <- sapply(seq_len(sdh$N), function(i)
    rowSums(ll_all[, idx[i, ], drop = FALSE]))
  r_eff_s <- relative_eff(exp(lls), chain_id = rep(1:4, each = nrow(lls) / 4))
  loos_subj[[mm]] <- loo(lls, r_eff = r_eff_s)

  rm(ll_all, ll, lls); gc()
  k <- loos_trial[[mm]]$diagnostics$pareto_k
  tab[[mm]] <- data.frame(
    model = mm,
    elpd_trial = loos_trial[[mm]]$estimates["elpd_loo", "Estimate"],
    se_trial   = loos_trial[[mm]]$estimates["elpd_loo", "SE"],
    elpd_per_trial = loos_trial[[mm]]$estimates["elpd_loo", "Estimate"] / NT,
    p_loo = loos_trial[[mm]]$estimates["p_loo", "Estimate"],
    k_gt_07 = sum(k > 0.7),
    elpd_subj = loos_subj[[mm]]$estimates["elpd_loo", "Estimate"],
    se_subj   = loos_subj[[mm]]$estimates["elpd_loo", "SE"])
}
tab <- do.call(rbind, tab)
tab$elpd_chance_per_trial <- log(1 / 4)
print(tab, digits = 4, row.names = FALSE)

cat("\n--- loo_compare, livello trial ---\n"); print(loo_compare(loos_trial))
cat("\n--- loo_compare, livello soggetto ---\n"); print(loo_compare(loos_subj))

write.csv(tab, "wcst_rl/loo_table.csv", row.names = FALSE)
saveRDS(list(trial = loos_trial, subj = loos_subj), "wcst_rl/fits/loos.RDS")
message("\nsalvato wcst_rl/loo_table.csv")
