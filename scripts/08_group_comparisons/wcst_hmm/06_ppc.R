# 06_ppc.R -----------------------------------------------------------------
# Posterior predictive check sulle firme comportamentali classiche del WCST.
# Non basta che il modello vinca il LOO: deve riprodurre le quantita' che i
# clinici leggono (accuratezza, errori perseverativi e non, win-stay/lose-shift,
# recupero dopo il cambio di regola).
#
# Si simula in avanti dai draw a posteriori: il feedback e' quello vero del
# compito, quindi la simulazione e' un test genuino del processo generativo.

suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })
source("wcst_rl/funs_simulate.R")

sdh <- readRDS("wcst_rl/wcst_stan_data_hier.RDS")
fit <- readRDS("wcst_rl/fits/hmm_hier.RDS")
N <- sdh$N

obs <- wcst_signatures(sdh, sdh$choice, sdh$rew)
SIG <- c("acc", "prop_pers_err", "prop_non_pers_err", "prop_rep_dim",
         "win_stay", "lose_shift", "acc_pos1_3", "acc_pos8_10")

dh <- fit$draws("h", format = "draws_matrix"); dd <- fit$draws("d", format = "draws_matrix")
dl <- as.vector(fit$draws("lapse", format = "draws_matrix"))
de <- as.vector(fit$draws("eta", format = "draws_matrix"))
NREP <- 200
set.seed(11)
idx <- sample(nrow(dh), NREP)

rep_group <- matrix(NA_real_, NREP, length(SIG), dimnames = list(NULL, SIG))
rep_subj  <- array(NA_real_, c(NREP, N, length(SIG)))
for (k in seq_along(idx)) {
  j <- idx[k]
  s <- sim_hmm(sdh, dh[j, ], dd[j, ], rep(dl[j], N), rep(de[j], N), seed = 1000 + k)
  sg <- wcst_signatures(sdh, s$choice, s$rew)
  rep_group[k, ] <- colMeans(sg[, SIG], na.rm = TRUE)
  rep_subj[k, , ] <- as.matrix(sg[, SIG])
  if (k %% 50 == 0) message("  replicazione ", k, "/", NREP)
}

ppc <- data.frame(
  indice = SIG,
  osservato = sapply(SIG, function(v) mean(obs[[v]], na.rm = TRUE)),
  predetto  = colMeans(rep_group, na.rm = TRUE),
  q05 = apply(rep_group, 2, quantile, .05, na.rm = TRUE),
  q95 = apply(rep_group, 2, quantile, .95, na.rm = TRUE))
ppc$p_bayes <- sapply(seq_along(SIG), function(k)
  mean(rep_group[, k] >= mean(obs[[SIG[k]]], na.rm = TRUE), na.rm = TRUE))
ppc$dentro_90 <- ppc$osservato >= ppc$q05 & ppc$osservato <= ppc$q95
cat("--- PPC a livello di gruppo ---\n"); print(ppc, digits = 3, row.names = FALSE)

# correlazione tra indice osservato e predetto a livello individuale:
# misura se il modello cattura le differenze individuali, non solo la media
cor_ind <- sapply(seq_along(SIG), function(k) {
  pred <- colMeans(rep_subj[, , k], na.rm = TRUE)
  suppressWarnings(cor(obs[[SIG[k]]], pred, use = "complete.obs"))
})
ind <- data.frame(indice = SIG, r_individuale = cor_ind)
cat("\n--- PPC a livello individuale (r osservato-predetto) ---\n")
print(ind, digits = 3, row.names = FALSE)

write.csv(ppc, "wcst_rl/ppc_table.csv", row.names = FALSE)
saveRDS(list(obs = obs, rep_group = rep_group, rep_subj = rep_subj,
             SIG = SIG, ppc = ppc, ind = ind), "wcst_rl/fits/ppc.RDS")
message("\nsalvato wcst_rl/ppc_table.csv")
