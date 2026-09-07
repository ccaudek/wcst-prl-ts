# Confronto tra strutture di covarianza: primitive di processo vs compito vs
# fattore generale vs indipendenza. Le medie (mu_p + b_p * gruppo) sono
# identiche in tutti i modelli, quindi LOO confronta solo la struttura latente.
library(cmdstanr); library(loo)

d <- read.csv("tre_compiti/dati_armonizzati.csv")
PARS <- c("prl_a","ts_a_0","ts_a_1","prl_v","ts_v_0","ts_v_1",
          "prl_t","ts_t_0","ts_t_1","prl_alpha","prl_pos_alpha","wcst_logit_h")
d <- d[complete.cases(d[, c(PARS, "wcst_logit_h_sd")]), ]
Y <- scale(as.matrix(d[, PARS]))
sd_h <- attr(Y, "scaled:scale")[["wcst_logit_h"]]
se_w <- d$wcst_logit_h_sd / sd_h
g <- as.numeric(d$grp == "AN")
N <- nrow(Y); P <- length(PARS); wcol <- which(PARS == "wcst_logit_h")
rel_h <- 1 - mean(se_w^2)
cat("N =", N, "| P =", P, "| pazienti =", sum(g), "| affidabilita' implicita di h =",
    round(rel_h, 3), "\n")

base <- list(N = N, P = P, y = Y, g = g, se_w = se_w, wcol = wcol)
mk <- function(F, anch, amap, fidx, fmap, fxidx = integer(0), fxmap = integer(0),
               fxval = numeric(0), res_free = rep(1L, P)) {
  c(base, list(F = F, Fa = length(anch), anch = as.array(anch), amap = as.array(amap),
               Nf = length(fidx), fidx = as.array(fidx), fmap = as.array(fmap),
               Nfx = length(fxidx), fxidx = as.array(fxidx), fxmap = as.array(fxmap),
               fxval = as.array(fxval), res_free = as.array(res_free)))
}
rf0 <- rep(1L, P); rf0[wcol] <- 0L   # per il fattore a indicatore singolo

M <- list(
  # 4 primitive di processo: soglia, drift, tempo non decisionale, apprendimento
  processo = mk(4, anch = c(1, 4, 7, 10), amap = c(1, 2, 3, 4),
                fidx = c(2, 3, 5, 6, 8, 9, 11, 12), fmap = c(1, 1, 2, 2, 3, 3, 4, 4)),
  # 3 fattori di compito; il WCST ha un solo indicatore: loading fissato a
  # sqrt(affidabilita') e residuo pari al solo errore di misura noto
  compito  = mk(3, anch = c(1, 2), amap = c(1, 2),
                fidx = c(4, 7, 10, 11, 3, 5, 6, 8, 9), fmap = c(1, 1, 1, 1, 2, 2, 2, 2, 2),
                fxidx = wcol, fxmap = 3, fxval = sqrt(rel_h), res_free = rf0),
  # fattore generale unico
  generale = mk(1, anch = 1, amap = 1, fidx = 2:P, fmap = rep(1, P - 1)),
  # indipendenza: tutti i loading nulli
  indip    = mk(1, anch = integer(0), amap = integer(0), fidx = integer(0), fmap = integer(0),
                fxidx = 1:P, fxmap = rep(1, P), fxval = rep(0, P))
)

mod <- cmdstan_model("tre_compiti/stan/cfa_marg.stan")
fits <- list(); loos <- list()
for (nm in names(M)) {
  message("\n=== ", nm, " ===")
  f <- mod$sample(data = M[[nm]], chains = 4, parallel_chains = 4,
                  iter_warmup = 1000, iter_sampling = 1000, seed = 2026,
                  adapt_delta = 0.95, refresh = 0, show_messages = FALSE)
  s <- f$summary(c("mu","b","lam_a","lam_f","psi_raw"))
  cat("divergenze:", sum(f$diagnostic_summary()$num_divergent),
      "| max R-hat:", round(max(s$rhat, na.rm = TRUE), 4),
      "| min ESS:", round(min(s$ess_bulk, na.rm = TRUE)), "\n")
  fits[[nm]] <- f
  loos[[nm]] <- f$loo(cores = 4)
  saveRDS(f, paste0("tre_compiti/fits/cfa_", nm, ".RDS"))
}
cmp <- loo_compare(loos)
print(cmp)
write.csv(data.frame(modello = rownames(cmp), cmp), "tre_compiti/loo_cfa.csv", row.names = FALSE)
saveRDS(loos, "tre_compiti/fits/loo_cfa.RDS")

# loading e correlazioni tra fattori del modello di processo
fp <- fits$processo
Ld <- fp$draws("Lam", format = "draws_matrix")
lab <- expand.grid(p = 1:P, f = 1:4)
keep <- which(colSums(abs(Ld)) > 0)
out <- data.frame(parametro = PARS[lab$p[keep]], fattore = lab$f[keep],
                  media = colMeans(Ld[, keep]),
                  lo = apply(Ld[, keep], 2, quantile, 0.025),
                  hi = apply(Ld[, keep], 2, quantile, 0.975))
out$primitiva <- c("soglia","drift","tempo non dec.","apprendimento")[out$fattore]
write.csv(out, "tre_compiti/loadings_processo.csv", row.names = FALSE)
print(round(out[, c("fattore","media","lo","hi")], 3))

Phi <- fp$draws("Phi", format = "draws_matrix")
pm <- matrix(colMeans(Phi), 4, 4)
dimnames(pm) <- list(c("soglia","drift","t_nondec","apprend"), c("soglia","drift","t_nondec","apprend"))
write.csv(round(pm, 3), "tre_compiti/phi_processo.csv")
print(round(pm, 3))

# quota di varianza di ciascun indicatore spiegata dal proprio fattore (comunalita')
com <- out
com$comunalita <- com$media^2
write.csv(com, "tre_compiti/comunalita_processo.csv", row.names = FALSE)
cat("\ncomunalita' mediana:", round(median(com$comunalita), 3),
    "| massima:", round(max(com$comunalita), 3), "\n")
message("\nsalvati loo_cfa.csv, loadings_processo.csv, phi_processo.csv")
