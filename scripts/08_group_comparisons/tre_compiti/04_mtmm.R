# Quanta varianza individuale e' condivisa dalla stessa primitiva TRA compiti
# (tratto), quanta dal compito (metodo), quanta dalla coppia primitiva-entro-compito?
# Quattro strutture confrontate con LOO.
library(cmdstanr); library(loo)

# percorso di CmdStan: si usa CMDSTAN se impostato, altrimenti quello
# dell'ambiente conda; adattare qui se l'installazione e' altrove.
cs <- Sys.getenv("CMDSTAN")
if (!nzchar(cs)) cs <- file.path(Sys.getenv("CONDA_PREFIX"), "bin", "cmdstan")
if (dir.exists(cs)) set_cmdstan_path(cs)

d <- read.csv("dati/dati_armonizzati.csv")
# stime WCST del modello a una sola popolazione: quelle del modello con il
# predittore di gruppo sono restringite verso la media del proprio gruppo e non
# possono essere usate per analisi a livello individuale.
# log d e kappa del WCST sono esclusi: l'errore di misura per soggetto e' 2.1 e
# 10.2 volte la varianza tra soggetti (affidabilita' individuale implicata 0.32 e
# 0.08), quindi le loro stime per soggetto non portano differenze individuali.
# Restano nell'analisi degli effetti di gruppo, dove il rumore attenua e non gonfia.
PARS <- c("prl_a","ts_a_0","ts_a_1","prl_v","ts_v_0","ts_v_1",
          "prl_t","ts_t_0","ts_t_1","prl_alpha","prl_pos_alpha","wcst_logit_h")
SE   <- c("wcst_logit_h" = "wcst_logit_h_sd")
d <- d[complete.cases(d[, c(PARS, unname(SE))]), ]
Y <- scale(as.matrix(d[, PARS]))
g <- as.numeric(d$grp == "AN")
P <- length(PARS)
se <- matrix(0, nrow(Y), P)
for (nm in names(SE)) se[, which(PARS == nm)] <- d[[SE[[nm]]]] / attr(Y, "scaled:scale")[[nm]]

#            prl_a ts_a_0 ts_a_1 prl_v ts_v_0 ts_v_1 prl_t ts_t_0 ts_t_1 alpha pos_a  h
# h del WCST e' raggruppato con le learning rate del PRL (apprendimento e
# flessibilita'): e' l'unica coppia di tratto che coinvolge il WCST. Non ha
# carico di metodo perche' e' l'unico indicatore del proprio compito.
trait <- c(      1,     1,     1,    2,     2,     2,    3,     3,     3,    4,    4,  4)
meth  <- c(      1,     2,     2,    1,     2,     2,    1,     2,     2,    1,    1,  0)
# coppie stessa primitiva entro lo stesso compito (le due condizioni del TS,
# le due learning rate del PRL): senza questa componente la varianza di tratto
# assorbirebbe correlazioni che non sono affatto tra compiti
tmix <- c(      0,     1,     1,    0,     2,     2,    0,     3,     3,    4,    4,  0)
TR <- c("soglia","drift","tempo non dec.","apprendimento/flessibilita'")
ME <- c("PRL","task switching")
KN <- c("tratto (tra compiti)","metodo (compito)","primitiva entro compito")
cat("N =", nrow(Y), "| pazienti =", sum(g),
    "| errore di misura medio sugli indicatori WCST (varianza standardizzata):",
    paste(round(colMeans(se^2)[se[1, ] > 0], 3), collapse = " "), "\n")

G <- cbind(trait, meth, tmix)
base <- list(N = nrow(Y), P = P, K = 3, y = Y, g = g, se = se, grp = G)
S <- list("completo"           = c(base, list(use = as.array(c(1L, 1L, 1L)))),
          "senza tratto"       = c(base, list(use = as.array(c(0L, 1L, 1L)))),
          "senza metodo"       = c(base, list(use = as.array(c(1L, 0L, 1L)))),
          "solo entro compito" = c(base, list(use = as.array(c(0L, 0L, 1L)))),
          "indipendenza"       = c(base, list(use = as.array(c(0L, 0L, 0L)))))

mod <- cmdstan_model("tre_compiti/stan/mtmm.stan")
dir.create("risultati/fits", recursive = TRUE, showWarnings = FALSE)
fits <- list(); loos <- list()
for (nm in names(S)) {
  message("\n=== ", nm, " ===")
  f <- mod$sample(data = S[[nm]], chains = 4, parallel_chains = 4,
                  iter_warmup = 1000, iter_sampling = 1000, seed = 2026,
                  adapt_delta = 0.95, refresh = 0, show_messages = FALSE)
  ss <- f$summary(c("mu","b","psi"))
  cat("divergenze:", sum(f$diagnostic_summary()$num_divergent),
      "| max R-hat:", round(max(ss$rhat, na.rm = TRUE), 4),
      "| min ESS:", round(min(ss$ess_bulk, na.rm = TRUE)), "\n")
  fits[[nm]] <- f; loos[[nm]] <- f$loo(cores = 4)
  saveRDS(f, paste0("risultati/fits/mtmm_", gsub("[ +]", "_", nm), ".RDS"))
}
cmp <- loo_compare(loos)
print(cmp)
cmpdf <- as.data.frame(unclass(cmp))
# a seconda della versione di loo, il nome del modello sta nei rownames o in
# una colonna "model"; si prende quello che c'e'.
nomi <- if (!is.null(cmpdf$model)) as.character(cmpdf$model) else rownames(cmp)
cmpdf <- cbind(modello = nomi, cmpdf)
stopifnot(!any(is.na(cmpdf$modello)), all(cmpdf$modello %in% names(S)))
write.csv(cmpdf, "risultati/loo_mtmm.csv", row.names = FALSE)

# decomposizione della varianza nel modello completo
f <- fits[["completo"]]
V <- f$draws("var_comp", format = "draws_matrix")     # colonne in ordine [p, k]
vu <- f$draws("var_unica", format = "draws_matrix")
cn <- function(k) V[, paste0("var_comp[", 1:P, ",", k, "]"), drop = FALSE]
V1 <- cn(1); V2 <- cn(2); V3 <- cn(3)
dec <- data.frame(parametro = PARS,
                  primitiva = ifelse(trait == 0, "inferenza di regola", TR[pmax(trait, 1)]),
                  compito = ifelse(meth == 0, "WCST", ME[pmax(meth, 1)]),
                  tratto = colMeans(V1), tratto_lo = apply(V1, 2, quantile, .025),
                  tratto_hi = apply(V1, 2, quantile, .975),
                  metodo = colMeans(V2), entro_compito = colMeans(V3), unica = colMeans(vu))
write.csv(dec, "risultati/decomposizione_varianza.csv", row.names = FALSE)
print(round(dec[, c("tratto","tratto_lo","tratto_hi","metodo","entro_compito","unica")], 3))

# contrasto decisivo: varianza di tratto (tra compiti) vs metodo, per draw
mt <- rowMeans(V1); mm <- rowMeans(V2); m3 <- rowMeans(V3)
dif <- mt - mm
res <- data.frame(
  quantita = c("varianza di tratto (tra compiti)","varianza di metodo (compito)",
               "varianza primitiva entro compito","tratto - metodo"),
  media = c(mean(mt), mean(mm), mean(m3), mean(dif)),
  lo = c(quantile(mt, .025), quantile(mm, .025), quantile(m3, .025), quantile(dif, .025)),
  hi = c(quantile(mt, .975), quantile(mm, .975), quantile(m3, .975), quantile(dif, .975)),
  P_positivo = c(NA, NA, NA, mean(dif > 0)))
write.csv(res, "risultati/contrasto_tratto_metodo.csv", row.names = FALSE)
print(round(res[, -1], 4))

# per primitiva: quanta varianza e' davvero condivisa tra compiti
per_tr <- do.call(rbind, lapply(1:4, function(k) {
  idx <- which(trait == k); v <- rowMeans(V1[, idx, drop = FALSE])
  data.frame(primitiva = TR[k], n_ind = length(idx), var_tratto = mean(v),
             lo = quantile(v, .025), hi = quantile(v, .975), P_gt_05 = mean(v > 0.05))
}))
write.csv(per_tr, "risultati/varianza_tratto_per_primitiva.csv", row.names = FALSE)
print(round(per_tr[, -1], 3)); print(per_tr$primitiva)
saveRDS(list(fits = fits, loos = loos), "risultati/fits/mtmm_all.RDS")
message("\nsalvati loo_mtmm.csv, decomposizione_varianza.csv, contrasto_tratto_metodo.csv")
