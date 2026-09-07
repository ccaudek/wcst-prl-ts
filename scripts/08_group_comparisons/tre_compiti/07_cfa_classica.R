# Complemento con la CFA classica (ML, indici di adattamento convenzionali):
# la struttura "per primitiva" contro la struttura "per compito".
library(lavaan)

d <- read.csv("dati/dati_armonizzati.csv")
PARS <- c("prl_a","ts_a_0","ts_a_1","prl_v","ts_v_0","ts_v_1",
          "prl_t","ts_t_0","ts_t_1","prl_alpha","prl_pos_alpha","wcst_logit_h")
D <- d[complete.cases(d[, PARS]), ]
D[PARS] <- scale(D[PARS])

M <- list(
  "per primitiva" = '
    soglia =~ prl_a + ts_a_0 + ts_a_1
    drift  =~ prl_v + ts_v_0 + ts_v_1
    tnd    =~ prl_t + ts_t_0 + ts_t_1
    appr   =~ prl_alpha + prl_pos_alpha + wcst_logit_h',
  "per compito" = '
    PRL =~ prl_a + prl_v + prl_t + prl_alpha + prl_pos_alpha
    TS  =~ ts_a_0 + ts_a_1 + ts_v_0 + ts_v_1 + ts_t_0 + ts_t_1',
  "fattore unico" = '
    gen =~ prl_a + ts_a_0 + ts_a_1 + prl_v + ts_v_0 + ts_v_1 +
           prl_t + ts_t_0 + ts_t_1 + prl_alpha + prl_pos_alpha + wcst_logit_h')

IX <- c("chisq","df","pvalue","cfi","tli","rmsea","rmsea.ci.lower","rmsea.ci.upper","srmr","aic","bic")
out <- list(); loads <- list()
for (nm in names(M)) {
  f <- try(cfa(M[[nm]], data = D, std.lv = TRUE, estimator = "ML"), silent = TRUE)
  if (inherits(f, "try-error")) { message(nm, ": errore"); next }
  ok <- lavInspect(f, "converged")
  if (!ok) {
    cat("\n---", nm, "| convergenza: FALSE (nessun indice di adattamento disponibile)\n")
    out[[nm]] <- data.frame(modello = nm, convergenza = FALSE,
                            t(as.data.frame(setNames(rep(NA_real_, length(IX)), IX))))
    next
  }
  ft <- fitMeasures(f, IX)
  out[[nm]] <- data.frame(modello = nm, convergenza = ok, t(as.data.frame(ft)))
  s <- standardizedSolution(f)
  s <- s[s$op == "=~", c("lhs","rhs","est.std","se","pvalue")]
  loads[[nm]] <- data.frame(modello = nm, s)
  cat("\n---", nm, "| convergenza:", ok, "\n")
  print(round(ft, 3))
}
fitab <- do.call(rbind, out); rownames(fitab) <- NULL
names(fitab) <- gsub("^X", "", names(fitab))
write.csv(fitab, "risultati/cfa_classica_fit.csv", row.names = FALSE)
ltab <- do.call(rbind, loads); rownames(ltab) <- NULL
write.csv(ltab, "risultati/cfa_classica_loadings.csv", row.names = FALSE)

cat("\n== range dei loading per modello convergente ==\n")
for (nm in names(loads)) cat(sprintf("%-16s min %.3f  max %.3f  mediana |.| %.3f\n", nm,
  min(loads[[nm]]$est.std), max(loads[[nm]]$est.std), median(abs(loads[[nm]]$est.std))))
message("\nsalvati cfa_classica_fit.csv, cfa_classica_loadings.csv")
