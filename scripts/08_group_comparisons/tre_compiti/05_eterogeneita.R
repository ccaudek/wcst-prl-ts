# ---------------------------------------------------------------------------
# 05. Eterogeneita' nello spazio dei parametri.
#
# Due domande distinte:
#   (a) i pazienti occupano una regione PIU' AMPIA dello spazio dei parametri?
#       -> distanza media dal centroide del proprio gruppo, test di permutazione
#          sull'etichetta di gruppo (la statistica e' invariante allo scambio di
#          etichette sotto l'ipotesi nulla)
#   (b) i due gruppi hanno centroidi diversi?
#       -> distanza euclidea tra centroidi, stesso schema di permutazione
#
# Indicatori: i 12 affidabili (analisi principale), con due analisi di
# sensibilita' - senza l'indicatore WCST (solo PRL e task switching, i parametri
# che non abbiamo ristimato) e con l'aggiunta di wcst_log_d.
#
# Output  tre_compiti/dispersione_test.csv, tre_compiti/centroidi_test.csv,
#         tre_compiti/eterogeneita.csv, tre_compiti/dispersione_loo.csv
# ---------------------------------------------------------------------------

set.seed(20260905)
NPERM <- 10000

d <- read.csv("dati/dati_armonizzati.csv")
D <- d[d$completo, ]
IND12 <- c("prl_a","ts_a_0","ts_a_1","prl_v","ts_v_0","ts_v_1",
           "prl_t","ts_t_0","ts_t_1","prl_alpha","prl_pos_alpha","wcst_logit_h")
gy <- D$grp == "AN"

# dispersione = distanza media dal centroide del proprio gruppo
disp <- function(Z, g) {
  sapply(c(TRUE, FALSE), function(m) {
    A <- Z[g == m, , drop = FALSE]
    mean(sqrt(rowSums((A - matrix(colMeans(A), nrow(A), ncol(A), byrow = TRUE))^2)))
  })
}
cdist <- function(Z, g) sqrt(sum((colMeans(Z[g, ]) - colMeans(Z[!g, ]))^2))

perm_test <- function(Z, g, stat) {
  obs <- stat(Z, g)
  nul <- replicate(NPERM, stat(Z, sample(g)))
  list(obs = obs, p = mean(abs(nul) >= abs(obs)), nul = nul)
}

test_disp <- function(cols, etichetta) {
  Z <- scale(as.matrix(D[, cols]))
  dd <- disp(Z, gy)
  r  <- perm_test(Z, gy, function(Zm, g) diff(rev(disp(Zm, g))))
  data.frame(insieme = etichetta, n_indicatori = length(cols),
             dispersione_AN = dd[1], dispersione_HC = dd[2],
             differenza = dd[1] - dd[2], p_permutazione = r$p)
}

dt <- rbind(
  test_disp(IND12, "12 indicatori affidabili"),
  test_disp(setdiff(IND12, "wcst_logit_h"), "11, solo PRL e task switching"),
  test_disp(c(IND12, "wcst_log_d"), "13, con wcst_log_d"))
write.csv(dt, "risultati/dispersione_test.csv", row.names = FALSE)

Z12 <- scale(as.matrix(D[, IND12]))
rc <- perm_test(Z12, gy, cdist)
ct <- data.frame(quantita = "distanza tra centroidi (12 dim., z)",
                 valore = rc$obs, p_permutazione = rc$p)
write.csv(ct, "risultati/centroidi_test.csv", row.names = FALSE)

# --- rapporti di varianza per singolo parametro -------------------------------
# Brown-Forsythe (test di Levene centrato sulla mediana): robusto alla non
# normalita', al contrario del test F sul rapporto di varianze, che qui si usa
# solo per l'intervallo di confidenza.
bf <- function(x, g) {
  z <- abs(x - tapply(x, g, median)[as.character(g)])
  a <- anova(lm(z ~ factor(g)))
  c(W = a[["F value"]][1], p = a[["Pr(>F)"]][1])
}
het <- do.call(rbind, lapply(seq_along(IND12), function(j) {
  a <- Z12[gy, j]; h <- Z12[!gy, j]
  vr <- var(a) / var(h); b <- bf(Z12[, j], gy)
  data.frame(parametro = IND12[j], sd_AN = sd(a), sd_HC = sd(h), rapporto_var = vr,
             lo = vr / qf(.975, length(a) - 1, length(h) - 1),
             hi = vr / qf(.025, length(a) - 1, length(h) - 1),
             W_BrownForsythe = b[["W"]], p = b[["p"]])
}))
het$p_BH <- p.adjust(het$p, "BH")
write.csv(het, "risultati/eterogeneita.csv", row.names = FALSE)

# --- robustezza: togliere un indicatore alla volta ----------------------------
loo <- do.call(rbind, lapply(IND12, function(p) {
  Z <- scale(as.matrix(D[, setdiff(IND12, p)]))
  dd <- disp(Z, gy)
  data.frame(escluso = p, differenza = dd[1] - dd[2],
             p_permutazione = perm_test(Z, gy,
               function(Zm, g) diff(rev(disp(Zm, g))))$p)
}))
write.csv(loo, "risultati/dispersione_loo.csv", row.names = FALSE)

cat("\n== dispersione ==\n"); print(round(dt[, -1], 4)); print(dt$insieme)
cat("\n== centroidi ==\n"); print(round(ct[, -1], 4))
cat(sprintf("\nrapporti di varianza > 1: %d/%d | minimo p_BH = %.3f\n",
            sum(het$rapporto_var > 1), nrow(het), min(het$p_BH)))
cat(sprintf("robustezza: p perm resta < 0.05 in %d/%d esclusioni singole\n",
            sum(loo$p_permutazione < .05), nrow(loo)))
