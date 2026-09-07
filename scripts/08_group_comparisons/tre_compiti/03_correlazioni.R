# ---------------------------------------------------------------------------
# 03. Struttura di correlazione tra i 12 indicatori affidabili, classificata per
#     tipo di coppia, e affidabilita' individuale implicata dei parametri WCST.
#
# La domanda: le coppie che misurano la STESSA primitiva in compiti DIVERSI sono
# piu' correlate delle coppie che non condividono nulla? Se le primitive fossero
# tratti individuali, si', e di molto.
#
# Output  tre_compiti/correlazioni_spearman.csv, tre_compiti/coppie_correlazioni.csv,
#         tre_compiti/struttura_cross_task.csv, tre_compiti/affidabilita_wcst.csv
# ---------------------------------------------------------------------------

d  <- read.csv("dati/dati_armonizzati.csv")
D  <- d[d$completo, ]

# --- affidabilita' individuale implicata dei parametri WCST -------------------
# Per ogni parametro, la varianza tra soggetti delle medie a posteriori e
# l'errore quadratico medio a posteriori. Il rapporto
#   var_tra / (var_tra + errore^2)
# e' la quota di varianza osservata attribuibile a differenze vere: sotto 0.5 la
# stima per soggetto e' dominata dall'incertezza e non e' usabile per
# differenze individuali.
WP <- c("wcst_logit_h", "wcst_log_d", "wcst_kappa")
rel <- do.call(rbind, lapply(WP, function(p) {
  v <- var(D[[p]]); e <- mean(D[[paste0(p, "_sd")]]^2)
  data.frame(parametro = p, var_tra_soggetti = v, errore_quad_medio = e,
             affidabilita_individuale = v / (v + e),
             usato_per_differenze_individuali = v / (v + e) > .5)
}))
write.csv(rel, "risultati/affidabilita_wcst.csv", row.names = FALSE)

# --- correlazioni sui 12 indicatori affidabili -------------------------------
ORD <- c("prl_a","ts_a_0","ts_a_1","prl_v","ts_v_0","ts_v_1",
         "prl_t","ts_t_0","ts_t_1","prl_alpha","prl_pos_alpha","wcst_logit_h")
PRIM <- c(prl_a = "soglia", ts_a_0 = "soglia", ts_a_1 = "soglia",
          prl_v = "drift", ts_v_0 = "drift", ts_v_1 = "drift",
          prl_t = "tempo", ts_t_0 = "tempo", ts_t_1 = "tempo",
          prl_alpha = "appr", prl_pos_alpha = "appr", wcst_logit_h = "inferenza")
CMP  <- c(prl_a = "PRL", ts_a_0 = "TS", ts_a_1 = "TS",
          prl_v = "PRL", ts_v_0 = "TS", ts_v_1 = "TS",
          prl_t = "PRL", ts_t_0 = "TS", ts_t_1 = "TS",
          prl_alpha = "PRL", prl_pos_alpha = "PRL", wcst_logit_h = "WCST")

R <- cor(D[, ORD], method = "spearman")
write.csv(R, "risultati/correlazioni_spearman.csv")

cp <- do.call(rbind, lapply(1:(length(ORD) - 1), function(i)
  do.call(rbind, lapply((i + 1):length(ORD), function(j) {
    a <- ORD[i]; b <- ORD[j]
    sp <- PRIM[[a]] == PRIM[[b]]; st <- CMP[[a]] == CMP[[b]]
    data.frame(par1 = a, par2 = b, rho = R[a, b], classe =
      if (sp && st) "stessa primitiva, entro compito"
      else if (sp)   "stessa primitiva, tra compiti"
      else if (st)   "primitive diverse, entro compito"
      else           "primitive diverse, tra compiti")
  }))))
write.csv(cp, "risultati/coppie_correlazioni.csv", row.names = FALSE)

tab <- do.call(rbind, lapply(split(cp, cp$classe), function(x)
  data.frame(classe = x$classe[1], n = nrow(x), media_abs = mean(abs(x$rho)),
             max_abs = max(abs(x$rho)), media = mean(x$rho))))
tab <- tab[order(-tab$media_abs), ]
write.csv(tab, "risultati/struttura_cross_task.csv", row.names = FALSE)

# --- correlazioni di h del WCST, con disattenuazione -------------------------
rh <- rel$affidabilita_individuale[rel$parametro == "wcst_logit_h"]
hc <- data.frame(altro = setdiff(ORD, "wcst_logit_h"))
hc$rho <- R["wcst_logit_h", hc$altro]
hc$rho_disattenuato <- hc$rho / sqrt(rh)
write.csv(hc, "risultati/correlazioni_h_wcst.csv", row.names = FALSE)

cat("\n== affidabilita' individuale, parametri WCST ==\n"); print(round(rel[, -1], 3))
cat("\n== correlazioni per classe di coppia ==\n"); print(round(tab[, -1], 3)); print(tab$classe)
cat(sprintf("\ncoppie stessa primitiva tra compiti: n = %d, media |rho| = %.3f, max |rho| = %.3f\n",
            sum(cp$classe == "stessa primitiva, tra compiti"),
            mean(abs(cp$rho[cp$classe == "stessa primitiva, tra compiti"])),
            max(abs(cp$rho[cp$classe == "stessa primitiva, tra compiti"]))))
