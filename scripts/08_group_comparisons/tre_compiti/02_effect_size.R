# ---------------------------------------------------------------------------
# 02. Effetti di gruppo per parametro: d di Cohen con IC bootstrap, e test di
#     concordanza di segno dentro ogni primitiva.
#
# Tutti e 14 i parametri entrano in questa analisi, compresi i due parametri
# WCST con bassa affidabilita' individuale (log d e kappa): il rumore delle
# stime per soggetto e' centrato, quindi ATTENUA il contrasto tra gruppi invece
# di gonfiarlo, e i d riportati sono conservativi.
#
# Output  tre_compiti/effect_sizes.csv, tre_compiti/concordanza_segni.csv
# ---------------------------------------------------------------------------

set.seed(20260905)
B <- 10000

d  <- read.csv("dati/dati_armonizzati.csv")
mt <- read.csv("dati/parametri_meta.csv")
PARS <- mt$parametro

cohen_d <- function(a, h) {
  sp <- sqrt(((length(a) - 1) * var(a) + (length(h) - 1) * var(h)) /
             (length(a) + length(h) - 2))
  (mean(a) - mean(h)) / sp
}

boot_d <- function(a, h, B = 10000) {
  bs <- replicate(B, cohen_d(sample(a, replace = TRUE), sample(h, replace = TRUE)))
  c(d = cohen_d(a, h), lo = quantile(bs, .025, names = FALSE),
    hi = quantile(bs, .975, names = FALSE),
    p_boot = 2 * min(mean(bs > 0), mean(bs < 0)))
}

# due campioni: "massimo" = tutti i casi disponibili per quel parametro,
# "completo" = i soli soggetti con dati nei tre compiti (n = 80)
out <- list()
for (p in PARS) {
  for (camp in c("massimo", "completo")) {
    sub <- if (camp == "completo") d[d$completo, ] else d
    ok  <- !is.na(sub[[p]]) & !is.na(sub$grp)
    a <- sub[[p]][ok & sub$grp == "AN"]
    h <- sub[[p]][ok & sub$grp != "AN"]
    r <- boot_d(a, h, B)
    out[[length(out) + 1]] <- data.frame(
      parametro = p, campione = camp, n_AN = length(a), n_HC = length(h),
      media_AN = mean(a), media_HC = mean(h), sd_AN = sd(a), sd_HC = sd(h),
      d = r[["d"]], lo = r[["lo"]], hi = r[["hi"]], p_boot = r[["p_boot"]])
  }
}
es <- do.call(rbind, out)
write.csv(es, "risultati/effect_sizes.csv", row.names = FALSE)

# --- concordanza di segno dentro ogni primitiva --------------------------------
# Se le primitive sono compromesse in modo coerente, dentro una primitiva tutti
# gli indicatori dei diversi compiti devono avere lo stesso segno. Il test
# principale usa gli 11 parametri con parametrizzazione DDM (PRL e task
# switching): il determinismo dell'HMM del WCST non ha la stessa orientazione di
# scala della soglia DDM, quindi entra solo nella versione di sensibilita'.
esc <- es[es$campione == "completo", ]
rownames(esc) <- esc$parametro
GRUPPI <- list(soglia        = c("prl_a","ts_a_0","ts_a_1","wcst_log_d"),
               drift         = c("prl_v","ts_v_0","ts_v_1"),
               tempo         = c("prl_t","ts_t_0","ts_t_1"),
               apprendimento = c("prl_alpha","prl_pos_alpha"))

conta <- function(gruppi) {
  cc <- tt <- 0
  for (g in gruppi) {
    s  <- sign(esc[g, "d"])
    mo <- sign(sum(s)); if (mo == 0) mo <- 1
    cc <- cc + sum(s == mo); tt <- tt + length(s)
  }
  c(concordi = cc, totale = tt,
    p = binom.test(cc, tt, .5, alternative = "greater")$p.value)
}
solo_ddm <- lapply(GRUPPI, function(g) g[!grepl("^wcst", g)])
cs <- rbind(
  data.frame(insieme = "11 parametri DDM (PRL + TS)", t(conta(solo_ddm))),
  data.frame(insieme = "12 con wcst_log_d",           t(conta(GRUPPI))))
write.csv(cs, "risultati/concordanza_segni.csv", row.names = FALSE)

cat("\n== effetti maggiori (campione completo) ==\n")
o <- esc[order(-abs(esc$d)), c("parametro","d","lo","hi")]
print(round(o[1:8, -1], 3)); print(o$parametro[1:8])
cat("\n== concordanza di segno ==\n"); print(cs)
