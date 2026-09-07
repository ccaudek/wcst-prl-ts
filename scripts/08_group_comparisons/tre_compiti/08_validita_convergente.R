# 08_validita_convergente.R -------------------------------------------------
# Validita' convergente del parametro WCST logit h contro gli indici
# comportamentali model-free, sulle stime SENZA gruppo.
#
# Perche' rifarla qui: in wcst_hmm/09_external_sticky.R la stessa verifica e'
# fatta sulle stime CON gruppo, dove la correlazione potrebbe essere gonfiata
# dal fatto che la differenza di gruppo e' presente in entrambi i termini.
# Qui usiamo le stime senza gruppo e riportiamo anche la correlazione
# calcolata entro gruppo (centrando entrambe le variabili sul proprio gruppo),
# che e' la forma non circolare della verifica.
#
# Output: risultati/validita_convergente_wcst.csv

set.seed(20250906)

bh <- read.csv("dati/wcst_behav_indices.csv")
ng <- read.csv("dati/wcst_params_nogroup.csv")

m <- merge(bh, ng, by.x = "subj_name", by.y = "user_id")
m$g <- ifelse(tolower(m$group) %in% c("an", "pazienti", "patients"), 1L, 0L)
stopifnot(all(m$g %in% c(0L, 1L)), nrow(m) > 50)

IDX <- c("prop_pers_err", "prop_non_pers_err", "prop_pers_resp")

centra <- function(x, g) x - ave(x, g, FUN = function(z) mean(z, na.rm = TRUE))

out <- do.call(rbind, lapply(IDX, function(i) {
  a  <- suppressWarnings(cor.test(m$wcst_logit_h, m[[i]], method = "spearman"))
  b  <- suppressWarnings(cor.test(centra(m$wcst_logit_h, m$g),
                                  centra(m[[i]], m$g), method = "spearman"))
  data.frame(indice = i, n = nrow(m),
             rho = unname(a$estimate), p = a$p.value,
             rho_entro_gruppo = unname(b$estimate), p_entro = b$p.value)
}))

write.csv(round_df <- data.frame(indice = out$indice, n = out$n,
                                 rho = round(out$rho, 3), p = signif(out$p, 3),
                                 rho_entro_gruppo = round(out$rho_entro_gruppo, 3),
                                 p_entro = signif(out$p_entro, 3)),
          "risultati/validita_convergente_wcst.csv", row.names = FALSE)

cat("== validita' convergente di logit h (stime senza gruppo) ==\n")
print(round_df)
cat(sprintf("| n = %d (%d pazienti, %d controlli)\n",
            nrow(m), sum(m$g == 1), sum(m$g == 0)))
