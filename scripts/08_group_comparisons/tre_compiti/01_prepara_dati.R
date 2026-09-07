# ---------------------------------------------------------------------------
# 01. Armonizzazione dei parametri dei tre compiti in un unico data frame.
#
# Input   prl_params.csv, task_switching_params.csv  (stime fornite, usate come sono)
#         wcst_rl/wcst_params_nogroup.csv            (rifit HMM senza predittore di gruppo)
# Output  tre_compiti/dati_armonizzati.csv
#         tre_compiti/parametri_meta.csv
#
# Nota sulle stime del WCST: le stime individuali provengono dal modello adattato
# SENZA il predittore di gruppo. Le stime del modello con il predittore sono
# contratte verso la media del proprio gruppo e non possono essere usate per
# analisi a livello individuale (vedi wcst_rl/10_nogroup.R).
# ---------------------------------------------------------------------------

prl <- read.csv("dati/prl_params.csv")
ts  <- read.csv("dati/task_switching_params.csv")
wc  <- read.csv("dati/wcst_params_nogroup.csv")

names(prl)[names(prl) != "user_id"] <- paste0("prl_", setdiff(names(prl), "user_id"))
names(ts)[names(ts)  != "user_id"] <- paste0("ts_",  setdiff(names(ts),  "user_id"))

# I due file usano vocabolari diversi per la stessa etichetta
# (patients/controls nel PRL, AN/HC nel task switching): si normalizza, si
# verifica che concordino dove il soggetto e' presente in entrambi, e si prende
# la prima disponibile.
norm_grp <- function(x) ifelse(is.na(x), NA,
  ifelse(as.character(x) %in% c("patients", "AN", "an"), "AN", "HC"))
d <- merge(prl, ts, by = "user_id", all = TRUE)
d$prl_group <- norm_grp(d$prl_group); d$ts_group <- norm_grp(d$ts_group)
ent <- !is.na(d$prl_group) & !is.na(d$ts_group)
stopifnot(all(as.character(d$prl_group[ent]) == as.character(d$ts_group[ent])))
d$grp <- ifelse(is.na(d$prl_group), as.character(d$ts_group), as.character(d$prl_group))
d$prl_group <- NULL; d$ts_group <- NULL
d <- merge(d, wc, by = "user_id", all = TRUE)

# i tassi di apprendimento del PRL sono forniti su scala logit: si aggiunge la
# versione su scala di probabilita' per le tabelle descrittive
d$prl_alpha_prob     <- plogis(d$prl_alpha)
d$prl_pos_alpha_prob <- plogis(d$prl_pos_alpha)

# i 12 indicatori usati nelle analisi di differenze individuali
IND12 <- c("prl_a","ts_a_0","ts_a_1","prl_v","ts_v_0","ts_v_1",
           "prl_t","ts_t_0","ts_t_1","prl_alpha","prl_pos_alpha","wcst_logit_h")
d$completo <- complete.cases(d[, c(IND12, "wcst_log_d", "wcst_kappa")]) & !is.na(d$grp)

meta <- data.frame(
  parametro = c("prl_a","ts_a_0","ts_a_1","prl_v","ts_v_0","ts_v_1",
                "prl_t","ts_t_0","ts_t_1","prl_alpha","prl_pos_alpha",
                "wcst_logit_h","wcst_log_d","wcst_kappa"),
  primitiva = c(rep("soglia", 3), rep("drift", 3), rep("tempo non dec.", 3),
                rep("apprendimento", 2), "inferenza di regola",
                "soglia/decisione", "perseverazione"),
  compito   = c("PRL","TS ripetizione","TS switch","PRL","TS ripetizione","TS switch",
                "PRL","TS ripetizione","TS switch","PRL","PRL","WCST","WCST","WCST"),
  affidabile_per_diff_individuali =
              c(rep(TRUE, 12), FALSE, FALSE))

dir.create("tre_compiti", showWarnings = FALSE)
write.csv(d, "dati/dati_armonizzati.csv", row.names = FALSE)
write.csv(meta, "dati/parametri_meta.csv", row.names = FALSE)

cat(sprintf("N totale = %d | completi nei tre compiti = %d (pazienti %d, controlli %d)\n",
            nrow(d), sum(d$completo), sum(d$completo & d$grp == "AN"),
            sum(d$completo & d$grp != "AN")))
