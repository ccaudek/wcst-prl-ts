# ---------------------------------------------------------------------------
# 06. Quanto valgono i tre compiti, insieme e separatamente, per distinguere i
#     due gruppi: AUC in validazione incrociata NIDIFICATA.
#
# Nidificata significa che la selezione del modello (quanti predittori tenere e
# quanta penalizzazione applicare) avviene DENTRO ogni fold di addestramento e
# non vede mai i dati di test di quel fold. Si riporta anche la versione non
# nidificata, per quantificare l'ottimismo che essa introduce: e' la quantita'
# che rende non confrontabili molte AUC pubblicate.
#
# Pipeline dentro ogni fold di addestramento:
#   standardizzazione -> selezione univariata dei k predittori con |t| maggiore
#   -> regressione logistica con penalizzazione ridge (glmnet, alpha = 0)
# Griglia: k in {2, 3, 5, tutti}, lambda sul percorso automatico di glmnet.
# Esterno: 5 fold stratificati x 20 ripetizioni = 100 fold.
#
# Qui entrano tutti e 14 i parametri. I due parametri WCST con bassa
# affidabilita' individuale sono esclusi dalle analisi di struttura (script 03 e
# 04) perche' li' l'attenuazione distorce le correlazioni, ma per la previsione
# il rumore non e' un problema di validita': un predittore rumoroso predice
# semplicemente peggio, e la validazione incrociata lo misura.
#
# Output  tre_compiti/auc_nested.csv, tre_compiti/auc_confronti.csv,
#         tre_compiti/auc_fold.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(glmnet))
set.seed(20260905)
NREP <- 20; NFOLD <- 5
KGRID <- c(2, 3, 5, Inf)

d <- read.csv("dati/dati_armonizzati.csv")
D <- d[d$completo, ]
y <- as.integer(D$grp == "AN")

P_PRL  <- c("prl_a","prl_v","prl_t","prl_alpha","prl_pos_alpha")
P_TS   <- c("ts_a_0","ts_a_1","ts_v_0","ts_v_1","ts_t_0","ts_t_1")
P_WCST <- c("wcst_logit_h","wcst_log_d","wcst_kappa")
SETS <- list("PRL" = P_PRL, "task switching" = P_TS, "WCST" = P_WCST,
             "PRL + TS" = c(P_PRL, P_TS), "PRL + WCST" = c(P_PRL, P_WCST),
             "TS + WCST" = c(P_TS, P_WCST),
             "tutti e tre" = c(P_PRL, P_TS, P_WCST))

auc <- function(y, s) {
  r <- rank(s); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# fold stratificati: si mescola dentro ogni classe, cosi' la proporzione di
# pazienti e' la stessa in ogni fold
strat_folds <- function(y, k) {
  f <- integer(length(y))
  for (cl in unique(y)) {
    i <- which(y == cl)
    f[sample(i)] <- rep_len(1:k, length(i))
  }
  f
}

# addestra la pipeline con k fissato; restituisce una funzione di punteggio
fit_pipe <- function(Xtr, ytr, k, lambda = NULL) {
  mu <- colMeans(Xtr); sg <- apply(Xtr, 2, sd); sg[sg == 0] <- 1
  Ztr <- sweep(sweep(Xtr, 2, mu), 2, sg, "/")
  tst <- abs(apply(Ztr, 2, function(x) {
    tt <- try(t.test(x[ytr == 1], x[ytr == 0])$statistic, silent = TRUE)
    if (inherits(tt, "try-error")) 0 else tt
  }))
  kk  <- min(k, ncol(Ztr))
  sel <- order(-tst)[1:kk]
  # con un solo predittore glmnet richiede almeno due colonne: si duplica
  Zs <- Ztr[, sel, drop = FALSE]
  if (ncol(Zs) == 1) Zs <- cbind(Zs, 0)
  g <- glmnet(Zs, ytr, family = "binomial", alpha = 0, lambda = lambda)
  list(mu = mu, sg = sg, sel = sel, g = g, dup = ncol(Zs) > kk)
}
score_pipe <- function(f, Xte, s = NULL) {
  Zte <- sweep(sweep(Xte, 2, f$mu), 2, f$sg, "/")[, f$sel, drop = FALSE]
  if (f$dup) Zte <- cbind(Zte, 0)
  predict(f$g, Zte, s = s, type = "link")
}

# lambda comune a tutti i fold, cosi' i percorsi sono allineati
LAM <- exp(seq(log(50), log(0.01), length.out = 30))

# I fold esterni si generano UNA volta e si riusano per ogni insieme: il
# confronto tra insiemi diventa appaiato fold per fold.
FOLDS <- lapply(1:NREP, function(i) strat_folds(y, NFOLD))

fold_rows <- list()
for (nm in names(SETS)) {
  X <- as.matrix(D[, SETS[[nm]]])
  for (rep in 1:NREP) {
    fo <- FOLDS[[rep]]
    for (v in 1:NFOLD) {
      tr <- fo != v; te <- !tr
      Xtr <- X[tr, , drop = FALSE]; ytr <- y[tr]

      # --- interno: scelta di (k, lambda) sui soli dati di addestramento ------
      fi <- strat_folds(ytr, NFOLD)
      sc <- array(NA_real_, c(length(KGRID), length(LAM)))
      for (ki in seq_along(KGRID)) {
        pr <- matrix(NA_real_, length(ytr), length(LAM))
        for (w in 1:NFOLD) {
          it <- fi != w
          if (length(unique(ytr[it])) < 2 || sum(!it) == 0) next
          f <- fit_pipe(Xtr[it, , drop = FALSE], ytr[it], KGRID[ki], LAM)
          pr[!it, ] <- score_pipe(f, Xtr[!it, , drop = FALSE])
        }
        ok <- !is.na(pr[, 1])
        sc[ki, ] <- apply(pr[ok, , drop = FALSE], 2, function(s) auc(ytr[ok], s))
      }
      bi <- which(sc == max(sc, na.rm = TRUE), arr.ind = TRUE)[1, ]
      kb <- KGRID[bi[1]]; lb <- LAM[bi[2]]

      # --- esterno: si valuta la configurazione scelta sul fold tenuto fuori --
      f <- fit_pipe(Xtr, ytr, kb, LAM)
      a_nid <- auc(y[te], as.numeric(score_pipe(f, X[te, , drop = FALSE], s = lb)))
      fold_rows[[length(fold_rows) + 1]] <-
        data.frame(insieme = nm, rip = rep, fold = v, auc = a_nid, k = kb, lambda = lb)
    }
  }
  cat("fatto:", nm, "\n")
}
fold <- do.call(rbind, fold_rows)

# --- versione NON nidificata: (k, lambda) scelti una volta su tutti i dati ----
flat_rows <- list()
for (nm in names(SETS)) {
  X <- as.matrix(D[, SETS[[nm]]])
  fi <- strat_folds(y, NFOLD)
  sc <- array(NA_real_, c(length(KGRID), length(LAM)))
  for (ki in seq_along(KGRID)) {
    pr <- matrix(NA_real_, length(y), length(LAM))
    for (w in 1:NFOLD) {
      it <- fi != w
      f <- fit_pipe(X[it, , drop = FALSE], y[it], KGRID[ki], LAM)
      pr[!it, ] <- score_pipe(f, X[!it, , drop = FALSE])
    }
    sc[ki, ] <- apply(pr, 2, function(s) auc(y, s))
  }
  bi <- which(sc == max(sc, na.rm = TRUE), arr.ind = TRUE)[1, ]
  kb <- KGRID[bi[1]]; lb <- LAM[bi[2]]
  aa <- c()
  for (rep in 1:NREP) {
    fo <- FOLDS[[rep]]
    for (v in 1:NFOLD) {
      tr <- fo != v
      f <- fit_pipe(X[tr, , drop = FALSE], y[tr], kb, LAM)
      aa <- c(aa, auc(y[!tr], as.numeric(score_pipe(f, X[!tr, , drop = FALSE], s = lb))))
    }
  }
  flat_rows[[nm]] <- data.frame(insieme = nm, auc_non_nidificata = mean(aa),
                                k_scelto = kb, lambda_scelto = lb)
}
flat <- do.call(rbind, flat_rows)

sm <- do.call(rbind, lapply(names(SETS), function(nm) {
  a <- fold$auc[fold$insieme == nm]
  data.frame(insieme = nm, n_par = length(SETS[[nm]]), auc_nidificata = mean(a),
             sd = sd(a), lo = quantile(a, .025, names = FALSE),
             hi = quantile(a, .975, names = FALSE))
}))
sm <- merge(sm, flat, by = "insieme", sort = FALSE)
sm$ottimismo <- sm$auc_non_nidificata - sm$auc_nidificata
sm <- sm[order(-sm$auc_nidificata), ]
write.csv(sm, "risultati/auc_nested.csv", row.names = FALSE)
write.csv(fold, "risultati/auc_fold.csv", row.names = FALSE)

# --- confronto appaiato per fold rispetto all'insieme completo ----------------
# I fold sono gli stessi per tutti gli insiemi dentro una ripetizione, quindi la
# differenza si puo' prendere fold per fold.
key <- function(nm) fold$auc[fold$insieme == nm]
cmp <- do.call(rbind, lapply(setdiff(names(SETS), "tutti e tre"), function(nm) {
  dd <- key("tutti e tre") - key(nm)
  data.frame(confronto = paste("tutti e tre -", nm), differenza = mean(dd),
             lo = quantile(dd, .025, names = FALSE),
             hi = quantile(dd, .975, names = FALSE),
             quota_fold_positivi = mean(dd > 0))
}))
write.csv(cmp, "risultati/auc_confronti.csv", row.names = FALSE)

cat("\n== AUC nidificata ==\n"); print(round(sm[, -1], 3)); print(sm$insieme)
cat("\n== confronti appaiati ==\n"); print(round(cmp[, -1], 3)); print(cmp$confronto)
