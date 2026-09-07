# funs_simulate.R ----------------------------------------------------------
# Simulazione generativa dell'HMM di inferenza della regola sul WCST reale.
# Il feedback e' quello del compito (deterministico: la scelta e' corretta se
# la pila coincide con quella indicata dalla regola vigente); eta entra solo
# nel modo in cui il partecipante *interpreta* il feedback.

sim_hmm <- function(sdh, h, d, lapse, eta, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  N <- sdh$N; T <- sdh$T
  ch <- matrix(0L, N, T); rw <- matrix(0L, N, T)
  K <- array(0L, c(N, T, 3))
  K[, , 1] <- sdh$k_color; K[, , 2] <- sdh$k_shape; K[, , 3] <- sdh$k_number
  for (i in 1:N) {
    b <- rep(1 / 3, 3)
    for (t in 1:T) {
      m <- matrix(0, 4, 3)
      m[K[i, t, 1], 1] <- 1; m[K[i, t, 2], 2] <- 1; m[K[i, t, 3], 3] <- 1
      z <- d[i] * as.vector(m %*% b); z <- z - max(z)
      p <- exp(z); p <- p / sum(p)
      p <- (1 - lapse[i]) * p + lapse[i] / 4
      a <- sample.int(4, 1, prob = p)
      correct_pile <- K[i, t, sdh$rule[i, t]]
      r <- as.integer(a == correct_pile)
      ch[i, t] <- a; rw[i, t] <- r
      matched <- m[a, ]
      pr1 <- matched * (1 - eta[i]) + (1 - matched) * eta[i]
      lik <- if (r == 1) pr1 else 1 - pr1
      b <- b * lik; b <- b / sum(b)
      b <- (1 - h[i]) * b + h[i] * (1 - b) / 2
    }
  }
  list(choice = ch, rew = rw)
}

sim_hmm_sticky <- function(sdh, h, d, kappa, lapse, eta, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  N <- sdh$N; T <- sdh$T
  ch <- matrix(0L, N, T); rw <- matrix(0L, N, T)
  K <- array(0L, c(N, T, 3))
  K[, , 1] <- sdh$k_color; K[, , 2] <- sdh$k_shape; K[, , 3] <- sdh$k_number
  for (i in 1:N) {
    b <- rep(1 / 3, 3); prev <- rep(0, 3)
    for (t in 1:T) {
      m <- matrix(0, 4, 3)
      m[K[i, t, 1], 1] <- 1; m[K[i, t, 2], 2] <- 1; m[K[i, t, 3], 3] <- 1
      z <- as.vector(m %*% (d[i] * b + kappa[i] * prev)); z <- z - max(z)
      p <- exp(z); p <- p / sum(p)
      p <- (1 - lapse[i]) * p + lapse[i] / 4
      a <- sample.int(4, 1, prob = p)
      r <- as.integer(a == K[i, t, sdh$rule[i, t]])
      ch[i, t] <- a; rw[i, t] <- r
      matched <- m[a, ]
      pr1 <- matched * (1 - eta[i]) + (1 - matched) * eta[i]
      lik <- if (r == 1) pr1 else 1 - pr1
      b <- b * lik; b <- b / sum(b)
      b <- (1 - h[i]) * b + h[i] * (1 - b) / 2
      prev <- matched
    }
  }
  list(choice = ch, rew = rw)
}

# Indici comportamentali classici del WCST, calcolati allo stesso modo su dati
# osservati e simulati.
wcst_signatures <- function(sdh, choice, rew) {
  N <- sdh$N; T <- sdh$T
  K <- array(0L, c(N, T, 3))
  K[, , 1] <- sdh$k_color; K[, , 2] <- sdh$k_shape; K[, , 3] <- sdh$k_number
  out <- data.frame(subj = 1:N)
  out$acc <- rowMeans(rew)
  # dimensione seguita a ogni trial: quella (o quelle) che la scelta soddisfa
  matched <- array(0L, c(N, T, 3))
  for (j in 1:3) matched[, , j] <- (choice == K[, , j])
  # errore perseverativo: errore in cui la scelta segue la regola del blocco precedente
  prev_rule <- cbind(NA_integer_, t(apply(sdh$rule, 1, function(r) {
    pr <- rep(NA_integer_, length(r)); cur <- r[1]; last <- NA_integer_
    for (t in seq_along(r)) { if (r[t] != cur) { last <- cur; cur <- r[t] }; pr[t] <- last }
    pr[-1]
  })))
  pers <- non_pers <- matrix(0L, N, T)
  for (i in 1:N) for (t in 1:T) {
    if (rew[i, t] == 0) {
      pr <- prev_rule[i, t]
      if (!is.na(pr) && matched[i, t, pr] == 1) pers[i, t] <- 1L else non_pers[i, t] <- 1L
    }
  }
  out$prop_pers_err <- rowSums(pers) / T
  out$prop_non_pers_err <- rowSums(non_pers) / T
  # perseverazione della dimensione: quanto spesso si ripete la dimensione seguita
  rep_dim <- matrix(NA, N, T - 1)
  for (i in 1:N) for (t in 2:T) rep_dim[i, t - 1] <- as.numeric(sum(matched[i, t, ] * matched[i, t - 1, ]) > 0)
  out$prop_rep_dim <- rowMeans(rep_dim)
  # win-stay / lose-shift a livello di dimensione
  ws <- ls_ <- matrix(NA, N, T - 1)
  for (i in 1:N) for (t in 2:T) {
    same <- sum(matched[i, t, ] * matched[i, t - 1, ]) > 0
    if (rew[i, t - 1] == 1) ws[i, t - 1] <- as.numeric(same) else ls_[i, t - 1] <- as.numeric(!same)
  }
  out$win_stay  <- rowMeans(ws, na.rm = TRUE)
  out$lose_shift <- rowMeans(ls_, na.rm = TRUE)
  # accuratezza per posizione nel blocco (costo dello switch)
  pos <- sdh$pos_in_block
  out$acc_pos1_3 <- sapply(1:N, function(i) mean(rew[i, pos[i, ] <= 3]))
  out$acc_pos8_10 <- sapply(1:N, function(i) mean(rew[i, pos[i, ] >= 8]))
  out
}
