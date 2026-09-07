// Modello a componenti di varianza per un disegno multitratto-multimetodo.
//   y_ip = mu_p + b_p * gruppo_i + sum_k a_kp * F_ik[grp(p,k)] + e_ip
// Ogni raggruppamento k definisce un insieme di fattori ortogonali: le primitive
// computazionali (tratto, condiviso TRA compiti), i compiti (metodo), e le coppie
// primitiva-entro-compito (necessarie perche' la varianza di tratto non assorba
// la correlazione tra due condizioni dello stesso compito).
// Loading non negativi: il modello puo' solo trovare varianza condivisa, mai
// inventarne di negativa, quindi il vincolo e' conservativo rispetto alla
// conclusione di assenza di struttura di tratto.
// Verosimiglianza marginalizzata: Sigma = sum_k Lk Lk' + diag(psi^2) + errore noto.
// Le medie sono identiche in tutte le strutture confrontate, quindi LOO isola
// la struttura di covarianza.
data {
  int<lower=1> N;
  int<lower=1> P;
  int<lower=1> K;                           // numero di raggruppamenti
  matrix[N, P] y;                           // indicatori standardizzati
  vector[N] g;
  matrix<lower=0>[N, P] se;                 // sd di misura nota per soggetto e
                                            // indicatore (0 = non disponibile)
  array[P, K] int<lower=0> grp;             // 0 = nessun carico su quel raggruppamento
  array[K] int<lower=0, upper=1> use;
}
transformed data {
  int na = 0;
  array[P, K] int ap;
  for (k in 1:K) for (p in 1:P) {
    if (grp[p, k] > 0 && use[k] == 1) { na += 1; ap[p, k] = na; } else ap[p, k] = 0;
  }
}
parameters {
  vector[P] mu;
  vector[P] b;
  vector<lower=0>[na] a_raw;
  vector<lower=0>[P] psi;
}
transformed parameters {
  matrix[P, K] a = rep_matrix(0, P, K);
  for (k in 1:K) for (p in 1:P) if (ap[p, k] > 0) a[p, k] = a_raw[ap[p, k]];
}
model {
  matrix[P, P] S0 = diag_matrix(square(psi) + 1e-6);
  mu ~ normal(0, 1);
  b ~ normal(0, 1);
  a_raw ~ normal(0, 0.5);
  psi ~ normal(0, 1);
  for (k in 1:K) for (p in 1:P) for (q in 1:P)
    if (grp[p, k] > 0 && grp[p, k] == grp[q, k]) S0[p, q] += a[p, k] * a[q, k];
  for (i in 1:N) {
    matrix[P, P] S = S0;
    for (p in 1:P) S[p, p] += square(se[i, p]);
    y[i] ~ multi_normal(mu + b * g[i], S);
  }
}
generated quantities {
  vector[N] log_lik;
  matrix[P, K] var_comp;
  vector[P] var_unica;
  {
    matrix[P, P] S0 = diag_matrix(square(psi) + 1e-6);
    for (k in 1:K) for (p in 1:P) for (q in 1:P)
      if (grp[p, k] > 0 && grp[p, k] == grp[q, k]) S0[p, q] += a[p, k] * a[q, k];
    for (i in 1:N) {
      matrix[P, P] S = S0;
      for (p in 1:P) S[p, p] += square(se[i, p]);
      log_lik[i] = multi_normal_lpdf(y[i] | mu + b * g[i], S);
    }
    for (p in 1:P) {
      real tot = square(psi[p]);
      for (k in 1:K) tot += square(a[p, k]);
      for (k in 1:K) var_comp[p, k] = square(a[p, k]) / tot;
      var_unica[p] = square(psi[p]) / tot;
    }
  }
}
