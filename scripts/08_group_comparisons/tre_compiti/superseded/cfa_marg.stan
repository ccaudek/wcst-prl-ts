// Modello a fattore latente confermativo, verosimiglianza marginalizzata.
// Sigma_i = Lambda Phi Lambda' + diag(psi^2) + errore di misura noto (colonna WCST).
// La struttura delle medie e' identica in tutti i modelli confrontati
// (mu_p + b_p * gruppo), cosi' il confronto LOO isola la struttura di covarianza.
data {
  int<lower=1> N;
  int<lower=1> P;
  int<lower=1> F;
  matrix[N, P] y;                       // indicatori standardizzati
  vector[N] g;                          // 0 = controlli, 1 = pazienti
  vector[N] se_w;                       // sd di misura nota, scala standardizzata
  int<lower=0, upper=P> wcol;           // colonna con errore di misura noto (0 = nessuna)
  int<lower=0> Fa;                      // loading ancora, vincolati positivi
  array[Fa] int<lower=1, upper=P> anch;
  array[Fa] int<lower=1, upper=F> amap;
  int<lower=0> Nf;                      // loading liberi
  array[Nf] int<lower=1, upper=P> fidx;
  array[Nf] int<lower=1, upper=F> fmap;
  int<lower=0> Nfx;                     // loading fissati a un valore noto
  array[Nfx] int<lower=1, upper=P> fxidx;
  array[Nfx] int<lower=1, upper=F> fxmap;
  vector[Nfx] fxval;
  array[P] int<lower=0, upper=1> res_free;  // 0 = varianza residua fissata a zero
}
transformed data {
  int Npsi = 0;
  array[P] int psi_pos;
  for (p in 1:P) {
    if (res_free[p] == 1) { Npsi += 1; psi_pos[p] = Npsi; } else psi_pos[p] = 0;
  }
}
parameters {
  vector[P] mu;
  vector[P] b;
  vector<lower=0>[Fa] lam_a;
  vector[Nf] lam_f;
  vector<lower=0>[Npsi] psi_raw;
  cholesky_factor_corr[F] Lphi;
}
transformed parameters {
  matrix[P, F] Lam = rep_matrix(0, P, F);
  vector[P] psi = rep_vector(0, P);
  for (k in 1:Fa)  Lam[anch[k], amap[k]]   = lam_a[k];
  for (k in 1:Nf)  Lam[fidx[k], fmap[k]]   = lam_f[k];
  for (k in 1:Nfx) Lam[fxidx[k], fxmap[k]] = fxval[k];
  for (p in 1:P) if (res_free[p] == 1) psi[p] = psi_raw[psi_pos[p]];
}
model {
  matrix[P, F] LP = Lam * Lphi;
  matrix[P, P] S0 = LP * LP';
  mu ~ normal(0, 1);
  b ~ normal(0, 1);
  lam_a ~ normal(0, 0.5);               // priore regolarizzante
  lam_f ~ normal(0, 0.5);
  psi_raw ~ normal(0, 1);
  Lphi ~ lkj_corr_cholesky(2);
  for (p in 1:P) S0[p, p] += square(psi[p]) + 1e-6;
  for (i in 1:N) {
    matrix[P, P] S = S0;
    if (wcol > 0) S[wcol, wcol] += square(se_w[i]);
    y[i] ~ multi_normal(mu + b * g[i], S);
  }
}
generated quantities {
  vector[N] log_lik;
  matrix[P, P] Sigma;
  matrix[F, F] Phi = multiply_lower_tri_self_transpose(Lphi);
  {
    matrix[P, F] LP = Lam * Lphi;
    Sigma = LP * LP';
    for (p in 1:P) Sigma[p, p] += square(psi[p]) + 1e-6;
    for (i in 1:N) {
      matrix[P, P] S = Sigma;
      if (wcol > 0) S[wcol, wcol] += square(se_w[i]);
      log_lik[i] = multi_normal_lpdf(y[i] | mu + b * g[i], S);
    }
  }
}
