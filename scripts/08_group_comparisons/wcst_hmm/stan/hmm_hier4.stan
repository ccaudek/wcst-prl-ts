// hmm_hier4.stan -----------------------------------------------------------
// Come hmm_hier.stan, ma TUTTI e quattro i parametri (h, d, lapse, eta)
// variano per soggetto. Serve a verificare se i dati sostengono differenze
// individuali anche in lapse ed eta, oppure se le loro SD a posteriori
// restano determinate dal prior.

data {
  int<lower=1> N;
  int<lower=1> T;
  array[N, T] int<lower=1, upper=4> choice;
  array[N, T] int<lower=1, upper=4> kc;
  array[N, T] int<lower=1, upper=4> ks;
  array[N, T] int<lower=1, upper=4> kn;
  array[N, T] int<lower=0, upper=1> rew;
  array[N, T] int<lower=0, upper=1> valid;
  vector<lower=0, upper=1>[N] grp;
}

transformed data {
  array[N, T] matrix[4, 3] m;
  for (i in 1:N) {
    for (t in 1:T) {
      m[i, t] = rep_matrix(0.0, 4, 3);
      m[i, t][kc[i, t], 1] = 1.0;
      m[i, t][ks[i, t], 2] = 1.0;
      m[i, t][kn[i, t], 3] = 1.0;
    }
  }
}

parameters {
  vector[4] mu;                    // logit h, log d, logit lapse, logit eta
  vector[4] bgrp;
  vector<lower=0>[4] sigma;
  cholesky_factor_corr[4] Lcorr;
  matrix[4, N] z;
}

transformed parameters {
  matrix[N, 4] theta;
  vector<lower=0, upper=1>[N] h;
  vector<lower=0>[N] d;
  vector<lower=0, upper=1>[N] lapse;
  vector<lower=0, upper=1>[N] eta;
  {
    matrix[4, N] r = diag_pre_multiply(sigma, Lcorr) * z;
    for (k in 1:4) {
      theta[, k] = mu[k] + bgrp[k] * grp + to_vector(r[k]);
    }
  }
  h     = inv_logit(theta[, 1]);
  d     = exp(theta[, 2]);
  lapse = inv_logit(theta[, 3]);
  // eta < 0.5 per identificabilita' (vedi hmm_hier.stan)
  eta   = 0.5 * inv_logit(theta[, 4]);
}

model {
  mu[1] ~ normal(-1.5, 1.0);
  mu[2] ~ normal(1.8, 0.7);
  mu[3] ~ normal(-4.0, 1.0);
  mu[4] ~ normal(-2.8, 1.0);
  bgrp  ~ normal(0, 0.5);
  sigma ~ normal(0, 1);
  Lcorr ~ lkj_corr_cholesky(2);
  to_vector(z) ~ std_normal();

  for (i in 1:N) {
    vector[3] b = rep_vector(1.0 / 3, 3);
    for (t in 1:T) {
      vector[4] p = softmax(d[i] * (m[i, t] * b));
      p = (1 - lapse[i]) * p + lapse[i] / 4;
      if (valid[i, t] == 1) {
        target += log(p[choice[i, t]]);
      }
      vector[3] matched = m[i, t][choice[i, t]]';
      vector[3] pr1 = matched * (1 - eta[i]) + (1 - matched) * eta[i];
      vector[3] lik = rew[i, t] == 1 ? pr1 : 1 - pr1;
      b = b .* lik;
      b /= sum(b);
      b = (1 - h[i]) * b + h[i] * (1 - b) / 2;
    }
  }
}

generated quantities {
  array[N, T] real log_lik;
  matrix[4, 4] Rho = multiply_lower_tri_self_transpose(Lcorr);
  for (i in 1:N) {
    vector[3] b = rep_vector(1.0 / 3, 3);
    for (t in 1:T) {
      vector[4] p = softmax(d[i] * (m[i, t] * b));
      p = (1 - lapse[i]) * p + lapse[i] / 4;
      log_lik[i, t] = valid[i, t] == 1 ? log(p[choice[i, t]]) : 0.0;
      vector[3] matched = m[i, t][choice[i, t]]';
      vector[3] pr1 = matched * (1 - eta[i]) + (1 - matched) * eta[i];
      vector[3] lik = rew[i, t] == 1 ? pr1 : 1 - pr1;
      b = b .* lik;
      b /= sum(b);
      b = (1 - h[i]) * b + h[i] * (1 - b) / 2;
    }
  }
}
