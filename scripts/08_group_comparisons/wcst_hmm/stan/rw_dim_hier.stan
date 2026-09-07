// rw_dim_hier.stan ---------------------------------------------------------
// Rescorla-Wagner gerarchico sulle DIMENSIONI della carta (colore, forma,
// numero) con perseverazione a livello di dimensione.
//
// E' il modello RL propriamente detto, usato come alternativa all'HMM di
// inferenza della regola. Differisce dal modello RW della pipeline precedente
// per un punto essenziale: i valori sono attaccati alle tre dimensioni, non
// alle quattro pile. L'identita' della pila non porta informazione sulla
// regola, quindi un RW sulle pile non puo' superare il caso.
//
// Parametri per soggetto: alpha_pos, alpha_neg (tassi di apprendimento dopo
// feedback positivo/negativo), beta (inverso della temperatura), kappa
// (tendenza a riscegliere la dimensione seguita al trial precedente).

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
  vector[4] mu;                    // scala non vincolata: logit ap, logit an, log beta, kappa
  vector[4] bgrp;
  vector<lower=0>[4] sigma;
  cholesky_factor_corr[4] Lcorr;
  matrix[4, N] z;
}

transformed parameters {
  matrix[N, 4] theta;
  vector<lower=0, upper=1>[N] alpha_pos;
  vector<lower=0, upper=1>[N] alpha_neg;
  vector<lower=0>[N] beta;
  vector[N] kappa;
  {
    matrix[4, N] r = diag_pre_multiply(sigma, Lcorr) * z;
    for (k in 1:4) {
      theta[, k] = mu[k] + bgrp[k] * grp + to_vector(r[k]);
    }
  }
  alpha_pos = inv_logit(theta[, 1]);
  alpha_neg = inv_logit(theta[, 2]);
  beta      = exp(theta[, 3]);
  kappa     = theta[, 4];
}

model {
  mu[1] ~ normal(0, 1.2);
  mu[2] ~ normal(0, 1.2);
  mu[3] ~ normal(1.0, 0.7);
  mu[4] ~ normal(0, 1.0);
  bgrp  ~ normal(0, 0.5);
  sigma ~ normal(0, 1);
  Lcorr ~ lkj_corr_cholesky(2);
  to_vector(z) ~ std_normal();

  for (i in 1:N) {
    vector[3] V = rep_vector(0.0, 3);
    vector[3] prev = rep_vector(0.0, 3);
    for (t in 1:T) {
      vector[4] val = m[i, t] * (V + kappa[i] * prev);
      vector[4] p = softmax(beta[i] * val);
      if (valid[i, t] == 1) {
        target += log(p[choice[i, t]]);
      }
      vector[3] matched = m[i, t][choice[i, t]]';
      real r_signed = rew[i, t] == 1 ? 1.0 : -1.0;
      real a = rew[i, t] == 1 ? alpha_pos[i] : alpha_neg[i];
      V += a * (r_signed - V) .* matched;
      prev = matched;
    }
  }
}

generated quantities {
  array[N, T] real log_lik;
  array[N, T] int choice_pred;
  matrix[4, 4] Rho = multiply_lower_tri_self_transpose(Lcorr);

  for (i in 1:N) {
    vector[3] V = rep_vector(0.0, 3);
    vector[3] prev = rep_vector(0.0, 3);
    for (t in 1:T) {
      vector[4] val = m[i, t] * (V + kappa[i] * prev);
      vector[4] p = softmax(beta[i] * val);
      log_lik[i, t] = valid[i, t] == 1 ? log(p[choice[i, t]]) : 0.0;
      choice_pred[i, t] = categorical_rng(p);
      vector[3] matched = m[i, t][choice[i, t]]';
      real r_signed = rew[i, t] == 1 ? 1.0 : -1.0;
      real a = rew[i, t] == 1 ? alpha_pos[i] : alpha_neg[i];
      V += a * (r_signed - V) .* matched;
      prev = matched;
    }
  }
}
