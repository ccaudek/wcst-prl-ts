// hmm_sticky.stan ----------------------------------------------------------
// hmm_hier.stan + perseverazione a livello di dimensione.
//
// Motivazione empirica: il PPC di hmm_hier sottostima il win-stay (0.978
// osservato vs 0.956 predetto) e la ripetizione della dimensione (0.866 vs
// 0.842), e sovrastima il lose-shift (0.671 vs 0.731) e la velocita' di
// recupero dopo il cambio di regola. Tutti segni della stessa cosa: i
// partecipanti ripetono la dimensione appena seguita piu' di quanto la sola
// credenza sulla regola giustifichi.
//
// kappa aggiunge un bonus alle pile che corrispondono alla dimensione seguita
// al trial precedente, indipendentemente dalla credenza. Varia per soggetto:
// nel WCST la perseverazione e' la misura clinica di riferimento, quindi la
// sua variabilita' individuale e' di interesse sostanziale.

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
  vector[3] mu;                    // logit h, log d, kappa
  vector[3] bgrp;
  vector<lower=0>[3] sigma;
  cholesky_factor_corr[3] Lcorr;
  matrix[3, N] z;
  real mu_lapse;
  real mu_eta;
}

transformed parameters {
  matrix[N, 3] theta;
  vector<lower=0, upper=1>[N] h;
  vector<lower=0>[N] d;
  vector[N] kappa;
  real<lower=0, upper=1> lapse = inv_logit(mu_lapse);
  real<lower=0, upper=0.5> eta = 0.5 * inv_logit(mu_eta);
  {
    matrix[3, N] r = diag_pre_multiply(sigma, Lcorr) * z;
    for (k in 1:3) {
      theta[, k] = mu[k] + bgrp[k] * grp + to_vector(r[k]);
    }
  }
  h     = inv_logit(theta[, 1]);
  d     = exp(theta[, 2]);
  kappa = theta[, 3];
}

model {
  mu[1] ~ normal(-1.5, 1.0);
  mu[2] ~ normal(1.8, 0.7);
  mu[3] ~ normal(0, 1.0);
  bgrp  ~ normal(0, 0.5);
  sigma ~ normal(0, 1);
  Lcorr ~ lkj_corr_cholesky(2);
  to_vector(z) ~ std_normal();
  mu_lapse ~ normal(-4.0, 1.0);
  mu_eta   ~ normal(-2.8, 1.0);

  for (i in 1:N) {
    vector[3] b = rep_vector(1.0 / 3, 3);
    vector[3] prev = rep_vector(0.0, 3);
    for (t in 1:T) {
      vector[4] p = softmax(m[i, t] * (d[i] * b + kappa[i] * prev));
      p = (1 - lapse) * p + lapse / 4;
      if (valid[i, t] == 1) {
        target += log(p[choice[i, t]]);
      }
      vector[3] matched = m[i, t][choice[i, t]]';
      vector[3] pr1 = matched * (1 - eta) + (1 - matched) * eta;
      vector[3] lik = rew[i, t] == 1 ? pr1 : 1 - pr1;
      b = b .* lik;
      b /= sum(b);
      b = (1 - h[i]) * b + h[i] * (1 - b) / 2;
      prev = matched;
    }
  }
}

generated quantities {
  array[N, T] real log_lik;
  matrix[3, 3] Rho = multiply_lower_tri_self_transpose(Lcorr);
  for (i in 1:N) {
    vector[3] b = rep_vector(1.0 / 3, 3);
    vector[3] prev = rep_vector(0.0, 3);
    for (t in 1:T) {
      vector[4] p = softmax(m[i, t] * (d[i] * b + kappa[i] * prev));
      p = (1 - lapse) * p + lapse / 4;
      log_lik[i, t] = valid[i, t] == 1 ? log(p[choice[i, t]]) : 0.0;
      vector[3] matched = m[i, t][choice[i, t]]';
      vector[3] pr1 = matched * (1 - eta) + (1 - matched) * eta;
      vector[3] lik = rew[i, t] == 1 ? pr1 : 1 - pr1;
      b = b .* lik;
      b /= sum(b);
      b = (1 - h[i]) * b + h[i] * (1 - b) / 2;
      prev = matched;
    }
  }
}
