// hmm_hier.stan ------------------------------------------------------------
// Modello di inferenza della regola per il WCST, versione GERARCHICA.
//
// Tre stati latenti (colore, forma, numero). A ogni trial la regola puo'
// cambiare con probabilita' h (hazard). Il partecipante mantiene una credenza
// b sulla regola vigente, la aggiorna col feedback (rumore eta) e scegle la
// pila con una softmax di consistenza d, piu' una quota di risposte casuali
// (lapse).
//
// Parametri variabili per soggetto: h (flessibilita') e d (consistenza).
//   lapse ed eta sono solo a livello di popolazione: sui 60 trial disponibili
//   i loro profili di verosimiglianza individuali sono monotoni verso il
//   bordo, quindi una loro stima per soggetto sarebbe determinata dal prior.
// Effetto di gruppo stimato dentro il modello (b_h, b_d), su grp = 1 per i
// pazienti.

data {
  int<lower=1> N;
  int<lower=1> T;
  array[N, T] int<lower=1, upper=4> choice;
  array[N, T] int<lower=1, upper=4> kc;   // pila che corrisponde al colore della carta
  array[N, T] int<lower=1, upper=4> ks;   // ... alla forma
  array[N, T] int<lower=1, upper=4> kn;   // ... al numero
  array[N, T] int<lower=0, upper=1> rew;
  array[N, T] int<lower=0, upper=1> valid; // 0 = trial senza risposta
  vector<lower=0, upper=1>[N] grp;
}

transformed data {
  // m[i,t] : matrice 4x3, m[k,j] = 1 se la pila k corrisponde alla carta sulla dimensione j
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
  real mu_h;                       // media di popolazione, scala logit
  real mu_d;                       // media di popolazione, scala log
  real b_h;                        // effetto di gruppo su h
  real b_d;                        // effetto di gruppo su d
  vector<lower=0>[2] sigma;        // SD tra soggetti di (logit h, log d)
  cholesky_factor_corr[2] Lcorr;
  matrix[2, N] z;                  // parametrizzazione non centrata
  real mu_lapse;
  real mu_eta;
}

transformed parameters {
  vector[N] logit_h;
  vector[N] log_d;
  vector<lower=0, upper=1>[N] h;
  vector<lower=0>[N] d;
  real<lower=0, upper=1> lapse = inv_logit(mu_lapse);
  // eta < 0.5 per identificabilita': con eta > 0.5 il feedback verrebbe letto
  // al contrario e il modello ha un secondo modo equivalente (simmetria).
  real<lower=0, upper=0.5> eta = 0.5 * inv_logit(mu_eta);
  {
    matrix[2, N] r = diag_pre_multiply(sigma, Lcorr) * z;
    logit_h = mu_h + b_h * grp + to_vector(r[1]);
    log_d   = mu_d + b_d * grp + to_vector(r[2]);
  }
  h = inv_logit(logit_h);
  d = exp(log_d);
}

model {
  mu_h  ~ normal(-1.5, 1.0);        // h ~ 0.18
  mu_d  ~ normal(1.8, 0.7);         // d ~ 6
  b_h   ~ normal(0, 0.5);
  b_d   ~ normal(0, 0.5);
  sigma ~ normal(0, 1);             // half-normal
  Lcorr ~ lkj_corr_cholesky(2);
  to_vector(z) ~ std_normal();
  mu_lapse ~ normal(-4.0, 1.0);     // lapse ~ 0.018
  mu_eta   ~ normal(-2.8, 1.0);     // eta   ~ 0.5*inv_logit(-2.8) = 0.029

  for (i in 1:N) {
    vector[3] b = rep_vector(1.0 / 3, 3);
    for (t in 1:T) {
      vector[4] p = softmax(d[i] * (m[i, t] * b));
      p = (1 - lapse) * p + lapse / 4;
      if (valid[i, t] == 1) {
        target += log(p[choice[i, t]]);
      }
      // aggiornamento della credenza dato il feedback
      vector[3] matched = m[i, t][choice[i, t]]';
      vector[3] pr1 = matched * (1 - eta) + (1 - matched) * eta;
      vector[3] lik = rew[i, t] == 1 ? pr1 : 1 - pr1;
      b = b .* lik;
      b /= sum(b);
      b = (1 - h[i]) * b + h[i] * (1 - b) / 2;
    }
  }
}

generated quantities {
  array[N, T] real log_lik;
  array[N, T] int choice_pred;
  array[N, T] real p_correct;       // prob. attesa di scelta corretta (per PPC)
  vector[2] tau;                    // SD tra soggetti, riportate a scala interpretabile
  real rho_hd = multiply_lower_tri_self_transpose(Lcorr)[1, 2];

  tau[1] = sigma[1];
  tau[2] = sigma[2];
  for (i in 1:N) {
    vector[3] b = rep_vector(1.0 / 3, 3);
    for (t in 1:T) {
      vector[4] p = softmax(d[i] * (m[i, t] * b));
      p = (1 - lapse) * p + lapse / 4;
      log_lik[i, t] = valid[i, t] == 1 ? log(p[choice[i, t]]) : 0.0;
      choice_pred[i, t] = categorical_rng(p);
      // pila corretta = quella indicata dalla dimensione vigente; ricavata dal feedback
      p_correct[i, t] = p[choice[i, t]];
      vector[3] matched = m[i, t][choice[i, t]]';
      vector[3] pr1 = matched * (1 - eta) + (1 - matched) * eta;
      vector[3] lik = rew[i, t] == 1 ? pr1 : 1 - pr1;
      b = b .* lik;
      b /= sum(b);
      b = (1 - h[i]) * b + h[i] * (1 - b) / 2;
    }
  }
}
