// rw_feature_based.stan
// Rescorla–Wagner feature-based per WCST
// - 3 valori per dimensione (color, shape, number)
// - Valore della pila = somma dei valori delle dimensioni che matchano la pila
// - Scelta: softmax su 4 pile (beta), opzionale sticky (kappa) e lapse
// - Update: V_j <- V_j + alpha_{sign(pe)} * pe * m[choice, j]
// - Emette log_lik_t[T] per PSIS-LOO
//
// Dati attesi:
//   choice in 1..4
//   rew, los binari (0/1) come nelle varianti RW
//   resp_* sono INDICI 1..4: la pila che matcha la dimensione al trial t

data {
  int<lower=1> T;
  array[T] int<lower=1, upper=4> choice;
  array[T] int<lower=0, upper=1> rew;       // 1 reward, 0 altrimenti
  array[T] int<lower=0, upper=1> los;       // 1 loss, 0 altrimenti
  array[T] int<lower=1, upper=4> resp_color;
  array[T] int<lower=1, upper=4> resp_shape;
  array[T] int<lower=1, upper=4> resp_number;
  int<lower=0, upper=1> use_kappa;
  int<lower=0, upper=1> use_lapse;
}

parameters {
  real alpha_pos_raw;
  real alpha_neg_raw;
  real beta_raw;
  real kappa;         // sticky bias (usato solo se use_kappa==1)
  real lapse_raw;     // in (0,1) via inv_logit
}

transformed parameters {
  real<lower=0, upper=1> alpha_pos = inv_logit(alpha_pos_raw);
  real<lower=0, upper=1> alpha_neg = inv_logit(alpha_neg_raw);
  real<lower=0>          beta      = exp(beta_raw);
  real<lower=0, upper=1> lapse     = inv_logit(lapse_raw);
}

model {
  // Priors deboli (adatta se preferisci)
  alpha_pos_raw ~ normal(0, 1);
  alpha_neg_raw ~ normal(0, 1);
  beta_raw      ~ normal(0, 0.7);
  kappa         ~ student_t(3, 0, 1);
  lapse_raw     ~ normal(-2.2, 1.0);

  {
    vector[3] V = rep_vector(0.0, 3);  // valori per dimensione (init 0)
    int prev = 0;

    for (t in 1:T) {
      // Matrice di match m (4x3): 1 se la pila k matcha la dimensione j al trial t
      matrix[4,3] m = rep_matrix(0, 4, 3);
      m[resp_color[t],  1] = 1;
      m[resp_shape[t],  2] = 1;
      m[resp_number[t], 3] = 1;

      // Valori per le 4 pile: prodotto matrice–vettore
      vector[4] pile_values = m * V;

      // Scelta: softmax (logit) con sticky opzionale
      vector[4] logits = beta * pile_values;
      if (use_kappa == 1 && prev > 0) logits[prev] += kappa;

      real lp_model = categorical_logit_lpmf(choice[t] | logits);
      if (use_lapse == 1)
        target += log_mix(1.0 - lapse, lp_model, -log(4.0));
      else
        target += lp_model;

      // Update TD feature-based solo sulle dimensioni della pila scelta
      int a = choice[t];
      real r = (rew[t] == 1 ? 1.0 : 0.0) - (los[t] == 1 ? 1.0 : 0.0); // +1/-1
      real vhat = m[a] * V;   // row_vector[3] * vector[3] -> real
      real pe = r - vhat;

      if (pe >= 0)
        V += alpha_pos * pe * m[a]';  // m[a]' è vector[3]
      else
        V += alpha_neg * pe * m[a]';

      prev = a;
    }
  }
}

generated quantities {
  array[T] real log_lik_t;

  {
    vector[3] V = rep_vector(0.0, 3);
    int prev = 0;

    for (t in 1:T) {
      matrix[4,3] m = rep_matrix(0, 4, 3);
      m[resp_color[t],  1] = 1;
      m[resp_shape[t],  2] = 1;
      m[resp_number[t], 3] = 1;

      vector[4] pile_values = m * V;

      vector[4] logits = beta * pile_values;
      if (use_kappa == 1 && prev > 0) logits[prev] += kappa;

      real lp_model = categorical_logit_lpmf(choice[t] | logits);
      if (use_lapse == 1)
        log_lik_t[t] = log_mix(1.0 - lapse, lp_model, -log(4.0));
      else
        log_lik_t[t] = lp_model;

      int a = choice[t];
      real r = (rew[t] == 1 ? 1.0 : 0.0) - (los[t] == 1 ? 1.0 : 0.0);
      real vhat = m[a] * V;
      real pe = r - vhat;

      if (pe >= 0)
        V += alpha_pos * pe * m[a]';
      else
        V += alpha_neg * pe * m[a]';

      prev = a;
    }
  }
}
