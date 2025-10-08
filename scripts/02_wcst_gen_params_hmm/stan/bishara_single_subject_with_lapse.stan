// bishara_single_subject_with_lapse.stan  (robusto, con rinormalizzazioni)

data {
  int<lower=1> T;
  array[T] int<lower=1, upper=4> resp_choice;
  array[T] int<lower=1, upper=4> resp_color;
  array[T] int<lower=1, upper=4> resp_shape;
  array[T] int<lower=1, upper=4> resp_number;
  array[T] int rew;  // 1 = corretto, altrimenti != 1
}

parameters {
  real<lower=0, upper=1> r;
  real<lower=0, upper=1> p;
  real<lower=0> d;
  real<lower=0, upper=1> lapse;  // mixing con uniforme su 4 scelte
}

model {
  // Priors (adatta se hai già definito altri)
  r ~ beta(2, 2);
  p ~ beta(2, 2);
  d ~ gamma(2, 1);
  lapse ~ beta(1, 20);

  {
    vector[3] a = rep_vector(1.0/3.0, 3);  // attenzione iniziale
    for (t in 1:T) {
      // Match matrix m (4x3)
      matrix[4,3] m = rep_matrix(0, 4, 3);
      m[resp_color[t],  1] = 1;
      m[resp_shape[t],  2] = 1;
      m[resp_number[t], 3] = 1;

      // Valori per pila = m * a^d
      vector[3] a_pow = pow(a, d);
      vector[4] pile_values = m * a_pow;

      // Probabilità base normalizzate in modo sicuro
      vector[4] p_base;
      real sum_pv = sum(pile_values);
      if (sum_pv > 0) p_base = pile_values / sum_pv;
      else            p_base = rep_vector(0.25, 4);  // fallback neutro

      // Lapse: mescola con uniforme
      vector[4] p_choice = (1 - lapse) * p_base + lapse * rep_vector(0.25, 4);

      target += categorical_lpmf(resp_choice[t] | p_choice);

      // Update attenzione (f = 1) con fallback se denom_s == 0
      int correct = (rew[t] == 1);
      vector[3] s;
      vector[3] mask;
      for (j in 1:3)
        mask[j] = correct ? m[resp_choice[t], j] : 1 - m[resp_choice[t], j];

      real denom_s = sum(a .* mask);
      if (denom_s > 0)
        s = (a .* mask) / denom_s;
      else
        s = a / sum(a);                 // fallback: nessun update effettivo

      if (correct == 1) a = (1 - r) * a + r * s;
      else              a = (1 - p) * a + p * s;

      // Rinormalizza sempre per stabilità numerica
      a = a / sum(a);
    }
  }
}

generated quantities {
  vector[T] log_lik_t;
  real log_lik;

  {
    vector[3] a = rep_vector(1.0/3.0, 3);
    log_lik = 0;
    for (t in 1:T) {
      matrix[4,3] m = rep_matrix(0, 4, 3);
      m[resp_color[t],  1] = 1;
      m[resp_shape[t],  2] = 1;
      m[resp_number[t], 3] = 1;

      vector[3] a_pow = pow(a, d);
      vector[4] pile_values = m * a_pow;

      vector[4] p_base;
      real sum_pv = sum(pile_values);
      if (sum_pv > 0) p_base = pile_values / sum_pv;
      else            p_base = rep_vector(0.25, 4);

      vector[4] p_choice = (1 - lapse) * p_base + lapse * rep_vector(0.25, 4);

      log_lik_t[t] = categorical_lpmf(resp_choice[t] | p_choice);
      log_lik += log_lik_t[t];

      int correct = (rew[t] == 1);
      vector[3] s;
      vector[3] mask;
      for (j in 1:3)
        mask[j] = correct ? m[resp_choice[t], j] : 1 - m[resp_choice[t], j];

      real denom_s = sum(a .* mask);
      if (denom_s > 0)
        s = (a .* mask) / denom_s;
      else
        s = a / sum(a);

      if (correct == 1) a = (1 - r) * a + r * s;
      else              a = (1 - p) * a + p * s;

      a = a / sum(a);
    }
  }
}
