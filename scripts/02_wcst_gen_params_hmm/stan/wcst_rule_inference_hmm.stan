
// wcst_rule_inference_hmm.stan
// HMM "rule-inference" per WCST con 3 regole (color, shape, number).
// - Stato latente z_t ∈ {1,2,3} = regola corrente
// - Belief b_t aggiornato da feedback (correct/incorrect) e hazard di switch h
// - Scelta: mixture-of-softmax sulle 4 pile pesata da b_t
// - Emissione feedback: parametrizzata con eta (rumore feedback)
// - Lapse: miscela con uniforme su 4 scelte
// - Emette log_lik_t[T] per PSIS-LOO
//
// Data: resp_* INDICI 1..4; correct = (rew==1) binario

data {
  int<lower=1> T;
  array[T] int<lower=1, upper=4> resp_choice;   // scelta 1..4
  array[T] int<lower=1, upper=4> resp_color;    // pile corrette per ogni regola
  array[T] int<lower=1, upper=4> resp_shape;
  array[T] int<lower=1, upper=4> resp_number;
  array[T] int<lower=0, upper=1> correct;       // 1 se feedback "corretto", 0 se "errato"
}

parameters {
  real<lower=0, upper=1> h;       // hazard di cambio regola
  real<lower=0>          d;       // consistenza decisionale (softmax)
  real<lower=0, upper=1> lapse;   // lapse su 4 scelte
  real<lower=0, upper=1> eta;     // rumore nel feedback (P(feedback errato))
}

model {
  // Priors
  h     ~ beta(1, 9);    // prior per switch rari
  d     ~ gamma(2, 1);   // moderata consistenza
  lapse ~ beta(1, 20);   // piccolo
  eta   ~ beta(2, 20);   // feedback quasi affidabile

  {
    vector[3] b = rep_vector(1.0/3.0, 3);  // belief iniziale
    for (t in 1:T) {
      // Indice della pila "corretta" sotto ciascuna regola
      array[3] int k_of_rule;
      k_of_rule[1] = resp_color[t];
      k_of_rule[2] = resp_shape[t];
      k_of_rule[3] = resp_number[t];

      // Politica: mixture of softmax sulle 4 pile
      vector[4] numer = rep_vector(0.0, 4);
      for (z in 1:3) {
        vector[4] u = rep_vector(0.0, 4);
        u[k_of_rule[z]] = 1.0;                     // utilità 1 sul mazzo “corretto” per quella regola
        numer += b[z] * exp(d * u);               // softmax locale, poi media pesata
      }
      vector[4] p_base = numer / sum(numer);
      vector[4] p_choice = (1 - lapse) * p_base + lapse * rep_vector(0.25, 4);

      target += categorical_lpmf(resp_choice[t] | p_choice);

      // Likelihood del feedback per ogni regola
      vector[3] L;
      for (z in 1:3) {
        int is_correct_under_z = (resp_choice[t] == k_of_rule[z]);
        if (correct[t] == 1)
          L[z] = is_correct_under_z ? (1 - eta) : eta;
        else
          L[z] = is_correct_under_z ? eta : (1 - eta);
      }

      // Bayes update (se la somma è 0 per numerica, fallback a uniforme)
      vector[3] post = b .* L;
      real Z = sum(post);
      if (Z > 0) post /= Z;
      else       post = rep_vector(1.0/3.0, 3);

      // Transizione con hazard h (primo ordine, passaggio al prossimo trial)
      b = (1 - h) * post + (h / 2.0) * (rep_vector(1.0, 3) - post);
    }
  }
}

generated quantities {
  vector[T] log_lik_t;
  real log_lik;

  {
    vector[3] b = rep_vector(1.0/3.0, 3);
    log_lik = 0;
    for (t in 1:T) {
      array[3] int k_of_rule;
      k_of_rule[1] = resp_color[t];
      k_of_rule[2] = resp_shape[t];
      k_of_rule[3] = resp_number[t];

      vector[4] numer = rep_vector(0.0, 4);
      for (z in 1:3) {
        vector[4] u = rep_vector(0.0, 4);
        u[k_of_rule[z]] = 1.0;
        numer += b[z] * exp(d * u);
      }
      vector[4] p_base = numer / sum(numer);
      vector[4] p_choice = (1 - lapse) * p_base + lapse * rep_vector(0.25, 4);

      log_lik_t[t] = categorical_lpmf(resp_choice[t] | p_choice);
      log_lik += log_lik_t[t];

      vector[3] L;
      for (z in 1:3) {
        int is_correct_under_z = (resp_choice[t] == k_of_rule[z]);
        if (correct[t] == 1)
          L[z] = is_correct_under_z ? (1 - eta) : eta;
        else
          L[z] = is_correct_under_z ? eta : (1 - eta);
      }
      vector[3] post = b .* L;
      real Z = sum(post);
      if (Z > 0) post /= Z;
      else       post = rep_vector(1.0/3.0, 3);

      b = (1 - h) * post + (h / 2.0) * (rep_vector(1.0, 3) - post);
    }
  }
}
