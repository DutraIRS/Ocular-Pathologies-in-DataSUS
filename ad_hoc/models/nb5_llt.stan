data {
  int<lower=1> N;
  array[N] int<lower=0> y;
  vector[N] log_pop;

  int<lower=1> N_states;
  int<lower=1> N_years;
  array[N] int state_idx;
  array[N] int year_idx;
  vector[N] age_centered;
}

parameters {
  matrix[N_states, N_years] alpha_st;   // level
  matrix[N_states, N_years] delta_st;   // trend/drift
  real beta_age;
  real<lower=0> phi;
  real<lower=0> sigma_alpha_rw;
  real<lower=0> sigma_delta_rw;
}

model {
  alpha_st[, 1] ~ normal(0, 2);
  delta_st[, 1] ~ normal(0, 0.5);
  phi ~ normal(0, 10);
  sigma_alpha_rw ~ normal(0, 1);
  sigma_delta_rw ~ normal(0, 0.5);
  beta_age ~ normal(0, 2);

  for (t in 2:N_years) {
    alpha_st[, t] ~ normal(alpha_st[, t-1] + delta_st[, t-1], sigma_alpha_rw);
    delta_st[, t] ~ normal(delta_st[, t-1], sigma_delta_rw);
  }

  {
    vector[N] log_lambda;
    for (n in 1:N)
      log_lambda[n] = log_pop[n] + alpha_st[state_idx[n], year_idx[n]]
                      + beta_age * age_centered[n];
    target += neg_binomial_2_log_lpmf(y | log_lambda, phi);
  }
}

generated quantities {
  real max_safe_log_lambda = 20.5;
  vector[N] log_lik;
  array[N] int y_rep;

  for (n in 1:N) {
    real log_lambda = log_pop[n] + alpha_st[state_idx[n], year_idx[n]]
                      + beta_age * age_centered[n];
    if (log_lambda > max_safe_log_lambda) log_lambda = max_safe_log_lambda;
    log_lik[n] = neg_binomial_2_log_lpmf(y[n] | log_lambda, phi);
    y_rep[n] = neg_binomial_2_log_rng(log_lambda, phi);
  }
}
