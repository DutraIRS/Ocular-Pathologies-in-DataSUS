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
  matrix[N_states, N_years] alpha_st;
  vector[N_states] mu_s;
  real mu_global;
  real<lower=0> sigma_mu;
  real<lower=0, upper=1> rho;
  real beta_age;
  real<lower=0> phi;
  real<lower=0> sigma_alpha;
}

model {
  mu_global ~ normal(0, 2);
  sigma_mu ~ exponential(1);
  mu_s ~ normal(mu_global, sigma_mu);
  rho ~ beta(8, 2);
  sigma_alpha ~ exponential(1);
  phi ~ normal(0, 10);
  beta_age ~ normal(0, 2);

  alpha_st[, 1] ~ normal(mu_s, sigma_alpha / sqrt(1 - square(rho)));

  for (t in 2:N_years)
    alpha_st[, t] ~ normal(mu_s + rho * (alpha_st[, t-1] - mu_s), sigma_alpha);

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
