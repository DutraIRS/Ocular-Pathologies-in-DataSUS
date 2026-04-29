data {
  int<lower=1> N;
  int<lower=1> N_age_groups;
  array[N] int<lower=0> y;
  vector[N] log_pop;

  int<lower=1> N_states;
  int<lower=1> N_years;
  array[N] int state_idx;
  array[N] int year_idx;
  array[N] int<lower=1, upper=7> age_idx;
}

parameters {
  matrix[N_states, N_years] alpha_st;
  vector[N_age_groups - 1] beta_age_free;     // first K-1 free
  real<lower=0> phi;
  real<lower=0> sigma_alpha_rw;
}

transformed parameters {
  vector[N_age_groups] beta_age_cat;
  for (j in 1:(N_age_groups - 1))
    beta_age_cat[j] = beta_age_free[j];
  beta_age_cat[N_age_groups] = -sum(beta_age_free);   // sum-to-zero
}

model {
  alpha_st[, 1] ~ normal(0, 2);
  phi ~ normal(0, 10);
  sigma_alpha_rw ~ normal(0, 1);
  beta_age_free ~ normal(0, 2);

  for (t in 2:N_years)
    alpha_st[, t] ~ normal(alpha_st[, t-1], sigma_alpha_rw);

  {
    vector[N] log_lambda;
    for (n in 1:N)
      log_lambda[n] = log_pop[n] + alpha_st[state_idx[n], year_idx[n]]
                      + beta_age_cat[age_idx[n]];
    target += neg_binomial_2_log_lpmf(y | log_lambda, phi);
  }
}

generated quantities {
  real max_safe_log_lambda = 20.5;
  vector[N] log_lik;
  array[N] int y_rep;

  for (n in 1:N) {
    real log_lambda = log_pop[n] + alpha_st[state_idx[n], year_idx[n]]
                      + beta_age_cat[age_idx[n]];
    if (log_lambda > max_safe_log_lambda) log_lambda = max_safe_log_lambda;
    log_lik[n] = neg_binomial_2_log_lpmf(y[n] | log_lambda, phi);
    y_rep[n] = neg_binomial_2_log_rng(log_lambda, phi);
  }
}
