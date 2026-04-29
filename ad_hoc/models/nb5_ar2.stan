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
  real<lower=0, upper=1> p1;              // partial autocorrelation lag 1 (positive)
  real<lower=-1, upper=1> p2;             // partial autocorrelation lag 2
  real beta_age;
  real<lower=0> phi;
  real<lower=0> sigma_alpha;
}

transformed parameters {
  real phi1 = p1 * (1 - p2);
  real phi2 = p2;
}

model {
  mu_s ~ normal(0, 2);
  p1 ~ beta(8, 2);                        // favors high persistence
  p2 ~ normal(0, 0.5);                    // weak prior on second lag
  sigma_alpha ~ exponential(1);
  phi ~ normal(0, 10);
  beta_age ~ normal(0, 2);

  // Initial conditions (use diffuse normals, AR(2) stationary distribution
  // is hard to write closed-form for both lags, but data dominates).
  alpha_st[, 1] ~ normal(mu_s, 2 * sigma_alpha);
  alpha_st[, 2] ~ normal(mu_s, 2 * sigma_alpha);

  for (t in 3:N_years)
    alpha_st[, t] ~ normal(
      mu_s + phi1 * (alpha_st[, t-1] - mu_s) + phi2 * (alpha_st[, t-2] - mu_s),
      sigma_alpha
    );

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
