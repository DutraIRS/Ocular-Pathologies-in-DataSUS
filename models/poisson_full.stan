data {
  int<lower=1> N;
  int<lower=1> N_states;
  int<lower=1> N_years;
  
  array[N] int<lower=1, upper=N_states> state_idx; 
  array[N] int<lower=1, upper=N_years> year_idx;
  
  vector[N] age_centered;
  vector[N] log_pop;
  array[N] int<lower=0> y;
  
  // Hyperparams
  real<lower=0> hyp_alpha;
  real<lower=0> hyp_beta;
}

parameters {
  matrix[N_states, N_years] beta_state;
  vector[N_years] beta_age;
  
  real<lower=0> sigma_state;
  real<lower=0> sigma_age;
}

model {
  # Priors for RW
  sigma_state ~ inv_gamma(hyp_alpha, hyp_beta);
  sigma_age ~ inv_gamma(hyp_alpha, hyp_beta);

  # Betas priors
  beta_state[, 1] ~ std_normal();
  beta_age[1] ~ std_normal();

  # RW
  for (t in 2:N_years) {
    beta_state[, t] ~ normal(beta_state[, t-1], sigma_state);
    beta_age[t] ~ normal(beta_age[t-1], sigma_age);
  }

  // Likelihood
  vector[N] log_lambda;
  
  for (n in 1:N) {
    int s = state_idx[n];
    int t = year_idx[n];
    
    log_lambda[n] = log_pop[n] 
                    + beta_state[s, t] 
                    + beta_age[t] * age_centered[n];
  }
  
  # force exact like calculation (needed for evidence later)
  target += poisson_log_lpmf(y | log_lambda);
}

generated quantities {
  # prevent exp overflow
  real max_safe_log_lambda = 20.7;
  
  # For LOO-CV and WAIC
  vector[N] log_lik;
  
  # For Posterior Predictive Checks
  array[N] int y_rep;
  
  for (n in 1:N) {
    int s = state_idx[n];
    int t = year_idx[n];
    real log_lambda = log_pop[n] + beta_state[s, t] + beta_age[t] * age_centered[n];

    log_lik[n] = poisson_log_lpmf(y[n] | log_lambda);
    
    # If value is too high, set to the max integer to prevent crash
    if (log_lambda > max_safe_log_lambda) {
      y_rep[n] = 2147483647; # Max 32-bit integer
    } else {
      y_rep[n] = poisson_log_rng(log_lambda);
    }
  }
}
