data {
  int<lower=1> N;                    // Number of observations
  array[N] int<lower=0> y;           // Count outcomes
  vector[N] log_pop;                 // Log population offset
  
  int<lower=1> N_states;             // Number of states
  int<lower=1> N_years;              // Number of years
  array[N] int state_idx;            // State index for each observation
  array[N] int year_idx;             // Year index for each observation
  vector[N] age_centered;            // Centered age covariate
}

parameters {
  matrix[N_states, N_years] alpha_st; // Time-varying state-specific intercepts
  matrix[N_states, N_years] beta_st; // Time-varying state-specific age effects
  real<lower=0> phi;                 // Overdispersion parameter
  real<lower=0> sigma_beta_rw;       // RW's variance for beta
  real<lower=0> sigma_alpha_rw;       // RW's variance for alpha
  
}

model {
  alpha_st[, 1] ~ normal(0, 2);
  phi ~ normal(0, 10);
  sigma_beta_rw ~ normal(0, 1);
  sigma_alpha_rw ~ normal(0, 1);
  beta_st[, 1] ~ normal(0, 2);
  for (t in 2:N_years) beta_st[, t] ~ normal(beta_st[, t-1], sigma_beta_rw);
  for (t in 2:N_years) alpha_st[, t] ~ normal(alpha_st[, t-1], sigma_alpha_rw);

  {
    vector[N] log_lambda;
    for (n in 1:N) log_lambda[n] = log_pop[n] + alpha_st[state_idx[n], year_idx[n]] + beta_st[state_idx[n], year_idx[n]] * age_centered[n];
    target += neg_binomial_2_log_lpmf(y | log_lambda, phi);
  }
}

generated quantities {
  real max_safe_log_lambda = 20.5;  // Prevent overflow in simulation
  vector[N] log_lik;                 // Pointwise log-likelihood
  array[N] int y_rep;                // Posterior predictive samples
  
  for (n in 1:N) {
    real log_lambda = log_pop[n] + alpha_st[state_idx[n], year_idx[n]] + beta_st[state_idx[n], year_idx[n]] * age_centered[n];
    
    // overflow protection
    if (log_lambda > max_safe_log_lambda) {
        log_lambda = max_safe_log_lambda;
    }
    
    // Log-likelihood for model comparison
    log_lik[n] = neg_binomial_2_log_lpmf(y[n] | log_lambda, phi);
    
    // Posterior predictive samples
    y_rep[n] = neg_binomial_2_log_rng(log_lambda, phi);
  }
}
