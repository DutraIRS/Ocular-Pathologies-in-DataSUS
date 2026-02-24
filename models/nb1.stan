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
  vector[N_states] alpha;            // State-specific intercepts
  real beta_age;                     // Global age effect
  real<lower=0> phi;                 // Overdispersion parameter
}

model {
  // Priors
  alpha ~ normal(0, 2);
  beta_age ~ normal(0, 2);
  phi ~ normal(0, 10);
  
  // Likelihood
  {
    vector[N] log_lambda;
    for (n in 1:N) {
      log_lambda[n] = log_pop[n] + alpha[state_idx[n]] + beta_age * age_centered[n];
    }
    target += neg_binomial_2_log_lpmf(y | log_lambda, phi);
  }
}

generated quantities {
  real max_safe_log_lambda = 20.5;  // Prevent overflow in simulation
  vector[N] log_lik;                 // Pointwise log-likelihood
  array[N] int y_rep;                // Posterior predictive samples
  
  for (n in 1:N) {
    real log_lambda = log_pop[n] + alpha[state_idx[n]] + beta_age * age_centered[n];
    
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
