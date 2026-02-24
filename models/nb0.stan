data {
  int<lower=1> N;                    // Number of observations
  array[N] int<lower=0> y;           // Count outcomes
  vector[N] log_pop;                 // Log population offset
  
  // Grouping structure (not used in this baseline model)
  int<lower=1> N_states;             // Number of states
  int<lower=1> N_years;              // Number of years
  array[N] int state_idx;            // State index for each observation
  array[N] int year_idx;             // Year index for each observation
  vector[N] age_centered;            // Centered age covariate
}

parameters {
  real beta_0;                       // Global intercept
  real<lower=0> phi;                 // Overdispersion parameter
}

model {
  // Priors
  beta_0 ~ normal(0, 5);
  phi ~ normal(0, 10);
  
  // Likelihood
  target += neg_binomial_2_log_lpmf(y | log_pop + beta_0, phi);
}

generated quantities {
  real max_safe_log_lambda = 20.5;  // Prevent overflow in simulation
  vector[N] log_lik;                 // Pointwise log-likelihood
  array[N] int y_rep;                // Posterior predictive samples
  
  for (n in 1:N) {
    real log_lambda = log_pop[n] + beta_0;
    
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
