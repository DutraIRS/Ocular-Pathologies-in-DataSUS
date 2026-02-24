library(loo)
library(rstan)
library(ggplot2)
library(bayesplot)
library(bridgesampling)

calc_loo_waic <- function(fit, output_dir) {
  log_lik <- extract_log_lik(fit, merge_chains = FALSE)
  r_eff <- relative_eff(exp(log_lik))
  
  loo_result <- loo(log_lik, r_eff = r_eff, cores = 1)
  waic_result <- waic(log_lik, cores = 1)
  
  sink(file.path(output_dir, "loo_waic_summary.txt"))
  print(loo_result)
  print(waic_result)
  sink()
  
  return(c(
    looic = loo_result$estimates["looic", "Estimate"],
    looic_se = loo_result$estimates["looic", "SE"],
    waic = waic_result$estimates["waic", "Estimate"],
    waic_se = waic_result$estimates["waic", "SE"]
  ))
}

check_ess_rhat <- function(fit, output_dir) {
  summary_fit <- summary(fit)$summary
  df_summary <- as.data.frame(summary_fit)
  df_summary$parameter <- rownames(summary_fit)
  
  # 1. Filter out the nuisance constant
  # We also remove 'lp__' (log probability) which often has lower ESS but isn't a parameter
  df_clean <- df_summary[!(df_summary$parameter %in% c("max_safe_log_lambda", "lp__")), ]
  
  # Save the FULL table (including the constant) to CSV for reference
  write.csv(df_summary, file.path(output_dir, "summary_stats.csv"), row.names = FALSE)
  
  # 2. Calculate stats only on the CLEAN data
  return(c(
    max_rhat = max(df_clean$Rhat, na.rm = TRUE),
    min_n_eff = min(df_clean$n_eff, na.rm = TRUE)
  ))
}

calc_marginal_likelihood <- function(fit, output_dir) {
  bridge <- bridge_sampler(fit, verbose = FALSE, silent = TRUE)

  err <- error_measures(bridge)
  log_ml_se <- err$re
  
  return(c(
    log_ml = bridge$logml, 
    log_ml_se = log_ml_se
  ))
}

plot_mcmc_diagnostics <- function(fit, output_dir, pars) {
  posterior <- as.array(fit)
  
  p_trace <- mcmc_trace(posterior, regex_pars = pars) 
  ggsave(file.path(output_dir, "trace.png"), p_trace, width = 10, height = 6, dpi = 300)
  
  p_acf <- mcmc_acf(posterior, regex_pars = pars)
  ggsave(file.path(output_dir, "acf.png"), p_acf, width = 10, height = 6, dpi = 300)
}

plot_ppc <- function(fit, data_y, output_dir) {
  y_rep <- as.matrix(fit, pars = "y_rep")
  
  p_dens <- ppc_dens_overlay(data_y, y_rep[1:50, ]) + 
    ggtitle("PPC: Density Overlay")
  ggsave(file.path(output_dir, "ppc_density.png"), p_dens, width = 8, height = 6, dpi = 300)
  
  p_stat <- ppc_stat(data_y, y_rep, stat = "mean") + 
    ggtitle("PPC: Mean Check")
  ggsave(file.path(output_dir, "ppc_mean.png"), p_stat, width = 8, height = 6, dpi = 300)
  
  p_sd <- ppc_stat(data_y, y_rep, stat = "sd") + 
    ggtitle("PPC: Standard Deviation (Overdispersion Check)")
  ggsave(file.path(output_dir, "ppc_stat_sd.png"), p_sd, width = 8, height = 6, dpi = 300)
  
  prop_zero <- function(x) mean(x == 0)
  
  p_zero <- ppc_stat(data_y, y_rep, stat = prop_zero) + 
    ggtitle("PPC: Proportion of Zeros")
  ggsave(file.path(output_dir, "ppc_stat_zeros.png"), p_zero, width = 8, height = 6, dpi = 300)
  
  p_max <- ppc_stat(data_y, y_rep, stat = "max") + 
    ggtitle("PPC: Maximum Value")
  ggsave(file.path(output_dir, "ppc_stat_max.png"), p_max, width = 8, height = 6, dpi = 300)
}

run_all_diagnostics <- function(fit, output_dir, data_y, pars = "beta") {
  res_loo <- calc_loo_waic(fit, output_dir)
  res_ess <- check_ess_rhat(fit, output_dir)
  res_ml  <- calc_marginal_likelihood(fit, output_dir)
  
  plot_mcmc_diagnostics(fit, output_dir, pars)
  plot_ppc(fit, data_y, output_dir)
  
  metrics <- c(res_loo, res_ml, res_ess)
  return(metrics)
}