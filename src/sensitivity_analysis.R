library(rstan)
library(loo)
library(ggplot2)
library(dplyr)
library(bayesplot)

source("src/data_processing.R")
source("src/model_selection.R")

options(mc.cores = parallel::detectCores())

DISEASES       <- c("glaucoma", "retina", "eye")
DISEASE_LABELS <- c(glaucoma = "Glaucoma", retina = "Retinopathy",
                     eye = "Eye & Appendage Diseases")

if (!dir.exists("sensitivity")) dir.create("sensitivity", recursive = TRUE)

# 1: Phi prior sensitivity
cat("=== Part 1: Phi Prior Sensitivity ===\n")

# Read the base NB5 model code
base_code <- paste(readLines("models/nb5.stan"), collapse = "\n")

# Define alternative phi priors
phi_priors <- list(
  baseline    = "phi ~ normal(0, 10);",   # current: half-N(0,10)
  tight_hn    = "phi ~ normal(0, 3);",    # tighter half-normal
  gamma_wide  = "phi ~ gamma(2, 0.1);",   # mean=20, sd~14
  exponential = "phi ~ exponential(0.1);", # mean=10
  half_cauchy = "phi ~ cauchy(0, 5);"      # heavy-tailed, recommended by Gelman
)

build_stan_data <- function(data) {
  list(
    N           = nrow(data),
    N_states    = length(unique(data$state_id)),
    N_years     = length(unique(data$year_id)),
    state_idx   = data$state_id,
    year_idx    = data$year_id,
    age_centered = data$age_relative,
    log_pop     = data$log_pop,
    y           = data$disease
  )
}

phi_results <- list()

for (d in DISEASES) {
  cat("\n--- Disease:", d, "---\n")
  data <- prepare_data(d)
  stan_data <- build_stan_data(data)

  for (prior_name in names(phi_priors)) {
    cat("  Prior:", prior_name, "... ")

    out_dir <- file.path("sensitivity", "phi_priors", d, prior_name)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    # Substitute the phi prior line
    mod_code <- sub("phi ~ normal\\(0, 10\\);",
                    phi_priors[[prior_name]],
                    base_code, fixed = FALSE)

    model <- stan_model(model_code = mod_code,
                        model_name = paste0("nb5_", prior_name))

    fit <- sampling(
      model, data = stan_data,
      iter = 4000, chains = 4, init = 0,
      control = list(adapt_delta = 0.8, max_treedepth = 10)
    )

    # --- Extract phi posterior ---
    phi_samp <- as.numeric(rstan::extract(fit, "phi")$phi)

    # --- LOO-IC ---
    ll <- extract_log_lik(fit, merge_chains = FALSE)
    r  <- relative_eff(exp(ll))
    loo_res <- loo(ll, r_eff = r, cores = 1)

    # --- Convergence ---
    ess_rhat <- check_ess_rhat(fit, out_dir)

    # --- PPC plots ---
    plot_ppc(fit, stan_data$y, out_dir)

    # --- Beta_age and sigma_rw posteriors (key parameters) ---
    beta_samp  <- as.numeric(rstan::extract(fit, "beta_age")$beta_age)
    sigma_samp <- as.numeric(rstan::extract(fit, "sigma_alpha_rw")$sigma_alpha_rw)

    saveRDS(fit, file.path(out_dir, "fit.rds"))

    phi_results[[length(phi_results) + 1]] <- data.frame(
      disease        = d,
      disease_label  = DISEASE_LABELS[d],
      prior          = prior_name,
      prior_spec     = phi_priors[[prior_name]],
      phi_mean       = round(mean(phi_samp), 3),
      phi_median     = round(median(phi_samp), 3),
      phi_sd         = round(sd(phi_samp), 3),
      phi_q025       = round(quantile(phi_samp, 0.025), 3),
      phi_q975       = round(quantile(phi_samp, 0.975), 3),
      beta_age_mean  = round(mean(beta_samp), 4),
      beta_age_sd    = round(sd(beta_samp), 4),
      sigma_rw_mean  = round(mean(sigma_samp), 4),
      sigma_rw_sd    = round(sd(sigma_samp), 4),
      looic          = round(loo_res$estimates["looic", "Estimate"], 1),
      looic_se       = round(loo_res$estimates["looic", "SE"], 1),
      max_rhat       = round(ess_rhat["max_rhat"], 4),
      min_n_eff      = round(ess_rhat["min_n_eff"], 0),
      stringsAsFactors = FALSE
    )

    rm(fit, model); gc()
    cat("done.\n")
  }
}

phi_table <- do.call(rbind, phi_results)
rownames(phi_table) <- NULL
write.csv(phi_table, "sensitivity/phi_prior_sensitivity.csv", row.names = FALSE)
cat("\nPhi sensitivity table saved.\n")

for (d in DISEASES) {
  phi_draws <- list()
  for (pn in names(phi_priors)) {
    fp <- file.path("sensitivity", "phi_priors", d, pn, "fit.rds")
    if (!file.exists(fp)) next
    fit <- readRDS(fp)
    phi_draws[[pn]] <- data.frame(
      prior = pn,
      phi   = as.numeric(rstan::extract(fit, "phi")$phi)
    )
    rm(fit); gc()
  }

  if (length(phi_draws) == 0) next
  phi_df <- do.call(rbind, phi_draws)
  phi_df$prior <- factor(phi_df$prior, levels = names(phi_priors))

  p <- ggplot(phi_df, aes(x = phi, color = prior, fill = prior)) +
    geom_density(alpha = 0.2, linewidth = 0.7) +
    labs(
      title    = paste0("Phi Posterior Under Alternative Priors: ",
                        DISEASE_LABELS[d]),
      subtitle = "If posteriors largely overlap, inference is robust to prior choice.",
      x = expression(phi ~ "(overdispersion)"),
      y = "Posterior Density"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(file.path("sensitivity", "phi_priors", d,
                    "phi_posterior_comparison.png"),
         p, width = 10, height = 6, dpi = 300)
  write.csv(phi_df,
            file.path("sensitivity", "phi_priors", d,
                       "phi_posterior_comparison.csv"),
            row.names = FALSE)
}

cat("Phi posterior comparison plots saved.\n")


# 2: Excluding 2024 to check for data completeness issues
cat("\n=== Part 2: Excluding 2024 ===\n")

# Modify prepare_data to also exclude 2024
prepare_data_no2024 <- function(disease = "glaucoma") {
  data <- read.csv("data/sus_ocular_data.csv")
  data[data == -1] <- NA
  data$state <- as.factor(data$state)
  data$age_relative <- (data$age - 60) / 10
  data <- subset(data, year != 2025 & year != 2024)
  data$year_relative <- data$year - min(data$year)
  data$log_pop <- log(data$population)
  data$disease <- data[[disease]]
  data <- data[c("state", "year_relative", "age_relative", "log_pop", "disease")]
  data <- na.omit(data)
  data$state_id <- as.integer(as.factor(data$state))
  data$year_id  <- data$year_relative + 1
  return(data)
}

year_results <- list()

model_nb5 <- stan_model(file = "models/nb5.stan", model_name = "nb5")

for (d in DISEASES) {
  cat("  Disease:", d, "... ")

  out_dir <- file.path("sensitivity", "exclude_2024", d)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  data <- prepare_data_no2024(d)
  stan_data <- build_stan_data(data)

  fit <- sampling(model_nb5, data = stan_data,
                  iter = 4000, chains = 4, init = 0,
                  control = list(adapt_delta = 0.8, max_treedepth = 10))

  phi_samp   <- as.numeric(rstan::extract(fit, "phi")$phi)
  beta_samp  <- as.numeric(rstan::extract(fit, "beta_age")$beta_age)
  sigma_samp <- as.numeric(rstan::extract(fit, "sigma_alpha_rw")$sigma_alpha_rw)

  ll <- extract_log_lik(fit, merge_chains = FALSE)
  r  <- relative_eff(exp(ll))
  loo_res <- loo(ll, r_eff = r, cores = 1)
  ess_rhat <- check_ess_rhat(fit, out_dir)
  plot_ppc(fit, stan_data$y, out_dir)
  saveRDS(fit, file.path(out_dir, "fit.rds"))

  year_results[[length(year_results) + 1]] <- data.frame(
    disease       = d,
    disease_label = DISEASE_LABELS[d],
    scenario      = "exclude_2024",
    phi_mean      = round(mean(phi_samp), 3),
    beta_age_mean = round(mean(beta_samp), 4),
    sigma_rw_mean = round(mean(sigma_samp), 4),
    looic         = round(loo_res$estimates["looic", "Estimate"], 1),
    max_rhat      = round(ess_rhat["max_rhat"], 4),
    min_n_eff     = round(ess_rhat["min_n_eff"], 0),
    stringsAsFactors = FALSE
  )

  rm(fit); gc()
  cat("done.\n")
}

year_table <- do.call(rbind, year_results)
write.csv(year_table, "sensitivity/exclude_2024_results.csv", row.names = FALSE)
cat("Exclude-2024 results saved.\n")


# 3: Treating missing data as zero instead of NA
cat("\n=== Part 3: Missing Data as Zero ===\n")

prepare_data_zero <- function(disease = "glaucoma") {
  data <- read.csv("data/sus_ocular_data.csv")
  data[data == -1] <- 0   # treat suppressed cells as zero
  data$state <- as.factor(data$state)
  data$age_relative <- (data$age - 60) / 10
  data <- subset(data, year != 2025)
  data$year_relative <- data$year - min(data$year)
  data$log_pop <- log(data$population)
  data$disease <- data[[disease]]
  data <- data[c("state", "year_relative", "age_relative", "log_pop", "disease")]
  data <- na.omit(data)
  data$state_id <- as.integer(as.factor(data$state))
  data$year_id  <- data$year_relative + 1
  return(data)
}

miss_results <- list()

for (d in DISEASES) {
  cat("  Disease:", d, "... ")

  out_dir <- file.path("sensitivity", "missing_as_zero", d)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  data <- prepare_data_zero(d)
  stan_data <- build_stan_data(data)

  fit <- sampling(model_nb5, data = stan_data,
                  iter = 4000, chains = 4, init = 0,
                  control = list(adapt_delta = 0.8, max_treedepth = 10))

  phi_samp   <- as.numeric(rstan::extract(fit, "phi")$phi)
  beta_samp  <- as.numeric(rstan::extract(fit, "beta_age")$beta_age)
  sigma_samp <- as.numeric(rstan::extract(fit, "sigma_alpha_rw")$sigma_alpha_rw)

  ll <- extract_log_lik(fit, merge_chains = FALSE)
  r  <- relative_eff(exp(ll))
  loo_res <- loo(ll, r_eff = r, cores = 1)
  ess_rhat <- check_ess_rhat(fit, out_dir)
  plot_ppc(fit, stan_data$y, out_dir)
  saveRDS(fit, file.path(out_dir, "fit.rds"))

  miss_results[[length(miss_results) + 1]] <- data.frame(
    disease       = d,
    disease_label = DISEASE_LABELS[d],
    scenario      = "missing_as_zero",
    phi_mean      = round(mean(phi_samp), 3),
    beta_age_mean = round(mean(beta_samp), 4),
    sigma_rw_mean = round(mean(sigma_samp), 4),
    looic         = round(loo_res$estimates["looic", "Estimate"], 1),
    max_rhat      = round(ess_rhat["max_rhat"], 4),
    min_n_eff     = round(ess_rhat["min_n_eff"], 0),
    stringsAsFactors = FALSE
  )

  rm(fit); gc()
  cat("done.\n")
}

miss_table <- do.call(rbind, miss_results)
write.csv(miss_table, "sensitivity/missing_as_zero_results.csv", row.names = FALSE)
cat("Missing-as-zero results saved.\n")


# summary
cat("\n=== Generating Combined Summary ===\n")

# Load baseline results from diagnostics for comparison
baseline_results <- list()
for (d in DISEASES) {
  comp <- read.csv(file.path("diagnostics", d, "model_comparison.csv"),
                   row.names = 1)
  nb5_row <- comp["nb5", ]

  # Load fit to get parameter estimates
  fit_path <- file.path("diagnostics", d, "nb5", "fit.rds")
  if (file.exists(fit_path)) {
    fit <- readRDS(fit_path)
    phi_samp   <- as.numeric(rstan::extract(fit, "phi")$phi)
    beta_samp  <- as.numeric(rstan::extract(fit, "beta_age")$beta_age)
    sigma_samp <- as.numeric(rstan::extract(fit, "sigma_alpha_rw")$sigma_alpha_rw)
    rm(fit); gc()

    baseline_results[[d]] <- data.frame(
      disease       = d,
      disease_label = DISEASE_LABELS[d],
      scenario      = "baseline",
      phi_mean      = round(mean(phi_samp), 3),
      beta_age_mean = round(mean(beta_samp), 4),
      sigma_rw_mean = round(mean(sigma_samp), 4),
      looic         = round(nb5_row$looic, 1),
      max_rhat      = round(nb5_row$max_rhat, 4),
      min_n_eff     = round(nb5_row$min_n_eff, 0),
      stringsAsFactors = FALSE
    )
  }
}

baseline_table <- do.call(rbind, baseline_results)
combined <- rbind(baseline_table, year_table, miss_table)
write.csv(combined, "sensitivity/scenario_comparison.csv", row.names = FALSE)

cat("\n=== Sensitivity analysis complete ===\n")
cat("All outputs in sensitivity/:\n")
cat("  phi_prior_sensitivity.csv      - phi prior robustness\n")
cat("  phi_priors/{disease}/           - fits, PPC, posterior plots\n")
cat("  exclude_2024_results.csv       - year completeness check\n")
cat("  missing_as_zero_results.csv    - alternative missing data\n")
cat("  scenario_comparison.csv        - combined baseline vs scenarios\n")
