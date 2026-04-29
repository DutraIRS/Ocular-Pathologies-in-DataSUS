library(rstan)
library(loo)
library(ggplot2)
library(dplyr)
library(tidyr)

source("src/data_processing.R")
source("src/model_selection.R")

options(mc.cores = parallel::detectCores() - 2)
rstan_options(auto_write = TRUE)

MODEL_VARIANTS <- c(
  "nb5", "nb5_expo", "nb5_llt", "nb5_llt_expo",
  "nb5_drift", "nb5_drift_expo",
  "nb5_llt_drift",
  "nb5_llt_drift_expo",
  # Narrower-interval experiments
  "nb5_tight_sigma",   # sigma ~ exp(5), no drift
  "nb5_drift_tight",   # drift + sigma ~ exp(5)
  "nb5_llt_damped",    # damped LLT (trend decays -> bounded forecast)
  "nb5_ar1",           # AR(1) stationary (forecast reverts to state mean)
  # AR family extensions
  "nb5_ar1_hier",      # AR(1) + hierarchical pooling on mu_s
  "nb5_ar1_trend",     # AR(1) around a deterministic linear trend
  "nb5_ar2",           # AR(2), stationary via PACF parameterization
  # Age-effect variants
  "nb5_age_quad",      # NB5 + quadratic age
  "nb5_age_cat",       # NB5 + categorical age (sum-to-zero deflections)
  "nb5_ar1_age_quad",  # AR(1) + quadratic age
  "nb5_ar1_age_cat"    # AR(1) + categorical age
)
DISEASES       <- c("glaucoma", "retina", "eye")
DISEASE_LABELS <- c(glaucoma = "Glaucoma", retina = "Retinopathy",
                     eye = "Eye & Appendage Diseases")
FORECAST_END   <- 2036
N_POST_DRAWS   <- 500
FORCE_REFIT    <- FALSE   # set TRUE to re-fit all models

AGE_GROUPS    <- c(25, 35, 45, 55, 65, 75, 90)
AGE_LABELS    <- c("25-34", "35-44", "45-54", "55-64",
                    "65-74", "75-89", "90+")
AGE_LABEL_MAP <- setNames(AGE_LABELS, AGE_GROUPS)

REGION_MAP <- c(
  RO = "North", AC = "North", AM = "North", RR = "North",
  PA = "North", AP = "North", TO = "North",
  MA = "Northeast", PI = "Northeast", CE = "Northeast", RN = "Northeast",
  PB = "Northeast", PE = "Northeast", AL = "Northeast", SE = "Northeast",
  BA = "Northeast",
  MG = "Southeast", ES = "Southeast", RJ = "Southeast", SP = "Southeast",
  PR = "South", SC = "South", RS = "South",
  MS = "Center-West", MT = "Center-West", GO = "Center-West",
  DF = "Center-West"
)
REGION_ORDER <- c("North", "Northeast", "Center-West", "Southeast", "South")

STATE_ORDER <- c(
  "AC", "AM", "AP", "PA", "RO", "RR", "TO",
  "AL", "BA", "CE", "MA", "PB", "PE", "PI", "RN", "SE",
  "DF", "GO", "MS", "MT",
  "ES", "MG", "RJ", "SP",
  "PR", "RS", "SC"
)

MODELS_DIR  <- "ad_hoc/models"
RESULTS_DIR <- "ad_hoc/results"
IBGE_PATH   <- "data/ibge_population_projections.csv"

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle    = element_text(size = base_size, color = "gray40", hjust = 0),
      strip.text       = element_text(face = "bold"),
      legend.position  = "bottom",
      panel.grid.minor = element_blank()
    )
}

raw <- read.csv("data/sus_ocular_data.csv")
raw[raw == -1] <- NA
raw <- subset(raw, year != 2025)
raw$region    <- factor(REGION_MAP[raw$state], levels = REGION_ORDER)
raw$age_label <- factor(AGE_LABEL_MAP[as.character(raw$age)], levels = AGE_LABELS)

min_year <- min(raw$year)
max_year <- max(raw$year)

# Variant structure detection
is_llt       <- function(variant) grepl("llt", variant)
has_drift    <- function(variant) grepl("drift", variant) && !grepl("damped", variant) && !grepl("ar1_trend", variant)
is_damped    <- function(variant) grepl("damped", variant)
# Plain AR(1) family — uses the same forecast logic (mu_s + rho)
is_ar1_plain <- function(variant) variant %in% c("nb5_ar1", "nb5_ar1_hier",
                                                  "nb5_ar1_age_quad", "nb5_ar1_age_cat")
is_ar1_trend <- function(variant) variant == "nb5_ar1_trend"
is_ar2       <- function(variant) variant == "nb5_ar2"
is_ar_family <- function(variant) is_ar1_plain(variant) || is_ar1_trend(variant) || is_ar2(variant)
# Age-effect detection
is_quad_age  <- function(variant) grepl("age_quad", variant)
is_cat_age   <- function(variant) grepl("age_cat", variant)

build_stan_data <- function(data, variant = NULL) {
  # age_idx: 1..7 corresponding to age groups 25, 35, 45, 55, 65, 75, 90
  age_levels <- c(-3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 3.0)
  age_idx <- match(round(data$age_relative, 1), round(age_levels, 1))

  ds <- list(
    N            = nrow(data),
    N_states     = length(unique(data$state_id)),
    N_years      = length(unique(data$year_id)),
    state_idx    = data$state_id,
    year_idx     = data$year_id,
    age_centered = data$age_relative,
    log_pop      = data$log_pop,
    y            = data$disease
  )
  # Categorical variants need age_idx and N_age_groups
  if (!is.null(variant) && is_cat_age(variant)) {
    ds$age_idx       <- age_idx
    ds$N_age_groups  <- length(AGE_GROUPS)
  }
  ds
}

compute_metrics <- function(y_obs, y_rep_sub, pop) {
  # y_rep_sub: [n_draws, n_obs]
  pp_mean <- colMeans(y_rep_sub)

  # Rate-scale metrics (per million)
  obs_rate  <- y_obs / pop * 1e6
  pred_rate <- pp_mean / pop * 1e6

  mae  <- mean(abs(obs_rate - pred_rate))
  rmse <- sqrt(mean((obs_rate - pred_rate)^2))
  bias <- mean(pred_rate - obs_rate)

  # Mean relative error (safe denominator: rate + 1 to avoid div by zero)
  mre <- mean(abs(obs_rate - pred_rate) / pmax(obs_rate, 1))

  # Coverage
  q025 <- apply(y_rep_sub, 2, quantile, 0.025)
  q975 <- apply(y_rep_sub, 2, quantile, 0.975)
  q10  <- apply(y_rep_sub, 2, quantile, 0.10)
  q90  <- apply(y_rep_sub, 2, quantile, 0.90)
  cov80 <- mean(y_obs >= q10 & y_obs <= q90)
  cov95 <- mean(y_obs >= q025 & y_obs <= q975)

  data.frame(
    N = length(y_obs),
    MAE = round(mae, 3),
    MRE = round(mre, 4),
    RMSE = round(rmse, 3),
    Bias = round(bias, 3),
    Cov_80 = round(cov80 * 100, 1),
    Cov_95 = round(cov95 * 100, 1)
  )
}

forecast_model <- function(fit, variant, n_forecast, n_use) {
  post <- rstan::extract(fit)
  alpha_full <- post$alpha_st
  n_draws <- dim(alpha_full)[1]
  n_states <- dim(alpha_full)[2]
  n_years  <- dim(alpha_full)[3]

  di <- sample(n_draws, min(n_draws, n_use))
  n_use <- length(di)

  alpha_st <- alpha_full[di, , , drop = FALSE]

  # Compute age_term: [n_use, N_age_groups] additive contribution by age
  # (handles linear, quadratic, and categorical age structures uniformly)
  age_x <- (AGE_GROUPS - 60) / 10
  if (is_cat_age(variant)) {
    age_term <- post$beta_age_cat[di, , drop = FALSE]
  } else if (is_quad_age(variant)) {
    b1 <- as.numeric(post$beta_age[di])
    b2 <- as.numeric(post$beta_age2[di])
    age_term <- outer(b1, age_x) + outer(b2, age_x^2)
  } else {
    b1 <- as.numeric(post$beta_age[di])
    age_term <- outer(b1, age_x)
  }

  # AR family uses sigma_alpha (no _rw suffix); RW/LLT use sigma_alpha_rw
  sigma_alpha <- if (is_ar_family(variant))
    as.numeric(post$sigma_alpha[di])
  else
    as.numeric(post$sigma_alpha_rw[di])

  alpha_last <- alpha_st[, , n_years]
  alpha_fwd  <- array(NA, dim = c(n_use, n_states, n_forecast))

  if (is_ar1_plain(variant)) {
    rho_draws <- as.numeric(post$rho[di])
    mu_s_mat  <- post$mu_s[di, , drop = FALSE]   # [n_use, n_states]
    for (h in seq_len(n_forecast)) {
      eps <- sweep(matrix(rnorm(n_use * n_states), n_use, n_states),
                   1, sigma_alpha, "*")
      centered <- alpha_last - mu_s_mat
      centered <- sweep(centered, 1, rho_draws, "*")
      alpha_last <- mu_s_mat + centered + eps
      alpha_fwd[, , h] <- alpha_last
    }
    return(list(alpha_st = alpha_st, alpha_fwd = alpha_fwd,
                age_term = age_term, di = di, n_use = n_use))
  }

  if (is_ar1_trend(variant)) {
    rho_draws   <- as.numeric(post$rho[di])
    delta_draws <- as.numeric(post$delta[di])
    mu_s_mat    <- post$mu_s[di, , drop = FALSE]
    # Last observed time index in original (1-based) coordinates
    t_last <- n_years
    # trend at t_last (column-recycle delta_draws across states):
    # trend_prev[i, s] = mu_s_mat[i, s] + delta_draws[i] * (t - 1)
    for (h in seq_len(n_forecast)) {
      t_now <- t_last + h
      trend_prev <- mu_s_mat + matrix(delta_draws * (t_now - 2),
                                      nrow = n_use, ncol = n_states)
      trend_now  <- mu_s_mat + matrix(delta_draws * (t_now - 1),
                                      nrow = n_use, ncol = n_states)
      eps <- sweep(matrix(rnorm(n_use * n_states), n_use, n_states),
                   1, sigma_alpha, "*")
      centered <- alpha_last - trend_prev
      centered <- sweep(centered, 1, rho_draws, "*")
      alpha_last <- trend_now + centered + eps
      alpha_fwd[, , h] <- alpha_last
    }
    return(list(alpha_st = alpha_st, alpha_fwd = alpha_fwd,
                age_term = age_term, di = di, n_use = n_use))
  }

  if (is_ar2(variant)) {
    phi1_draws <- as.numeric(post$phi1[di])
    phi2_draws <- as.numeric(post$phi2[di])
    mu_s_mat   <- post$mu_s[di, , drop = FALSE]
    # alpha_last is alpha[T] (last obs); we also need alpha[T-1]
    alpha_lag1 <- alpha_last
    alpha_lag2 <- alpha_st[, , n_years - 1]
    for (h in seq_len(n_forecast)) {
      eps <- sweep(matrix(rnorm(n_use * n_states), n_use, n_states),
                   1, sigma_alpha, "*")
      c1 <- sweep(alpha_lag1 - mu_s_mat, 1, phi1_draws, "*")
      c2 <- sweep(alpha_lag2 - mu_s_mat, 1, phi2_draws, "*")
      alpha_new <- mu_s_mat + c1 + c2 + eps
      # Shift lags
      alpha_lag2 <- alpha_lag1
      alpha_lag1 <- alpha_new
      alpha_fwd[, , h] <- alpha_new
    }
    return(list(alpha_st = alpha_st, alpha_fwd = alpha_fwd,
                age_term = age_term, di = di, n_use = n_use))
  }

  if (is_damped(variant)) {
    phi_damp <- as.numeric(post$phi_damp[di])
    delta_full  <- post$delta_st
    sigma_delta <- as.numeric(post$sigma_delta_rw[di])
    delta_last  <- delta_full[di, , n_years]
    for (h in seq_len(n_forecast)) {
      eps_a <- sweep(matrix(rnorm(n_use * n_states), n_use, n_states),
                     1, sigma_alpha, "*")
      eps_d <- sweep(matrix(rnorm(n_use * n_states), n_use, n_states),
                     1, sigma_delta, "*")
      damped_delta <- sweep(delta_last, 1, phi_damp, "*")
      alpha_last   <- alpha_last + damped_delta + eps_a
      delta_last   <- damped_delta + eps_d
      alpha_fwd[, , h] <- alpha_last
    }
    return(list(alpha_st = alpha_st, alpha_fwd = alpha_fwd,
                age_term = age_term, di = di, n_use = n_use))
  }

  mu_draws <- if (has_drift(variant)) as.numeric(post$mu[di]) else NULL

  if (is_llt(variant)) {
    delta_full  <- post$delta_st
    sigma_delta <- as.numeric(post$sigma_delta_rw[di])
    delta_last  <- delta_full[di, , n_years]
  }

  for (h in seq_len(n_forecast)) {
    eps_a <- sweep(matrix(rnorm(n_use * n_states), n_use, n_states),
                   1, sigma_alpha, "*")
    if (is_llt(variant)) {
      eps_d <- sweep(matrix(rnorm(n_use * n_states), n_use, n_states),
                     1, sigma_delta, "*")
      alpha_last <- alpha_last + delta_last + eps_a
      delta_last <- delta_last + eps_d
      if (has_drift(variant)) {
        delta_last <- sweep(delta_last, 1, mu_draws, "+")
      }
    } else {
      alpha_last <- alpha_last + eps_a
      if (has_drift(variant)) {
        alpha_last <- sweep(alpha_last, 1, mu_draws, "+")
      }
    }
    alpha_fwd[, , h] <- alpha_last
  }

  list(alpha_st = alpha_st, alpha_fwd = alpha_fwd,
       age_term = age_term, di = di, n_use = n_use)
}

pop_obs <- raw %>% select(state, year, age, population) %>% distinct()
forecast_years <- (max_year + 1):FORECAST_END

if (file.exists(IBGE_PATH)) {
  pop_ibge <- read.csv(IBGE_PATH)
  names(pop_ibge) <- c("state", "year", "age", "population")
  pop_future <- pop_ibge %>% filter(year %in% forecast_years)
} else {
  stop("IBGE projections not found. Run src/fetch_ibge_projections.R first.")
}
pop_all <- rbind(pop_obs, pop_future)
all_years <- min_year:FORECAST_END
age_cent_vals <- (AGE_GROUPS - 60) / 10

overall_comparison <- list()
subgroup_comparison <- list()

for (variant in MODEL_VARIANTS) {
  cat("\n\n############################################################\n")
  cat("MODEL VARIANT:", variant, "\n")
  cat("############################################################\n")

  for (d in DISEASES) {
    cat("\n---- Disease:", d, "----\n")

    out_dir <- file.path(RESULTS_DIR, variant, d)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    fit_path <- file.path(out_dir, "fit.rds")

    data <- prepare_data(d)
    stan_data <- build_stan_data(data, variant)

    if (!FORCE_REFIT && file.exists(fit_path)) {
      cat("  Loading existing fit...\n")
      fit <- readRDS(fit_path)
    } else {
      cat("  Compiling & sampling...\n")
      model <- stan_model(file = file.path(MODELS_DIR, paste0(variant, ".stan")),
                          model_name = variant)
      # LLT, damped, and all AR-family models benefit from higher adapt_delta
      ad <- if (is_llt(variant) || is_damped(variant) || is_ar_family(variant)) 0.95 else 0.8
      fit <- sampling(
        model, data = stan_data,
        iter = 4000, chains = 4, init = 0,
        control = list(adapt_delta = ad, max_treedepth = 12)
      )
      saveRDS(fit, fit_path)
      rm(model)
    }

    cat("  Running diagnostics...\n")
    ess_rhat <- check_ess_rhat(fit, out_dir)
    loo_wa   <- calc_loo_waic(fit, out_dir)
    plot_ppc(fit, stan_data$y, out_dir)

    y_rep <- as.matrix(fit, pars = "y_rep")
    n_use <- min(nrow(y_rep), N_POST_DRAWS)
    y_rep <- y_rep[sample(nrow(y_rep), n_use), , drop = FALSE]

    # Metadata for subgroups
    meta <- data.frame(
      y_obs = data$disease,
      state = as.character(data$state),
      age   = round(data$age_relative * 10 + 60),
      pop   = exp(data$log_pop),
      stringsAsFactors = FALSE
    )
    meta$age_label <- factor(AGE_LABEL_MAP[as.character(meta$age)],
                             levels = AGE_LABELS)
    meta$region    <- factor(REGION_MAP[meta$state], levels = REGION_ORDER)

    # Overall
    overall_m <- compute_metrics(meta$y_obs, y_rep, meta$pop)
    overall_comparison[[length(overall_comparison) + 1]] <- cbind(
      variant = variant,
      disease = DISEASE_LABELS[d],
      overall_m,
      LOOIC     = round(loo_wa["looic"], 1),
      LOOIC_SE  = round(loo_wa["looic_se"], 1),
      WAIC      = round(loo_wa["waic"], 1),
      max_rhat  = round(ess_rhat["max_rhat"], 4),
      min_n_eff = round(ess_rhat["min_n_eff"], 0)
    )

    # By age group
    for (ag in AGE_LABELS) {
      ii <- which(meta$age_label == ag)
      if (length(ii) == 0) next
      m <- compute_metrics(meta$y_obs[ii], y_rep[, ii, drop = FALSE], meta$pop[ii])
      subgroup_comparison[[length(subgroup_comparison) + 1]] <- cbind(
        variant  = variant,
        disease  = DISEASE_LABELS[d],
        group    = "Age",
        subgroup = ag,
        m
      )
    }
    # By region
    for (rg in REGION_ORDER) {
      ii <- which(meta$region == rg)
      if (length(ii) == 0) next
      m <- compute_metrics(meta$y_obs[ii], y_rep[, ii, drop = FALSE], meta$pop[ii])
      subgroup_comparison[[length(subgroup_comparison) + 1]] <- cbind(
        variant  = variant,
        disease  = DISEASE_LABELS[d],
        group    = "Region",
        subgroup = rg,
        m
      )
    }

    cat("  Generating fitted chart...\n")

    # Per-observation predicted mean
    pp_mean <- colMeans(y_rep)
    q025 <- apply(y_rep, 2, quantile, 0.025)
    q975 <- apply(y_rep, 2, quantile, 0.975)

    meta$pred     <- pp_mean
    meta$pred_lo  <- q025
    meta$pred_hi  <- q975
    meta$year     <- NA_integer_

    # Recover year from data (data has year_relative; year = year_relative + min_year)
    meta$year <- data$year_relative + min_year

    fitted_sy <- meta %>%
      group_by(state, year) %>%
      summarise(
        obs_rate  = sum(y_obs) / sum(pop) * 1e6,
        pred_rate = sum(pred)  / sum(pop) * 1e6,
        pred_lo   = sum(pred_lo) / sum(pop) * 1e6,
        pred_hi   = sum(pred_hi) / sum(pop) * 1e6,
        .groups = "drop"
      ) %>%
      mutate(
        region  = factor(REGION_MAP[state], levels = REGION_ORDER),
        state_f = factor(state, levels = STATE_ORDER)
      )

    # Aggregate regionally for clarity
    fitted_reg <- meta %>%
      mutate(region = REGION_MAP[state]) %>%
      group_by(region, year) %>%
      summarise(
        obs_rate  = sum(y_obs) / sum(pop) * 1e6,
        pred_rate = sum(pred)  / sum(pop) * 1e6,
        pred_lo   = sum(pred_lo) / sum(pop) * 1e6,
        pred_hi   = sum(pred_hi) / sum(pop) * 1e6,
        .groups = "drop"
      ) %>%
      mutate(region = factor(region, levels = REGION_ORDER))

    fit_plot <- ggplot(fitted_reg, aes(x = year)) +
      geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi, fill = region), alpha = 0.2) +
      geom_line(aes(y = pred_rate, color = region), linewidth = 0.8) +
      geom_point(aes(y = obs_rate, color = region), size = 1.5) +
      facet_wrap(~ region, scales = "free_y", ncol = 3) +
      scale_color_brewer(palette = "Set1", guide = "none") +
      scale_fill_brewer(palette = "Set1", guide = "none") +
      labs(
        title = paste0("Fitted vs Observed: ", DISEASE_LABELS[d],
                        " (", variant, ")"),
        subtitle = "Lines = posterior mean, ribbons = 95% PI, points = observed.",
        x = "Year", y = "Rate per 1,000,000"
      ) +
      theme_pub(base_size = 10)

    ggsave(file.path(out_dir, "fitted_values.png"), fit_plot,
           width = 12, height = 8, dpi = 300)
    write.csv(fitted_reg, file.path(out_dir, "fitted_values.csv"),
              row.names = FALSE)

    actual_plot <- ggplot(fitted_reg, aes(x = year, y = obs_rate, color = region)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.5) +
      facet_wrap(~ region, scales = "free_y", ncol = 3) +
      scale_color_brewer(palette = "Set1", guide = "none") +
      labs(
        title = paste0("Actual Observed Rates: ", DISEASE_LABELS[d]),
        subtitle = "Raw hospitalization rate by region.",
        x = "Year", y = "Rate per 1,000,000"
      ) +
      theme_pub(base_size = 10)

    ggsave(file.path(out_dir, "actual_values.png"), actual_plot,
           width = 12, height = 8, dpi = 300)
    write.csv(fitted_reg %>% select(region, year, obs_rate),
              file.path(out_dir, "actual_values.csv"), row.names = FALSE)

    cat("  Forecasting to", FORECAST_END, "...\n")
    n_forecast <- length(forecast_years)
    fc <- forecast_model(fit, variant, n_forecast, N_POST_DRAWS)

    state_names <- levels(data$state)

    # Build long forecast data frame
    fc_rows <- list()
    for (s in seq_along(state_names)) {
      st <- state_names[s]
      for (yr in all_years) {
        pop_sy <- pop_all %>%
          filter(state == st, year == yr) %>%
          arrange(age) %>%
          pull(population)
        if (length(pop_sy) != length(AGE_GROUPS)) next
        total_pop <- sum(pop_sy)
        if (total_pop <= 0) next

        if (yr <= max_year) {
          t_idx <- yr - min_year + 1
          a_draws <- fc$alpha_st[, s, t_idx]
        } else {
          h <- yr - max_year
          a_draws <- fc$alpha_fwd[, s, h]
        }

        total_cases <- rep(0, fc$n_use)
        for (a_idx in seq_along(AGE_GROUPS)) {
          lam <- log(pop_sy[a_idx]) + a_draws + fc$age_term[, a_idx]
          lam <- pmin(lam, 20.5)
          total_cases <- total_cases + exp(lam)
        }
        rate_draws <- total_cases / total_pop * 1e6
        fc_rows[[length(fc_rows) + 1]] <- data.frame(
          state = st, year = yr,
          mean_rate = as.numeric(mean(rate_draws)),
          lower_95 = as.numeric(quantile(rate_draws, 0.025)),
          upper_95 = as.numeric(quantile(rate_draws, 0.975)),
          lower_80 = as.numeric(quantile(rate_draws, 0.10)),
          upper_80 = as.numeric(quantile(rate_draws, 0.90)),
          is_forecast = yr > max_year
        )
      }
    }
    fc_df <- do.call(rbind, fc_rows)
    fc_df$region  <- factor(REGION_MAP[fc_df$state], levels = REGION_ORDER)
    fc_df$state_f <- factor(fc_df$state, levels = STATE_ORDER)
    write.csv(fc_df, file.path(out_dir, "forecast_data.csv"), row.names = FALSE)

    # Observed rates overlay
    obs_df <- raw %>%
      group_by(state, year) %>%
      summarise(cases = sum(.data[[d]], na.rm = TRUE),
                pop   = sum(population, na.rm = TRUE), .groups = "drop") %>%
      mutate(obs_rate = cases / pop * 1e6)
    obs_df$state_f <- factor(obs_df$state, levels = STATE_ORDER)

    # Per-state forecast charts
    fc_dir <- file.path(out_dir, "forecasts")
    dir.create(fc_dir, recursive = TRUE, showWarnings = FALSE)

    for (st in STATE_ORDER) {
      df_s <- fc_df %>% filter(state == st)
      if (nrow(df_s) == 0) next
      obs_s <- obs_df %>% filter(state == st)
      rg <- as.character(df_s$region[1])

      p <- ggplot() +
        geom_ribbon(data = df_s %>% filter(is_forecast),
                    aes(x = year, ymin = lower_95, ymax = upper_95),
                    fill = "firebrick", alpha = 0.15) +
        geom_ribbon(data = df_s %>% filter(is_forecast),
                    aes(x = year, ymin = lower_80, ymax = upper_80),
                    fill = "firebrick", alpha = 0.25) +
        geom_line(data = df_s %>% filter(!is_forecast),
                  aes(x = year, y = mean_rate),
                  color = "steelblue", linewidth = 0.8) +
        geom_line(data = df_s %>% filter(is_forecast),
                  aes(x = year, y = mean_rate),
                  color = "firebrick", linewidth = 0.8) +
        geom_point(data = obs_s, aes(x = year, y = obs_rate),
                   size = 1.5, alpha = 0.7, color = "black") +
        geom_vline(xintercept = max_year, linetype = "dashed",
                   color = "gray50", linewidth = 0.4) +
        scale_x_continuous(breaks = seq(2010, FORECAST_END, 2)) +
        labs(
          title = paste0(st, " (", rg, ") - ", DISEASE_LABELS[d],
                          " [", variant, "]"),
          subtitle = paste0("Forecast to ", FORECAST_END,
                            " with 80%/95% credible intervals."),
          x = "Year", y = "Rate per 1,000,000"
        ) +
        theme_pub()

      # Save chart data: model trajectory + observed dots
      st_csv_data <- bind_rows(
        df_s %>% mutate(layer = "model"),
        obs_s %>% transmute(state, year, obs_rate, layer = "observed")
      )
      write.csv(st_csv_data, file.path(fc_dir, paste0(st, ".csv")),
                row.names = FALSE)

      ggsave(file.path(fc_dir, paste0(st, ".png")), p,
             width = 8, height = 5, dpi = 300)
    }
    cat("  Saved 27 forecast charts.\n")

    rm(fit, y_rep, fc); gc()
  }
}

cat("\n\n=== Writing comparison tables ===\n")

overall_df <- do.call(rbind, overall_comparison)
rownames(overall_df) <- NULL
write.csv(overall_df, file.path(RESULTS_DIR, "overall_comparison.csv"),
          row.names = FALSE)

subgroup_df <- do.call(rbind, subgroup_comparison)
rownames(subgroup_df) <- NULL
write.csv(subgroup_df, file.path(RESULTS_DIR, "subgroup_comparison.csv"),
          row.names = FALSE)

# Compact cross-model summary (one row per variant x disease with key metrics)
summary_df <- overall_df %>%
  select(variant, disease, LOOIC, LOOIC_SE, MAE, MRE, RMSE,
         Cov_80, Cov_95, max_rhat, min_n_eff)
write.csv(summary_df, file.path(RESULTS_DIR, "cross_model_summary.csv"),
          row.names = FALSE)

cat("\n=== Cross-model summary ===\n")
print(summary_df, row.names = FALSE)

cat("\n=== All outputs in", RESULTS_DIR, "===\n")
cat("Files:\n")
cat("  overall_comparison.csv    - full metrics per variant x disease\n")
cat("  subgroup_comparison.csv   - metrics by age group / region\n")
cat("  cross_model_summary.csv   - compact comparison table\n")
cat("  {variant}/{disease}/      - fits, diagnostics, PPC, charts\n")
cat("  {variant}/{disease}/forecasts/  - per-state forecast charts\n")
