library(rstan)
library(loo)
library(ggplot2)
library(dplyr)
library(tidyr)

source("src/data_processing.R")
source("src/model_selection.R")

options(mc.cores = parallel::detectCores() - 2)
rstan_options(auto_write = TRUE)

VALIDATION_VARIANTS <- c("nb5_age_cat")   # selected model
DISEASES       <- c("glaucoma", "retina", "eye")
DISEASE_LABELS <- c(glaucoma = "Glaucoma", retina = "Retinopathy",
                     eye = "Eye & Appendage Diseases")
AGE_GROUPS     <- c(25, 35, 45, 55, 65, 75, 90)
HOLDOUT_YEARS  <- 2023:2024
TRAIN_END_YEAR <- min(HOLDOUT_YEARS) - 1                # 2022
N_POST_DRAWS   <- 500
FORCE_REFIT_HOLDOUT <- FALSE

MODELS_DIR  <- "ad_hoc/models"
HOLDOUT_DIR <- "ad_hoc/results_holdout"
dir.create(HOLDOUT_DIR, recursive = TRUE, showWarnings = FALSE)

is_llt       <- function(v) grepl("llt", v)
is_damped    <- function(v) grepl("damped", v)
has_drift    <- function(v) grepl("drift", v) && !grepl("damped", v) && !grepl("ar1_trend", v)
is_ar1_plain <- function(v) v %in% c("nb5_ar1", "nb5_ar1_hier",
                                       "nb5_ar1_age_quad", "nb5_ar1_age_cat")
is_ar1_trend <- function(v) v == "nb5_ar1_trend"
is_ar2       <- function(v) v == "nb5_ar2"
is_ar_family <- function(v) is_ar1_plain(v) || is_ar1_trend(v) || is_ar2(v)
is_quad_age  <- function(v) grepl("age_quad", v)
is_cat_age   <- function(v) grepl("age_cat", v)

prepare_data_train <- function(disease, train_end = TRAIN_END_YEAR) {
  data <- read.csv("data/sus_ocular_data.csv")
  data[data == -1] <- NA
  data$state <- as.factor(data$state)
  data$age_relative <- (data$age - 60) / 10
  data <- subset(data, year <= train_end)
  data$year_relative <- data$year - min(data$year)
  data$log_pop <- log(data$population)
  data$disease <- data[[disease]]
  data <- data[c("state", "year_relative", "age_relative", "log_pop", "disease")]
  data <- na.omit(data)
  data$state_id <- as.integer(as.factor(data$state))
  data$year_id  <- data$year_relative + 1
  return(data)
}

build_stan_data <- function(data, variant) {
  age_levels <- c(-3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 3.0)
  age_idx    <- match(round(data$age_relative, 1), round(age_levels, 1))
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
  if (is_cat_age(variant)) {
    ds$age_idx      <- age_idx
    ds$N_age_groups <- length(AGE_GROUPS)
  }
  ds
}

forecast_holdout <- function(fit, variant, n_holdout) {
  post <- rstan::extract(fit)
  alpha_full <- post$alpha_st
  n_draws  <- dim(alpha_full)[1]
  n_states <- dim(alpha_full)[2]
  n_years  <- dim(alpha_full)[3]

  di <- sample(n_draws, min(n_draws, N_POST_DRAWS))
  n_use <- length(di)

  alpha_st <- alpha_full[di, , , drop = FALSE]

  # Age term
  age_x <- (AGE_GROUPS - 60) / 10
  if (is_cat_age(variant)) {
    age_term <- post$beta_age_cat[di, , drop = FALSE]
  } else if (is_quad_age(variant)) {
    age_term <- outer(as.numeric(post$beta_age[di]),  age_x) +
                outer(as.numeric(post$beta_age2[di]), age_x^2)
  } else {
    age_term <- outer(as.numeric(post$beta_age[di]), age_x)
  }

  sigma_alpha <- if (is_ar_family(variant))
    as.numeric(post$sigma_alpha[di])
  else
    as.numeric(post$sigma_alpha_rw[di])

  alpha_last <- alpha_st[, , n_years]
  alpha_fwd  <- array(NA, dim = c(n_use, n_states, n_holdout))

  if (is_ar1_plain(variant)) {
    rho_draws <- as.numeric(post$rho[di])
    mu_s_mat  <- post$mu_s[di, , drop = FALSE]
    for (h in seq_len(n_holdout)) {
      eps <- sweep(matrix(rnorm(n_use * n_states), n_use, n_states),
                   1, sigma_alpha, "*")
      centered <- alpha_last - mu_s_mat
      centered <- sweep(centered, 1, rho_draws, "*")
      alpha_last <- mu_s_mat + centered + eps
      alpha_fwd[, , h] <- alpha_last
    }
  } else if (is_ar1_trend(variant)) {
    rho_draws   <- as.numeric(post$rho[di])
    delta_draws <- as.numeric(post$delta[di])
    mu_s_mat    <- post$mu_s[di, , drop = FALSE]
    for (h in seq_len(n_holdout)) {
      t_now <- n_years + h
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
  } else if (has_drift(variant)) {
    mu_draws <- as.numeric(post$mu[di])
    for (h in seq_len(n_holdout)) {
      eps <- sweep(matrix(rnorm(n_use * n_states), n_use, n_states),
                   1, sigma_alpha, "*")
      alpha_last <- sweep(alpha_last + eps, 1, mu_draws, "+")
      alpha_fwd[, , h] <- alpha_last
    }
  } else {
    # Plain RW (and quadratic-age RW)
    for (h in seq_len(n_holdout)) {
      eps <- sweep(matrix(rnorm(n_use * n_states), n_use, n_states),
                   1, sigma_alpha, "*")
      alpha_last <- alpha_last + eps
      alpha_fwd[, , h] <- alpha_last
    }
  }

  list(alpha_fwd = alpha_fwd, age_term = age_term, n_use = n_use)
}

raw <- read.csv("data/sus_ocular_data.csv")
raw[raw == -1] <- NA
raw <- subset(raw, year != 2025)

holdout_summary <- list()

for (variant in VALIDATION_VARIANTS) {
  cat("\n##############################\n# VARIANT:", variant, "\n##############################\n")

  for (d in DISEASES) {
    cat("\n--- Disease:", d, "---\n")

    out_dir  <- file.path(HOLDOUT_DIR, variant, d)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    fit_path <- file.path(out_dir, "fit_train.rds")

    train_data <- prepare_data_train(d, TRAIN_END_YEAR)
    stan_train <- build_stan_data(train_data, variant)

    if (!FORCE_REFIT_HOLDOUT && file.exists(fit_path)) {
      cat("  Loading existing training fit...\n")
      fit <- readRDS(fit_path)
    } else {
      cat("  Fitting on 2010-", TRAIN_END_YEAR, "...\n", sep = "")
      model <- stan_model(file = file.path(MODELS_DIR, paste0(variant, ".stan")),
                          model_name = paste0(variant, "_train"))
      ad <- if (is_llt(variant) || is_damped(variant) || is_ar_family(variant)) 0.95 else 0.8
      fit <- sampling(model, data = stan_train,
                      iter = 4000, chains = 4, init = 0,
                      control = list(adapt_delta = ad, max_treedepth = 12))
      saveRDS(fit, fit_path)
      rm(model)
    }

    # Forecast hold-out years
    n_holdout <- length(HOLDOUT_YEARS)
    fc <- forecast_holdout(fit, variant, n_holdout)

    # Held-out denominators (observed)
    pop_obs_h <- raw %>%
      select(state, year, age, population) %>%
      distinct() %>%
      filter(year %in% HOLDOUT_YEARS)

    state_names <- levels(as.factor(raw$state))
    rows <- list()
    for (s in seq_along(state_names)) {
      st <- state_names[s]
      for (h in seq_len(n_holdout)) {
        yr <- TRAIN_END_YEAR + h
        pop_sy <- pop_obs_h %>%
          filter(state == st, year == yr) %>%
          arrange(age) %>%
          pull(population)
        if (length(pop_sy) != length(AGE_GROUPS)) next
        total_pop <- sum(pop_sy)
        if (total_pop <= 0) next
        a_draws <- fc$alpha_fwd[, s, h]
        total_cases <- rep(0, fc$n_use)
        for (a_idx in seq_along(AGE_GROUPS)) {
          lam <- log(pop_sy[a_idx]) + a_draws + fc$age_term[, a_idx]
          lam <- pmin(lam, 20.5)
          total_cases <- total_cases + exp(lam)
        }
        rate_draws <- total_cases / total_pop * 1e6
        rows[[length(rows) + 1]] <- data.frame(
          state = st, year = yr,
          mean_rate = mean(rate_draws),
          lo_80 = as.numeric(quantile(rate_draws, 0.10)),
          hi_80 = as.numeric(quantile(rate_draws, 0.90)),
          lo_95 = as.numeric(quantile(rate_draws, 0.025)),
          hi_95 = as.numeric(quantile(rate_draws, 0.975))
        )
      }
    }
    proj <- do.call(rbind, rows)

    # Observed hold-out rates
    obs_h <- raw %>%
      filter(year %in% HOLDOUT_YEARS, !is.na(.data[[d]])) %>%
      group_by(state, year) %>%
      summarise(cases = sum(.data[[d]], na.rm = TRUE),
                pop   = sum(population, na.rm = TRUE), .groups = "drop") %>%
      mutate(obs_rate = cases / pop * 1e6) %>%
      select(state, year, obs_rate)

    eval_df <- proj %>%
      inner_join(obs_h, by = c("state", "year")) %>%
      mutate(
        in_80   = obs_rate >= lo_80 & obs_rate <= hi_80,
        in_95   = obs_rate >= lo_95 & obs_rate <= hi_95,
        abs_err = abs(obs_rate - mean_rate),
        rel_err = abs(obs_rate - mean_rate) / pmax(obs_rate, 1)
      )

    write.csv(eval_df, file.path(out_dir, "holdout_predictions.csv"),
              row.names = FALSE)

    holdout_summary[[length(holdout_summary) + 1]] <- data.frame(
      variant      = variant,
      disease      = DISEASE_LABELS[d],
      MAE          = round(mean(eval_df$abs_err), 3),
      MRE          = round(mean(eval_df$rel_err), 4),
      RMSE         = round(sqrt(mean((eval_df$obs_rate - eval_df$mean_rate)^2)), 3),
      Bias         = round(mean(eval_df$mean_rate - eval_df$obs_rate), 3),
      Cov_80       = round(mean(eval_df$in_80) * 100, 1),
      Cov_95       = round(mean(eval_df$in_95) * 100, 1),
      N_obs        = nrow(eval_df),
      stringsAsFactors = FALSE
    )

    rm(fit); gc()
  }
}

holdout_df <- do.call(rbind, holdout_summary)
write.csv(holdout_df, file.path(HOLDOUT_DIR, "holdout_summary.csv"),
          row.names = FALSE)
          
naive_rows <- list()
for (d in DISEASES) {
  long <- raw %>%
    filter(!is.na(.data[[d]])) %>%
    group_by(state, year) %>%
    summarise(cases = sum(.data[[d]], na.rm = TRUE),
              pop   = sum(population, na.rm = TRUE), .groups = "drop") %>%
    mutate(rate = cases / pop * 1e6)

  lvcf <- long %>% filter(year == 2022) %>%
    rename(pred = rate) %>% select(state, pred)

  obs_h <- long %>% filter(year %in% HOLDOUT_YEARS) %>%
    inner_join(lvcf, by = "state")

  lin_pred <- long %>%
    filter(year >= 2018, year <= 2022) %>%
    group_by(state) %>%
    do({
      df <- .
      if (nrow(df) < 2) return(data.frame(year = HOLDOUT_YEARS, pred_lin = NA))
      m <- lm(rate ~ year, data = df)
      data.frame(year = HOLDOUT_YEARS,
                 pred_lin = as.numeric(predict(m, newdata =
                                                  data.frame(year = HOLDOUT_YEARS))))
    }) %>% ungroup()

  obs_h <- obs_h %>% inner_join(lin_pred, by = c("state", "year"))

  naive_rows[[length(naive_rows) + 1]] <- data.frame(
    variant = "naive_LVCF", disease = DISEASE_LABELS[d],
    MAE  = round(mean(abs(obs_h$rate - obs_h$pred), na.rm = TRUE), 3),
    MRE  = round(mean(abs(obs_h$rate - obs_h$pred) / pmax(obs_h$rate, 1), na.rm = TRUE), 4),
    RMSE = round(sqrt(mean((obs_h$rate - obs_h$pred)^2, na.rm = TRUE)), 3),
    Bias = round(mean(obs_h$pred - obs_h$rate, na.rm = TRUE), 3),
    Cov_80 = NA, Cov_95 = NA, N_obs = sum(!is.na(obs_h$pred)),
    stringsAsFactors = FALSE
  )
  naive_rows[[length(naive_rows) + 1]] <- data.frame(
    variant = "naive_linear", disease = DISEASE_LABELS[d],
    MAE  = round(mean(abs(obs_h$rate - obs_h$pred_lin), na.rm = TRUE), 3),
    MRE  = round(mean(abs(obs_h$rate - obs_h$pred_lin) / pmax(obs_h$rate, 1), na.rm = TRUE), 4),
    RMSE = round(sqrt(mean((obs_h$rate - obs_h$pred_lin)^2, na.rm = TRUE)), 3),
    Bias = round(mean(obs_h$pred_lin - obs_h$rate, na.rm = TRUE), 3),
    Cov_80 = NA, Cov_95 = NA, N_obs = sum(!is.na(obs_h$pred_lin)),
    stringsAsFactors = FALSE
  )
}
naive_df <- do.call(rbind, naive_rows)

full_df <- rbind(holdout_df, naive_df)
write.csv(full_df, file.path(HOLDOUT_DIR, "holdout_with_naive.csv"),
          row.names = FALSE)

cat("\n=== Out-of-sample (2023-2024) holdout, with naive baselines ===\n")
print(full_df, row.names = FALSE)

cat("\nOutputs:\n")
cat("  ", file.path(HOLDOUT_DIR, "holdout_summary.csv"), "\n")
cat("  ", file.path(HOLDOUT_DIR, "holdout_with_naive.csv"), "\n")
cat("  ", file.path(HOLDOUT_DIR, "{variant}/{disease}/holdout_predictions.csv"), "\n")
