library(rstan)
library(ggplot2)
library(dplyr)
library(tidyr)

source("src/data_processing.R")

BEST_MODEL     <- "nb5"
DISEASES       <- c("glaucoma", "retina", "eye")
DISEASE_LABELS <- c(glaucoma = "Glaucoma", retina = "Retinopathy",
                     eye = "Eye & Appendage Diseases")
FORECAST_END   <- 2050
N_POST_DRAWS   <- 500   # posterior draws for forecasting / PPC (speed)

AGE_GROUPS     <- c(25, 35, 45, 55, 65, 75, 90)
AGE_LABELS     <- c("25-34", "35-44", "45-54", "55-64",
                     "65-74", "75-89", "90+")
AGE_LABEL_MAP  <- setNames(AGE_LABELS, AGE_GROUPS)

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

theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 2,
                                      hjust = 0),
      plot.subtitle    = element_text(size = base_size, color = "gray40",
                                      hjust = 0),
      strip.text       = element_text(face = "bold"),
      legend.position  = "bottom",
      panel.grid.minor = element_blank()
    )
}

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

raw <- read.csv("data/sus_ocular_data.csv")
raw[raw == -1] <- NA
raw <- subset(raw, year != 2025)
raw$region    <- factor(REGION_MAP[raw$state], levels = REGION_ORDER)
raw$age_label <- factor(AGE_LABEL_MAP[as.character(raw$age)],
                        levels = AGE_LABELS)
raw$state_f   <- factor(raw$state, levels = STATE_ORDER)

min_year <- min(raw$year)
max_year <- max(raw$year)
cat("Data loaded:", nrow(raw), "obs | years", min_year, "-", max_year,
    "| states:", length(unique(raw$state)),
    "| age groups:", length(AGE_GROUPS), "\n")


cat("\n=== Figure 1: Temporal Trends ===\n")

for (d in DISEASES) {
  d_label <- DISEASE_LABELS[d]

  # Aggregate by region (not state) so there's one line per region
  regional <- raw %>%
    filter(!is.na(.data[[d]])) %>%
    group_by(region, year) %>%
    summarise(cases = sum(.data[[d]], na.rm = TRUE),
              pop   = sum(population, na.rm = TRUE), .groups = "drop") %>%
    mutate(rate = cases / pop * 1e6)

  national <- raw %>%
    filter(!is.na(.data[[d]])) %>%
    group_by(year) %>%
    summarise(cases = sum(.data[[d]], na.rm = TRUE),
              pop   = sum(population, na.rm = TRUE), .groups = "drop") %>%
    mutate(rate = cases / pop * 1e6)

  fig1 <- ggplot() +
    geom_line(data = regional,
              aes(x = year, y = rate, color = region),
              linewidth = 0.8) +
    geom_point(data = regional,
               aes(x = year, y = rate, color = region),
               size = 1.2) +
    geom_line(data = national, aes(x = year, y = rate),
              color = "black", linewidth = 1.3) +
    geom_point(data = national, aes(x = year, y = rate),
               color = "black", size = 1.5) +
    annotate("text",
             x = max_year - 0.3,
             y = max(national$rate, na.rm = TRUE) * 1.08,
             label = "National", fontface = "bold", size = 3.5, hjust = 1) +
    scale_x_continuous(breaks = seq(2010, 2024, 2)) +
    scale_color_brewer(palette = "Set1", name = "Region") +
    labs(
      title    = paste0("Temporal Trends of ", d_label, " Prevalence by Region"),
      subtitle = paste0("Hospitalization rate per 1,000,000 inhabitants (",
                        min_year, "-", max_year,
                        "). Black line = national average."),
      x = "Year", y = "Rate per 1,000,000"
    ) +
    theme_pub() +
    guides(color = guide_legend(nrow = 1))

  ggsave(paste0("results/fig1_temporal_trends_", d, ".png"), fig1,
         width = 12, height = 7, dpi = 300)
  ggsave(paste0("results/fig1_temporal_trends_", d, ".pdf"), fig1,
         width = 12, height = 7)
  cat("  Saved:", d, "\n")
}


cat("\n=== Figure 2: Age-Stratified ===\n")

for (d in DISEASES) {
  d_label <- DISEASE_LABELS[d]

  # Aggregate by region (not state)
  regional_age <- raw %>%
    filter(!is.na(.data[[d]])) %>%
    group_by(region, year, age_label) %>%
    summarise(cases = sum(.data[[d]], na.rm = TRUE),
              pop   = sum(population, na.rm = TRUE), .groups = "drop") %>%
    mutate(rate = cases / pop * 1e6)

  national_age <- raw %>%
    filter(!is.na(.data[[d]])) %>%
    group_by(year, age_label) %>%
    summarise(cases = sum(.data[[d]], na.rm = TRUE),
              pop   = sum(population, na.rm = TRUE), .groups = "drop") %>%
    mutate(rate = cases / pop * 1e6)

  fig2 <- ggplot() +
    geom_line(data = regional_age,
              aes(x = year, y = rate, color = region),
              linewidth = 0.6) +
    geom_line(data = national_age, aes(x = year, y = rate),
              color = "black", linewidth = 1.0) +
    facet_wrap(~ age_label, scales = "free_y", ncol = 4) +
    scale_x_continuous(breaks = seq(2010, 2024, 4)) +
    scale_color_brewer(palette = "Set1", name = "Region") +
    labs(
      title    = paste0("Age-Stratified Temporal Dynamics: ", d_label),
      subtitle = "Rate per 1,000,000 by region and age group. Black line = national average.",
      x = "Year", y = "Rate per 1,000,000"
    ) +
    theme_pub(base_size = 10) +
    guides(color = guide_legend(nrow = 1))

  ggsave(paste0("results/fig2_age_stratified_", d, ".png"), fig2,
         width = 16, height = 10, dpi = 300)
  ggsave(paste0("results/fig2_age_stratified_", d, ".pdf"), fig2,
         width = 16, height = 10)
  cat("  Saved:", d, "\n")
}


cat("\n=== Table 1a: Model Comparison ===\n")

all_comp <- list()
for (d in DISEASES) {
  path <- file.path("diagnostics", d, "model_comparison.csv")
  if (!file.exists(path)) { cat("  Missing:", path, "\n"); next }
  comp <- read.csv(path, row.names = 1)
  comp$disease <- DISEASE_LABELS[d]
  comp$model   <- toupper(rownames(comp))
  rownames(comp) <- NULL
  all_comp[[d]] <- comp
}
comparison_df <- do.call(rbind, all_comp)

table1a <- comparison_df %>%
  transmute(
    Disease     = disease,
    Model       = model,
    LOOIC       = sprintf("%.1f", looic),
    LOOIC_SE    = sprintf("%.1f", looic_se),
    WAIC        = sprintf("%.1f", waic),
    Log_ML      = sprintf("%.1f", log_ml),
    Max_Rhat    = sprintf("%.4f", max_rhat),
    Min_ESS     = sprintf("%.0f", min_n_eff),
    Converged   = ifelse(max_rhat < 1.01 & min_n_eff > 400, "Yes", "No")
  )

write.csv(table1a, "results/table1a_model_comparison.csv", row.names = FALSE)
cat("  Saved.\n")


cat("\n=== Table 1b: NB5 Demographic Performance ===\n")

subgroup_results <- list()

for (d in DISEASES) {
  fit_path <- file.path("diagnostics", d, BEST_MODEL, "fit.rds")
  if (!file.exists(fit_path)) {
    cat("  Skipping", d, "- fit.rds not found.\n")
    next
  }

  cat("  Loading", d, "...")
  fit   <- readRDS(fit_path)
  y_rep <- as.matrix(fit, pars = "y_rep")
  rm(fit); gc()

  data_prep <- prepare_data(d)
  N <- nrow(data_prep)
  stopifnot(ncol(y_rep) == N)

  # Reconstruct observation metadata
  meta <- data.frame(
    y_obs = data_prep$disease,
    state = as.character(data_prep$state),
    age   = round(data_prep$age_relative * 10 + 60),
    pop   = exp(data_prep$log_pop),
    stringsAsFactors = FALSE
  )
  meta$age_label <- factor(AGE_LABEL_MAP[as.character(meta$age)],
                           levels = AGE_LABELS)
  meta$region <- factor(REGION_MAP[meta$state], levels = REGION_ORDER)

  # Subsample draws for speed
  n_use <- min(nrow(y_rep), N_POST_DRAWS)
  idx   <- sample(nrow(y_rep), n_use)
  y_rep <- y_rep[idx, ]

  # Helper: compute performance metrics for a subset of observations
  eval_group <- function(obs_idx, grp_type, grp_val) {
    y   <- meta$y_obs[obs_idx]
    pop <- meta$pop[obs_idx]
    yr  <- y_rep[, obs_idx, drop = FALSE]

    pp_mean <- colMeans(yr)
    obs_r   <- y / pop * 1e6
    pred_r  <- pp_mean / pop * 1e6

    q025 <- apply(yr, 2, quantile, 0.025)
    q975 <- apply(yr, 2, quantile, 0.975)
    q10  <- apply(yr, 2, quantile, 0.10)
    q90  <- apply(yr, 2, quantile, 0.90)

    data.frame(
      Disease = DISEASE_LABELS[d],
      Group   = grp_type,
      Subgroup = grp_val,
      N       = length(obs_idx),
      MAE     = round(mean(abs(obs_r - pred_r)), 2),
      RMSE    = round(sqrt(mean((obs_r - pred_r)^2)), 2),
      Bias    = round(mean(pred_r - obs_r), 2),
      Cov_80  = round(mean(y >= q10 & y <= q90) * 100, 1),
      Cov_95  = round(mean(y >= q025 & y <= q975) * 100, 1),
      stringsAsFactors = FALSE
    )
  }

  # By age group
  for (ag in AGE_LABELS) {
    ii <- which(meta$age_label == ag)
    if (length(ii) > 0)
      subgroup_results[[length(subgroup_results) + 1]] <-
        eval_group(ii, "Age", ag)
  }
  # By region
  for (rg in REGION_ORDER) {
    ii <- which(meta$region == rg)
    if (length(ii) > 0)
      subgroup_results[[length(subgroup_results) + 1]] <-
        eval_group(ii, "Region", rg)
  }
  # Overall
  subgroup_results[[length(subgroup_results) + 1]] <-
    eval_group(seq_len(N), "Overall", "All")

  rm(y_rep); gc()
  cat(" done.\n")
}

table1b <- do.call(rbind, subgroup_results)
write.csv(table1b, "results/table1b_nb5_demographic_performance.csv",
          row.names = FALSE)
cat("  Saved.\n")


cat("\n=== Figure 3: Forecasting ===\n")

forecast_years <- (max_year + 1):FORECAST_END
n_forecast     <- length(forecast_years)
all_years      <- min_year:FORECAST_END
age_cent_vals  <- (AGE_GROUPS - 60) / 10

pop_obs <- raw %>%
  select(state, year, age, population) %>%
  distinct()

IBGE_PATH <- "data/ibge_population_projections.csv"

cat("  Loading IBGE official population projections...\n")
pop_ibge <- read.csv(IBGE_PATH)
names(pop_ibge) <- c("state", "year", "age", "population")
pop_future <- pop_ibge %>%
  filter(year %in% forecast_years)
cat("  IBGE projections loaded:", nrow(pop_future), "rows.\n")

  
pop_all <- rbind(pop_obs, pop_future)
cat("  Population ready:", nrow(pop_all), "total rows.\n")

forecast_all <- list()

for (d in DISEASES) {
  fit_path <- file.path("diagnostics", d, BEST_MODEL, "fit.rds")
  if (!file.exists(fit_path)) {
    cat("  Skipping", d, "- fit not found.\n")
    next
  }

  cat("  Forecasting", d, "...")
  fit  <- readRDS(fit_path)
  post <- rstan::extract(fit)
  rm(fit); gc()

  alpha_full <- post$alpha_st        # [draws, states, years]
  beta_full  <- post$beta_age        # [draws] (global age effect in NB5)
  sigma_full <- post$sigma_alpha_rw  # [draws]
  rm(post); gc()

  n_draws  <- dim(alpha_full)[1]
  n_states <- dim(alpha_full)[2]
  n_years  <- dim(alpha_full)[3]

  # State name ordering (matches Stan indices)
  dp <- prepare_data(d)
  state_names <- levels(dp$state)
  rm(dp)

  # Subsample posterior draws
  n_use <- min(n_draws, N_POST_DRAWS)
  di    <- sample(n_draws, n_use)

  alpha_st <- alpha_full[di, , , drop = FALSE]  # [n_use, states, years]
  beta_age <- as.numeric(beta_full[di])          # [n_use] plain vector
  sigma_rw <- as.numeric(sigma_full[di])         # [n_use] plain vector
  rm(alpha_full, beta_full, sigma_full); gc()

  # Forward-simulate random walk for alpha (vectorized over draws & states)
  alpha_last <- alpha_st[, , n_years]  # [n_use, n_states]
  alpha_fwd  <- array(NA, dim = c(n_use, n_states, n_forecast))

  for (h in seq_len(n_forecast)) {
    eps <- sweep(matrix(rnorm(n_use * n_states), n_use, n_states),
                 1, sigma_rw, "*")
    alpha_last <- alpha_last + eps
    alpha_fwd[, , h] <- alpha_last
  }

  # Compute predicted rates for all state-years
  results <- vector("list", n_states * length(all_years))
  k <- 0

  for (s in seq_len(n_states)) {
    st <- state_names[s]

    for (yr in all_years) {
      pop_sy <- pop_all %>%
        filter(state == st, year == yr) %>%
        arrange(age) %>%
        pull(population)

      if (length(pop_sy) != length(AGE_GROUPS)) next
      total_pop <- sum(pop_sy)
      if (total_pop <= 0) next

      # Get alpha draws for this state-year
      if (yr <= max_year) {
        t_idx    <- yr - min_year + 1
        a_draws  <- alpha_st[, s, t_idx]
      } else {
        h        <- yr - max_year
        a_draws  <- alpha_fwd[, s, h]
      }

      # Sum predicted cases across age groups (vectorized over draws)
      total_cases <- rep(0, n_use)
      for (a_idx in seq_along(AGE_GROUPS)) {
        log_lam <- log(pop_sy[a_idx]) + a_draws +
                   beta_age * age_cent_vals[a_idx]
        log_lam <- pmin(log_lam, 20.5)
        total_cases <- total_cases + exp(log_lam)
      }

      rate_draws <- total_cases / total_pop * 1e6

      k <- k + 1
      results[[k]] <- data.frame(
        state       = st,
        year        = yr,
        disease     = d,
        median_rate = as.numeric(median(rate_draws)),
        lower_95    = as.numeric(quantile(rate_draws, 0.025)),
        upper_95    = as.numeric(quantile(rate_draws, 0.975)),
        lower_80    = as.numeric(quantile(rate_draws, 0.10)),
        upper_80    = as.numeric(quantile(rate_draws, 0.90)),
        is_forecast = yr > max_year,
        stringsAsFactors = FALSE
      )
    }
  }

  forecast_all[[d]] <- do.call(rbind, results[seq_len(k)])
  rm(alpha_st, alpha_fwd, beta_age, sigma_rw); gc()
  cat(" done.\n")
}

forecast_df <- do.call(rbind, forecast_all)
forecast_df$region        <- factor(REGION_MAP[forecast_df$state],
                                    levels = REGION_ORDER)
forecast_df$disease_label <- DISEASE_LABELS[forecast_df$disease]
forecast_df$state_f       <- factor(forecast_df$state, levels = STATE_ORDER)

write.csv(forecast_df, "results/forecast_data.csv", row.names = FALSE)

obs_rates <- list()
for (d in DISEASES) {
  obs_d <- raw %>%
    group_by(state, year) %>%
    summarise(cases = sum(.data[[d]], na.rm = TRUE),
              pop   = sum(population, na.rm = TRUE), .groups = "drop") %>%
    mutate(obs_rate = cases / pop * 1e6, disease = d)
  obs_rates[[d]] <- obs_d
}
obs_df <- do.call(rbind, obs_rates)
obs_df$state_f <- factor(obs_df$state, levels = STATE_ORDER)

if (!dir.exists("results/forecasts")) dir.create("results/forecasts", recursive = TRUE)

for (d in DISEASES) {
  d_label <- DISEASE_LABELS[d]
  df_fc_d <- forecast_df %>% filter(disease == d)
  df_obs_d <- obs_df %>% filter(disease == d)

  for (st in STATE_ORDER) {
    df_fc_s  <- df_fc_d %>% filter(state == st)
    df_obs_s <- df_obs_d %>% filter(state == st)
    if (nrow(df_fc_s) == 0) next

    st_region <- as.character(df_fc_s$region[1])

    fig3 <- ggplot() +
      geom_ribbon(
        data = df_fc_s %>% filter(is_forecast),
        aes(x = year, ymin = lower_95, ymax = upper_95),
        fill = "firebrick", alpha = 0.15
      ) +
      geom_ribbon(
        data = df_fc_s %>% filter(is_forecast),
        aes(x = year, ymin = lower_80, ymax = upper_80),
        fill = "firebrick", alpha = 0.25
      ) +
      geom_line(
        data = df_fc_s %>% filter(!is_forecast),
        aes(x = year, y = median_rate),
        color = "steelblue", linewidth = 0.8
      ) +
      geom_line(
        data = df_fc_s %>% filter(is_forecast),
        aes(x = year, y = median_rate),
        color = "firebrick", linewidth = 0.8
      ) +
      geom_point(
        data = df_obs_s,
        aes(x = year, y = obs_rate),
        size = 1.5, alpha = 0.7, color = "black"
      ) +
      geom_vline(xintercept = max_year, linetype = "dashed",
                 color = "gray50", linewidth = 0.4) +
      scale_x_continuous(breaks = seq(2010, 2050, 5)) +
      labs(
        title    = paste0(st, " (", st_region, ") - ", d_label, " Forecast"),
        subtitle = "NB5 posterior projection with 80%/95% credible intervals.",
        x = "Year",
        y = "Rate per 1,000,000"
      ) +
      theme_pub()

    ggsave(paste0("results/forecasts/fig3_", d, "_", st, ".png"), fig3,
           width = 8, height = 5, dpi = 300)
  }
  cat("  Saved 27 charts for", d, "\n")
}

cat("\n=== All outputs saved to results/ ===\n")
cat("Files:\n")
cat("  results/fig1_temporal_trends_{glaucoma,retina,eye}.{png,pdf}\n")
cat("  results/fig2_age_stratified_{glaucoma,retina,eye}.{png,pdf}\n")
cat("  results/table1a_model_comparison.csv\n")
cat("  results/table1b_nb5_demographic_performance.csv\n")
cat("  results/forecasts/fig3_{disease}_{state}.png  (27 x 3 = 81 charts)\n")
cat("  results/forecast_data.csv\n")
