
library(dplyr)
library(ggplot2)
library(tidyr)
library(rstan)

source("src/data_processing.R")

BEST_MODEL <- "nb5_age_cat"
DISEASES <- c("glaucoma", "retina", "eye")
DISEASE_LABELS <- c(glaucoma = "Glaucoma", retina = "Retinopathy", eye = "Eye & Appendage Diseases")

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

all_alpha <- list()

for (d in DISEASES) {
  cat(paste("Processing", d, "...\n"))
  fit_file <- file.path("ad_hoc", "results", BEST_MODEL, d, "fit.rds")

  if (!file.exists(fit_file)) {
      cat("File not found:", fit_file, "\n")
      next
  }
  fit <- readRDS(fit_file)
  post <- rstan::extract(fit)
  
  alpha_full <- post$alpha_st
  
  dp <- prepare_data(d)
  state_names <- levels(factor(dp$state))
  
  years <- 2010:(2010 + dim(alpha_full)[3] - 1)
  
  for(s_idx in 1:length(state_names)) {
    s_name <- state_names[s_idx]
    region <- REGION_MAP[s_name]
    
    for(t_idx in 1:length(years)) {
      y <- years[t_idx]
      vals <- alpha_full[, s_idx, t_idx]
      
      all_alpha[[length(all_alpha) + 1]] <- data.frame(
        disease = DISEASE_LABELS[d],
        state = s_name,
        region = region,
        year = y,
        mean = mean(vals),
        lo_95 = quantile(vals, 0.025),
        hi_95 = quantile(vals, 0.975)
      )
    }
  }
  rm(fit, post, alpha_full)
  gc()
}

if(length(all_alpha) > 0) {
    df_alpha <- bind_rows(all_alpha)
    write.csv(df_alpha, "results/state_effects_evolution.csv", row.names = FALSE)
    
    # Load population data to calculate weighted averages
    sus_data <- read.csv("data/sus_ocular_data.csv")
    pop_data <- sus_data %>%
      group_by(state, year) %>%
      summarise(pop = sum(population, na.rm = TRUE), .groups = "drop")
      
    # Join pop to df_alpha
    df_alpha_pop <- df_alpha %>%
      left_join(pop_data, by = c("state", "year")) %>%
      mutate(pop = ifelse(is.na(pop), 1, pop))
      
    # Calculate regional aggregates (population-weighted mean)
    df_regional <- df_alpha_pop %>%
      group_by(disease, region, year) %>%
      summarise(
        mean = weighted.mean(mean, pop, na.rm = TRUE), 
        lo_95 = weighted.mean(lo_95, pop, na.rm = TRUE), 
        hi_95 = weighted.mean(hi_95, pop, na.rm = TRUE), 
        .groups="drop"
      ) %>%
      mutate(state = paste("Region Average:", region))
      
    # Calculate national aggregate (population-weighted mean)
    df_national <- df_alpha_pop %>%
      group_by(disease, year) %>%
      summarise(
        mean = weighted.mean(mean, pop, na.rm = TRUE), 
        lo_95 = weighted.mean(lo_95, pop, na.rm = TRUE), 
        hi_95 = weighted.mean(hi_95, pop, na.rm = TRUE), 
        .groups="drop"
      ) %>%
      mutate(state = "National Average", region = "National")

    df_all_plots <- bind_rows(df_alpha, df_regional, df_national)

    out_dir <- "results/state_effects"
    dir.create(out_dir, showWarnings = FALSE)
    
    unique_entities <- unique(df_all_plots$state)
    for(s in unique_entities) {
        df_state <- df_all_plots %>% filter(state == s)
        
        p <- ggplot(df_state, aes(x = year, y = mean, color = disease, fill = disease)) +
            geom_line(linewidth = 1) +
            geom_ribbon(aes(ymin = lo_95, ymax = hi_95), alpha = 0.2, color = NA) +
            theme_minimal() +
            labs(
                title = paste("Evolution of State Effects -", s),
                subtitle = paste("Region:", df_state$region[1]),
                x = "Year",
                y = expression("Posterior Mean of "*alpha),
                color = "Disease",
                fill = "Disease"
            )

        safe_name <- gsub("[: ]", "_", s)
        ggsave(file.path(out_dir, paste0("state_effects_", safe_name, ".png")), plot = p, width = 8, height = 6)
    }
    cat("State and aggregate plots saved to", out_dir, "\n")
} else {
    cat("No data parsed.\n")
}

