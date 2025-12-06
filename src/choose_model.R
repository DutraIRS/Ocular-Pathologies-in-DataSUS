source("src/data_processing.R")
source("src/model_selection.R")

library(rstan)

# Detect available cores
options(mc.cores=parallel::detectCores())

data <- prepare_data("glaucoma")

# define models
models <- list(
  "poisson_full" = list(
                                N = nrow(data),
                                N_states = length(unique(data$state_id)),
                                N_years = length(unique(data$year_id)),
                                
                                state_idx = data$state_id,
                                year_idx = data$year_id,
                                
                                age_centered = data$age_relative,
                                log_pop = data$log_pop,
                                y = data$glaucoma,
                                
                                # tight prior to avoid exp explosion
                                hyp_alpha = 5,
                                hyp_beta = 0.5
                              )
)

model_results <- list()

# run each one
for (model_name in names(models)) {
  stan_data <- models[[model_name]]
  
  output_dir <- file.path("diagnostics", model_name)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  model <- stan_model(file = paste0("models/",
                      model_name, ".stan"),
                      model_name = model_name)
  
  fit <- sampling(
    model, 
    data = stan_data, 
    iter = 2000, 
    chains = 4
  )
  
  saveRDS(fit, file.path(output_dir, "fit.rds"))
  
  metrics_vec <- run_all_diagnostics(
    fit = fit, 
    output_dir = output_dir, 
    data_y = stan_data$y, 
    pars = "sigma"
  )
  
  model_results[[model_name]] <- metrics_vec
}

comparison_table <- do.call(rbind, model_results)

write.csv(comparison_table, "diagnostics/model_comparison.csv")
print(comparison_table)