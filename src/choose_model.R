source("src/data_processing.R")
source("src/model_selection.R")

library(rstan)

# Detect available cores
options(mc.cores=parallel::detectCores())

disease_names <- c("glaucoma", "retina", "eye")

for (disease_name in disease_names) {
  message("Preparing data for disease: ", disease_name)
  
  data <- prepare_data(disease_name)
  
  stan_data <- list(
    N = nrow(data),
    N_states = length(unique(data$state_id)),
    
    N_years = length(unique(data$year_id)),
    
    state_idx = data$state_id,
    year_idx = data$year_id,
    
    age_centered = data$age_relative,
    log_pop = data$log_pop,
    y = data$disease
  )
  
  # list models to evaluate
  models <- c("nb0", "nb1" , "nb2", "nb3", "nb4", "nb5", "nb6", "nb7", "nb8")
  
  model_results <- list()
  
  # run each model
  for (model_name in models) {
    message("Fitting model: ", model_name)
    
    # create folder
    output_dir <- file.path("diagnostics", disease_name, model_name)
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    # stan magic
    model <- stan_model(file = paste0("models/",
                                      model_name, ".stan"),
                        model_name = model_name)
    
    fit <- sampling(
      model, 
      data = stan_data, 
      iter = 4000, 
      chains = 4,
      init = 0,
      control = list(
        adapt_delta = 0.8,  # default: 0.8, if bigger, smaller steps
        max_treedepth = 10   # default: 10, if bigger, deeper exploration
      )
    )
    
    # save for ad hoc
    saveRDS(fit, file.path(output_dir, "fit.rds"))
    
    message("Evaluating diagnostics...")
    # eval
    metrics_vec <- run_all_diagnostics(
      fit = fit, 
      output_dir = output_dir, 
      data_y = stan_data$y, 
      pars = "beta"
    )
    # bookkepping
    model_results[[model_name]] <- metrics_vec
    
    rm(fit)
    gc()
  }
  
  comparison_table <- do.call(rbind, model_results)
  
  write.csv(comparison_table,
            paste0("diagnostics/", disease_name, "/model_comparison.csv"),
            row.names = TRUE)
  print(comparison_table)
}
