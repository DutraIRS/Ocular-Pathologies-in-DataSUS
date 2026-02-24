# Ocular Pathologies in DataSUS

A complete reproduction of results involves:

-   Set working directory to the root of this repository
-   Run `src/choose_model.R` to compare different possible models for this problem
-   Evaluate MCMC convergence and diagnostics
-   Check `diagnostics/{disease_name}/model_comparison.csv` to evaluate the best one
-   Run `src/main_analysis.R` to fit the selected model and generate results

## Model Comparison Strategy

The full spectrum of candidate models ranges from "Null" models (intercept only) to complex dynamic models (spatiotemporal interactions). In a initial step, we tested two distributional families: **Poisson (P)** and **Negative Binomial (NB)**. Having verified the overdispersion of the data, we proceed with the following negative binomial models:

| Model ID     | Intercept Structure         | Slope (Age Effect) Structure |
|:-------------|:---------------|:---------------|
| **NB0** | Global Constant ($\beta_0$) | None (Null)                  |
| **NB1** | Static State ($\alpha_s$)   | Global Constant ($\beta$)    |
| **NB2** | Static State ($\alpha_s$)   | State-Specific ($\beta_s$)   |
| **NB3** | Static State ($\alpha_s$)   | Time-Varying (RW $\beta_t$)     |
| **NB4** | Static State ($\alpha_s$)   | Interaction (RW $\beta_{s,t}$)  |
| **NB5** | Dynamic (RW $\alpha_{s,t}$) | Global Constant ($\beta$)    |
| **NB6** | Dynamic (RW $\alpha_{s,t}$) | State-Specific ($\beta_s$)   |
| **NB7** | Dynamic (RW $\alpha_{s,t}$) | Time-Varying (RW $\beta_t$)     |
| **NB8** | Dynamic (RW $\alpha_{s,t}$) | Interaction (RW $\beta_{s,t}$)  |
