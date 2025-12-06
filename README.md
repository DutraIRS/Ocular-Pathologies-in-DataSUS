# Ocular Pathologies in DataSUS

A complete reproduction of results involves:

-   Set working directory to the root of this repository
-   Run `src/choose_model.R` to compare different possible models for this problem
-   Evaluate MCMC convergence and diagnostics
-   Check `diagnostics/model_comparison.csv` to evaluate the best one
-   Run `src/main_analysis.R` to fit the selected model and generate results
