# Ocular Pathologies in DataSUS

Bayesian hierarchical analysis of ocular disease hospitalization rates across all 27 Brazilian states (2010–2024), using DATASUS/SIH data. Three disease categories are analyzed: glaucoma, retinopathy, and general eye & appendage diseases.

## Reproducing the analysis

All scripts assume the working directory is the repo root.

### Step 1 — Fetch population denominators
```r
source("src/fetch_ibge_projections.R")
# Writes: data/ibge_population_projections.csv
```

### Step 2 — Initial model comparison (NB0–NB8)
```r
source("src/choose_model.R")
# Writes: diagnostics/{disease}/model_comparison.csv + fit.rds per model
```

### Step 3 — Extended model comparison (19 variants)
```r
source("ad_hoc/src/run_all.R")
# Writes: ad_hoc/results/{variant}/{disease}/ (fits, diagnostics, charts)
#         ad_hoc/results/cross_model_summary.csv
```

### Step 4 — Apply model selection rule
```r
source("ad_hoc/src/select_best_model.R")
# Writes: ad_hoc/results/selected_model.csv
```

### Step 5 — Final analysis (selected model: NB5-AgeCat)
```r
source("src/main_analysis.R")
# Writes: results/ (all manuscript figures, tables, posterior plots)
```

### Step 6 — Sensitivity analysis
```r
source("src/sensitivity_analysis.R")
# Writes: sensitivity/ (phi prior robustness, exclude-2024, missing-as-zero)
```

### Step 7 — Out-of-sample validation (hold out 2023–2024)
```r
source("ad_hoc/src/holdout_validation.R")
# Writes: ad_hoc/results_holdout/ (predictions vs observed, naive baselines)
```

## Selected model

**NB5-AgeCat**: Negative Binomial with dynamic state intercepts (Gaussian random walk) and categorical age effects (sum-to-zero deflection coding across 7 age strata).

$$\log\lambda_{s,t,a} = \log P_{s,t,a} + \alpha_{s,t} + \beta_{\text{age},a}$$

$$\alpha_{s,t} \mid \alpha_{s,t-1} \sim \mathcal{N}(\alpha_{s,t-1},\, \sigma_\alpha), \quad \sum_a \beta_{\text{age},a} = 0$$

Selected from 19 candidate models by the prespecified rule: filter on convergence ($\hat{R} < 1.01$, ESS $> 400$), then lowest LOO-IC among converged models, with parsimony as a tiebreaker. The categorical age structure improved LOO-IC by **1700–4300 points** and reduced MAE by **65–77%** over all linear-age models.

## Repository structure

```
data/                          Raw data + IBGE population projections
models/                        Base Stan models NB0–NB8
src/                           Main analysis pipeline
  choose_model.R               Model selection (NB0–NB8)
  data_processing.R            Data preparation helper
  model_selection.R            Diagnostic functions
  fetch_ibge_projections.R     IBGE API download
  main_analysis.R              Final analysis (runs after model selection)
  sensitivity_analysis.R       Prior and data sensitivity checks
ad_hoc/
  models/                      Extended Stan models (19 variants)
  src/
    run_all.R                  Fit all 19 variants
    select_best_model.R        Apply prespecified selection rule
    holdout_validation.R       Out-of-sample validation (2023–2024)
  MODELS.md                    Full model descriptions with LaTeX notation
diagnostics/                   Per-model fits and diagnostics (gitignored)
sensitivity/                   Sensitivity analysis outputs
results/                       Final manuscript outputs (figures, tables)
eda/                           Exploratory data analysis notebook
```

## Key outputs

| File | Description |
|------|-------------|
| `results/fig1_temporal_trends_{disease}.png` | Regional hospitalization rate trends 2010–2024 |
| `results/fig2_age_stratified_{disease}.png` | Age-stratified trends by region |
| `results/forecasts/fig3_{disease}_{state}.png` | State-level forecasts 2025–2036 |
| `results/posterior_age_effects.png` | Posterior of categorical age offsets β_age |
| `results/posterior_sigma_rw.png` | Posterior of RW variance σ_α |
| `results/posterior_phi.png` | Posterior of NegBin overdispersion φ |
| `results/table1a_model_comparison.csv` | Full LOO-IC comparison across 19 models |
| `results/table1b_demographic_performance.csv` | Model performance by age and region |
| `ad_hoc/results/cross_model_summary.csv` | Compact summary of all 19 variants |
| `sensitivity/phi_prior_sensitivity.csv` | Prior robustness for φ |
| `sensitivity/scenario_comparison.csv` | Year-exclusion and missingness sensitivity |
