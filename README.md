# Herbivore assemblages and nutrient regime jointly shape plant coexistence

------

This repository contains all code and data used for the analysis.

Each R script uses `here::i_am()` to declare its location within the project and `here::here()` for all input and output paths. The scripts can therefore be run from the project root or from any subdirectory within the project without changing their behavior.

* `data_raw`: contains the raw Excel sheet for the duckweed competition, herbivore preference, nutrient, and biomass measurements
* `data_processed`: contains R scripts for processing the raw data and the resulting analysis-ready CSV files
  * `data_processed_area.R`: processes duckweed surface-area measurements from the single-herbivore, herbivore-combination, and sequential-arrival experiments
  * `data_processed_eta_pref.R`: summarizes the feeding and nesting preferences of *Elophila turbata* (ETA)
* `script`: contains shared R code for preparing the processed datasets and calculating Shannon diversity
  * `script_diversity.R`: combines the experimental datasets and calculates diversity for each treatment, nutrient condition, replicate, and observation day
* `stanmodel`: contains the Stan model for estimating intrinsic growth rates, density-dependent competition coefficients, and pairwise coexistence metrics
  * `model_default.stan`: Stan model used for the population-dynamic analysis
* `stanfit`: contains the R script for fitting the population-dynamic model
  * `stanfit_default_all.R`: prepares the model data, fits the Stan model, and saves the fitted model as `stanfit_default_all.rda`
* `figure`: contains R scripts for generating the main and supplementary figures
  * `figure_data.R`: plots the observed duckweed surface-area and Shannon-diversity trajectories
  * `figure_regression_arrival_v2.R`: estimates and plots herbivore main effects, interaction effects, sequential-arrival effects, and ETA preference comparisons
  * `figure_stanfit_default_all.R`: plots posterior estimates of intrinsic growth rates and competition coefficients
  * `figure_stanfit_default_niche.R`: translates posterior model estimates into niche differences, fitness differences, and pairwise competition outcomes
* `figure_stanfit`: contains R scripts for posterior diagnostics and model-fit figures
  * `figure_stanfit_default_all_fit.R`: compares observed and model-predicted growth rates
  * `figure_stanfit_default_all_post.R`: plots posterior distributions by sampling chain

------

* Run `data_processed/data_processed_area.R` and `data_processed/data_processed_eta_pref.R` to generate the processed CSV files.
* Run `stanfit/stanfit_default_all.R` after the processed data have been generated. This script sources `script/script_diversity.R`, fits the Stan model, and saves the fitted model as `stanfit/stanfit_default_all.rda`.
* `figure/figure_data.R` and `figure/figure_regression_arrival_v2.R` can be run after data processing is complete.
* Stan-based scripts in `figure` and `figure_stanfit` require `stanfit/stanfit_default_all.rda` and should be run after the model has been fitted.
