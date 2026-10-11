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

------

```
R version 4.3.2 (2023-10-31)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS 26.2

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Asia/Seoul
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] readxl_1.4.3       here_1.0.1         ggdensity_1.0.0    egg_0.4.5          gridExtra_2.3      ggplot2_3.5.2     
[7] dplyr_1.1.4        rstan_2.32.6       StanHeaders_2.32.9

loaded via a namespace (and not attached):
 [1] generics_0.1.4     magrittr_2.0.3     grid_4.3.2         RColorBrewer_1.1-3 cellranger_1.1.0   rprojroot_2.0.4   
 [7] jsonlite_2.0.0     pkgbuild_1.4.8     Formula_1.2-5      viridisLite_0.4.2  QuickJSR_1.2.2     scales_1.4.0      
[13] textshaping_0.3.7  codetools_0.2-19   abind_1.4-8        cli_3.6.5          rlang_1.1.6        withr_3.0.2       
[19] tools_4.3.2        inline_0.3.19      parallel_4.3.2     colorspace_2.1-1   curl_5.2.1         vctrs_0.6.5       
[25] R6_2.6.1           matrixStats_1.3.0  stats4_4.3.2       lifecycle_1.0.4    V8_4.4.2           car_3.1-3         
[31] MASS_7.3-60        ragg_1.2.7         pkgconfig_2.0.3    RcppParallel_5.1.7 pillar_1.11.0      gtable_0.3.6      
[37] loo_2.7.0          glue_1.8.0         Rcpp_1.1.0         systemfonts_1.2.3  tibble_3.3.0       tidyselect_1.2.1  
[43] rstudioapi_0.15.0  farver_2.1.2       labeling_0.4.3     carData_3.0-5      compiler_4.3.2    
```
