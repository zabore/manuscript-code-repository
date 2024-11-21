###################################################################################
# This code is used for the simulation studies
# It conducts 1000 simulations all 9 missing data patterns and 2 sample size
# It generates original data, introduces missing data, bootstraps, imputes, 
## summarizes imputation success and regression success
## calculates performance metrics, predicted values and bias
###################################################################################


# Load needed libraries------------------------------------------------------------
library(mvtnorm)
library(dplyr)
library(tibble)
library(purrr)
library(survival)
library(riskRegression)
library(future)
library(furrr)
library(tidyverse)
library(stats)


# Source simulation functions -----------------------------------------------------
# Loading it from my local Github clone
source("D:/manuscript-code-repository/Mi-Zabor_bootstrap-impute-predict-tutorial/simulation-functions.R")


# Set simulation parameters -------------------------------------------------------
# Number of simulated datasets
nsim <- 1000

# List of missing variable patterns
missing_var_list <- list(c('x1'),
                         c('x1'),
                         c('x1'),
                         c('x1','x3','x4'),
                         c('x1','x3','x4'),
                         c('x1','x3','x4'),
                         c('x1','x3','x4','x7','x10','x11'),
                         c('x1','x3','x4','x7','x10','x11'),
                         c('x1','x3','x4','x7','x10','x11'))

# List of missing percentages
missing_pct_list <- list(c(0.05),
                         c(0.15),
                         c(0.6),
                         c(0.05,0.05,0.05),
                         c(0.05,0.15,0.3),
                         c(0.15,0.3,0.6),
                         c(0.05,0.05,0.05,0.05,0.05,0.05),
                         c(0.05,0.05,0.15,0.15,0.3,0.3),
                         c(0.15,0.15,0.3,0.3,0.6,0.6))


# For nsim datasets for multiple scenarios ----------------------------------------
set.seed(20241003)

# generate the full datasets - separately by sample size
dat_750_nsim <- 
  map(
    1:nsim,
    ~ sim_sample(750)
  )

dat_3500_nsim <- 
  map(
    1:nsim,
    ~ sim_sample(3500)
  )

# then generate missing data sets into a list - separately by sample size
missing_list_750_nsim <-
  map(
    dat_750_nsim,
    ~ map2(
      missing_var_list,
      missing_pct_list,
      function(a, b) missing_dat(.x, a, b)
    )
  )

missing_list_3500_nsim <-
  map(
    dat_3500_nsim,
    ~ map2(
      missing_var_list,
      missing_pct_list,
      function(a, b) missing_dat(.x, a, b)
    )
  )

# then combine the two sample sizes 
# this is a nested list of length nsim
# and the second level is a list of 18 missing datasets
missing_list_nsim <- 
  map2(
    missing_list_750_nsim,
    missing_list_3500_nsim,
    ~ c(.x, .y)
  )

# Running 1000 simulations for all scenarios all together might take too much memory
# Could consider split it by scenarios
# e.g. missing_list_nsim <- lapply(missing_list_nsim, function(x) x[1])

# impute missing data to get 54 single imputation lists
# need to map over nsim (3), missing_list (length 18) and list(1, 2, 3) (length 3)
# this results in a nested list of length nsim
# the second level is length 18 (missing patterns)
# the third level is length 3 (imputation methods)
# Then we flatten it so that: the first level is length nsim
## the second level is length 54 (18 missing patterns x 3 imputation patterns = 54 scenarios)
impute_list_nsim <- 
  missing_list_nsim |> 
  map(
    ~ .x |>  map(function(b) map(list(1, 2, 3), function(a) imp_dat(b, a)))
  ) |> 
  map(~ list_flatten(.x))

# summarize imputation success

## Apply count_all_na_pred_unified to each method within the current simulation
impute_list_nsim_success <- map_dfr(impute_list_nsim, function(a) {
  map_dfr(a, count_all_na_pred_unified, .id = "method_id")
}, .id = "sim_id")

## Summarize the results to get the count of datasets (method1-3) with fully NA columns across all simulations
impute_list_nsim_success_res <- impute_list_nsim_success %>%
  pivot_longer(cols = -c(sim_id, method_id), names_to = "variable", values_to = "count") %>%
  group_by(method_id, variable) %>%
  summarize(total_na_count = sum(count), .groups = "drop")

## Save the imputation success results
save(impute_list_nsim_success_res, file = "D:/impute-bootstrap-predict/impute_list_nsim_success.rda")


# bootstrap and impute missing data to get 54 bootstrapped lists
# this results in a nested list of length nsim
# the second level is length 18 (missing patterns)
# the third level is length 3 (imputation methods)
# the fourth level is length 500 (bootstrap datasets)
# Then we flatten it so that:
## the first level is length nsim
## the second level is length 54 (18 missing patterns x 3 imputation patterns = 54 scenarios)
## the third is length B (bootstrap datasets)

B <- 500

# Parallelization
plan(multisession, workers = 80)

# Bootstrapping with future_map
boot_dat <- 
  missing_list_nsim |> 
  future_map(
    ~ .x |> map(function(df)
      map(1:B, ~ slice_sample(df, prop = 1, replace = TRUE))
    ),
    .options = furrr_options(seed = TRUE)  
  ) |> 
  map(~ list_flatten(.x))

# impute the 500 bootstrapped data
boot_impute_list_nsim <- boot_dat |> 
  future_map(
    function(sim_data) {
      map(1:3, function(method_index) {
        map(sim_data, function(boot_sample) {
          # Perform transformation and imputation in a single step
          imp_dat(boot_sample, method_index)
          })
        })
      },
    .options = furrr_options(seed = TRUE) 
    )


# summarize bootstrap imputation success
## Apply count_all_na_pred_unified to each method within the current simulation
boot_impute_list_nsim_success <- map_dfr(boot_impute_list_nsim, function(a) {
  map_dfr(a, count_all_na_pred_unified, .id = "method_id")
}, .id = "sim_id")

## Summarize the results to get the count of datasets (method1-3) with fully NA columns across all simulations
boot_impute_list_nsim_success_res <- boot_impute_list_nsim_success %>%
  pivot_longer(cols = -c(sim_id, method_id), names_to = "variable", values_to = "count") %>%
  group_by(method_id, variable) %>%
  summarize(total_na_count = sum(count), .groups = "drop")

## Save the bootstrap imputation success results
save(boot_impute_list_nsim_success_res, file = "D:/impute-bootstrap-predict/boot_impute_list_nsim_success.rda")

# calculate performance metrics and predicted values
# this results in a nested list of length nsim
# the second level is length 54
# the third is length 36 - these are the performance metrics and predictions
res_nsim <- 
  future_map2(
    impute_list_nsim,
    boot_impute_list_nsim,
    ~ map2(.x, .y, calc_perf_pred),
    .options = furrr_options(seed = TRUE)
    )

# finally calculate the bias
# this results in a nested list of length nsim
# the second level is length 54
# the third is 3 - these are the results summaries for each setting for each sim
bias_nsim <- map(
  res_nsim,
  ~ map(.x, bias_summary)
)

# Save the performance and bias results
save(bias_nsim, file = "D:/impute-bootstrap-predict/sim-bias-res.rda")
save(res_nsim, file = "D:/impute-bootstrap-predict/sim-res.rda")