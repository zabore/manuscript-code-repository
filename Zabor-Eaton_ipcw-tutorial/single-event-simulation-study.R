library(tibble)
library(tidyr)
library(dplyr)
library(survival)
library(purrr)
library(broom)
library(future)
library(furrr)

# Scenario 1


# For now this script just generates the data and then creates the weights
# I'll have a separate script to fit the various models


# Source the functions
source(here::here("code", "fns-sim-ipcw.R"))


# Now run the simulation ----
set.seed(20240425) 

M <- 500
n <- 500 
alpha <- 0.05
x_prop <- 0.5
a <- 2
sigma <- 500
beta <- log(0.25)
lambda <- 0.01
phi <- -5

# This returns a list of tibbles
sim_dat <- 
  map(
    1:M, 
    ~ gen_sim_data(n, alpha, x_prop, a, sigma, beta, lambda, phi)
  )

plan(multicore(workers = 100))

# This returns a list of tibbles
ipcw_dat <-
  furrr::future_map(
    sim_dat,
    ~get_ipcw(.x)
  )

# Real survival probability
# Returns a long tibble with columns m (simulation number), time, surv
real_surv_prob_x <-
  map(
    sim_dat,
    ~get_real_surv_prob_x(.x, pre_times = seq(0, 1000, 1))
  ) |>
  list_rbind(names_to = "m")

# Standard Kaplan-Meier
# Returns a long tibble with columns m, time, surv
km_prob_x <-
  map(
    ipcw_dat,
    ~get_km_prob_x(.x, pre_times = seq(0, 1000, 1))
  ) |>
  list_rbind(names_to = "m")

# IPCW Kaplan-Meier
# Returns a long tibble with columns m, time, surv
# This is only for wgt
ipcw_km_prob_x <-
  map(
    ipcw_dat,
    ~get_ipcw_km_prob_x(.x, pre_times = seq(0, 1000, 1))
  ) |>
  list_rbind(names_to = "m")

# Standard Cox regression
# Returns a tibble with one row per simulation number
cox_res <-
  map_df(
    ipcw_dat,
    ~get_cox_fit(.x),
    .id = "m"
  )

# IPCW Cox regression
# Returns a tibble with one row per simulation number
ipcw_cox_res_wgt <-
  map_df(
    ipcw_dat,
    ~get_ipcw_cox_fit(.x, "wgt"),
    .id = "m"
  )

# Conduct the bootstrapping ----
plan(multicore(workers = 100))
B <- 500
set.seed(20240123)

# Bootstrap the data and re-estimate the IPCW weights
# This results in a nested list
# The first level of the list is length M for the number of simulated datasets
# Nested within each m is a list of dataframes that is length B for the bootstrap samples
# Right now this is only for the unstabilized untruncated weight based on the pooled data - i.e. the original form of the weight

ipcw_boot_dat <-
  future_map(
    sim_dat, 
    function(x) 
      map(
        1:B, 
        function(y) 
          slice_sample(x, prop = 1, replace = T)
      )
  ) |> 
  future_map(
    ~map(.x, ~get_ipcw_wgt(.))
  )