# Load the basket package
library(basket)

# Load some additional packages
library(dplyr)
library(purrr)
library(furrr)
library(tictoc)

# Load the function
source(here::here("code", "heterogeneous_response_sim_fn.R"))

# Set a seed for this session
set.seed(19158)

# Make a tibble of input parameters - include the prior PEP in this!
ins <- 
  tibble(
    pi1 = rep(seq(0.5, 0.1, -0.1), 27),
    pi2 = rep(rep(0.5, 5), 27),
    pi3 = rep(seq(0.5, 0.9, 0.1), 27),
    sampsize = rep(rep(c(10, 20, 30), each = 5), 9),
    prior_pep = rep(seq(0.1, 0.9, 0.1), each = 15)
  )

# map over the number of simulations with pmap to run many times across each set of inputs
nsim <- 500

# make this parallel for speed using the furrr package
plan(multiprocess)

tic()
sim_res <-
  future_map(1:nsim, 
        ~ pmap(
          ins,
          exact_3basket
          ))
toc()

save(
  sim_res, 
  file = here::here(
    "results", 
    "heterogeneous_response_sim_res.rda"
    )
  )
  
