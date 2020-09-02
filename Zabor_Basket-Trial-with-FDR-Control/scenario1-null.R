# Note that this file is setup to run on the HPC, submitted using scenario1-null-batch-script.txt. 

# Load needed libraries
library(purrr)
library(dplyr)
library(basket)
library(furrr)
library(binom)

# Fix the sample size vector
ss <- c(26, 25, 18, 17, 11, 8, 5, 7, 5, 19)

# Fix the true response rate vector
trr <- rep(0.1, 10)

# Fix the basket names
subs <- paste0("Basket", seq(1:10))

# Source the function to run a single simulation
source("./do_single_sim.R")

# Fix the number of simulated datasets
nsim <- 5000

# make this parallel for speed using the furrr package
plan(multiprocess(workers = 20))

# wrap the function in possibly and return NA if there is an error
try_do_single_sim <- possibly(do_single_sim, otherwise = NA)

# set a seed for reproducibility
set.seed(20190925)

# cycle of nsim and save all results to list
sim_res <- 
  future_map(1:nsim,
      ~ try_do_single_sim(ss, trr, subs)
  )

# save resulting list object to file
save(
  sim_res, 
  file = "./scenario1-null.rda"
)