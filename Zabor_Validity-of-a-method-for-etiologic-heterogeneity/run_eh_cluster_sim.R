library(dplyr)
library(purrr)

# Source the simulation function
setwd("D:/GitHub/manuscript-code-repository/Zabor_Validity-of-a-method-for-etiologic-heterogeneity")
source("./eh_cluster_sim_fun.R")

# Fix the simulation parameters
M <- 3
N <- 2000 
pi <- c(0.4, 0.2, 0.2, 0.2) 
P <- 2
mu_m <- matrix(c(1.5, 0.75, 0, 0, 0.75, 1.5), ncol = 3, byrow = T) 
K_A <- 15 
a <- 2.1
lambda_Am <- matrix(rep(c(a, 0, a, 0, a), times = c(5, 15, 5, 15, 5)), ncol = 3) 
J <- 3
K_B <- 15 
b <- 1.275
lambda_Bj <- matrix(rep(c(b, 0, b, 0, b), times = c(5, 15, 5, 15, 5)), ncol = 3) 
K_C <- 70
nstart <- 1000
reducedim <- TRUE


# iterate over the function nsim times using the specified parameters
nsim <- 1000
set.seed(20171116)
simres <- rerun(nsim, eh_cluster_sim(M, N, pi, P, mu_m, K_A, lambda_Am, 
                                     J, K_B, lambda_Bj, K_C, nstart, reducedim))


# The result in simres is a list of length nsim containing: 
# max_d: the maximum D 
# min_misclass: the minimum misclassification rate
# unique_g: the number of unique solutions from k-means clustering
# max_d_univd: the maximum D for sequentially increasing subsets of tumor markers after ranking by univariate D
# min_misclass_univd: the minimum misclassification rate for sequentially increasing subsets of tumor markers after ranking by univariate D
# n_ka_sel: number of tumor markers in Y_A selected in dimension reduction
# n_kb_sel: number of tumor markers in Y_B selected in dimension reduction
# n_kc_sel: number of tumor markers in Y_C selected in dimension reduction
