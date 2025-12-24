# Performs simulations and saves results in results.txt which includes one row per simulation. 
# Each row includes the simulations number;
# naive estimates of cumulative incidence at times seq(from = 0, to =10, length.out = 100);
# weighted average estimates of cumulative incidence at times seq(from = 0, to =10, length.out = 100);
# Cox model-based IPCW estimates of cumulative incidence at times seq(from = 0, to =10, length.out = 100);
# naive FG estimates of beta_11, beta_12, beta_13, the corresponding SE estimates, 
# 95% confidence interval lower bounds, and 95% confidence interval upper bounds;
# naive FG estimates of beta_11, beta_12, beta_13, the corresponding SE estimates, 
# 95% confidence interval lower bounds, and 95% confidence interval upper bounds;
# stratified KM-based IPCW FG estimates of beta_11, beta_12, beta_13, the corresponding SE estimates, 
# 95% confidence interval lower bounds, and 95% confidence interval upper bounds;
# Cox model-based IPCW FG estimates of beta_11, beta_12, beta_13, the corresponding SE estimates, 
# 95% confidence interval lower bounds, and 95% confidence interval upper bounds.

library("foreach"); library("doParallel")
source("competing-risks-functions.R")

beta11_sh <- log (1.5)
beta12_sh <- log (2.25)
beta13_sh <- log (3.4)
p <- 0.3

ncores <- 16

cl <- makeCluster(ncores) # creates a cluster with <ncore> cores
registerDoParallel(cl) # register the cluster
clusterCall(cl, function() {library(survival); library(purrr); library(dplyr)})

seed<-567
set.seed(seed); seeds <- seed + sample(-2^25:2^25, 5000)
iter<- 1:5000

#Run simulation in parallel on <ncores> cores
res <- foreach(iter = iter, .combine = rbind) %dopar% { 
  set.seed(seeds[iter])
  
  dat <- sim_data_CR(500, censoring = 'baseline')
  dat_long  <-  wide_to_long_CR(dat)
  
  # Get bootstrap samples to use for SE estimation 
  ipcw_boot_dat <- map(1:500 , ~ slice_sample(dat , prop = 1, replace = T))
  ipcw_boot_dat <- map(ipcw_boot_dat, function(x) {x[,'id'] <- 1:nrow(x); return(x)})
  ipcw_boot_dat_long <- map(ipcw_boot_dat, wide_to_long_CR)
  
  # Estimate cuminc naive
  naive <- cuminc_naive(dat, seq(from = 0, to =10, length.out = 100))
  
  # Estimate cumulative incidence with weighted average
  w_average <- cuminc_waverage(dat, seq(from = 0, to =10, length.out = 100))
  
  # Estimate cumulative incidence with weights from a Cox model, standard error 
  cox <- cuminc_ipcw(add_ipcw_weights(dat_long, strat = 'no'), seq(from = 0, to =10, length.out = 100))
  
  # Fine gray estimates 
  dat_long_fg <- fg_split(dat_long)
  ipcw_boot_dat_fg <- map(ipcw_boot_dat_long, fg_split)
  
  ## Naive
  fg_naive_est <- fg_naive(dat)
  fg_naive_est <- c(fg_naive_est[,1],fg_naive_est[,2], fg_naive_est[,1]-1.96*fg_naive_est[,2],fg_naive_est[,1]+1.96*fg_naive_est[,2])
    
  ## Stratified
  fg_strat <- fg_weighted(add_fg_weights(dat_long_fg, strat = 'yes'), extend = F)
  fg_strat_boot <- map(ipcw_boot_dat_fg, ~fg_weighted(add_fg_weights(., strat = 'yes'), extend = F)[,1])
  fg_strat_boot <- matrix(unlist(fg_strat_boot), ncol = 3, byrow = T)
  sd_strat_boot <- apply(fg_strat_boot, 2, sd)
  lower_fg_strat <- apply(fg_strat_boot, 2, function(x) ifelse(sum(is.na(x))>0,NA, quantile(x, prob = .025)))
  upper_fg_strat <- apply(fg_strat_boot, 2, function(x) ifelse(sum(is.na(x))>0,NA, quantile(x, prob = .975)))
  fg_strat <- c(fg_strat[,1], sd_strat_boot, lower_fg_strat, upper_fg_strat)
    
  # Cox
  fg_cox <- fg_weighted(add_fg_weights(dat_long_fg, strat = 'no'))
  fg_cox_boot <- map(ipcw_boot_dat_fg, ~fg_weighted(add_fg_weights(., strat = 'no'))[,1])
  fg_cox_boot <- matrix(unlist(fg_cox_boot), ncol = 3, byrow = T)
  sd_cox_boot <- apply(fg_cox_boot, 2, sd)
  lower_fg_cox <- apply(fg_cox_boot, 2, function(x) ifelse(sum(is.na(x))>0,NA, quantile(x, prob = .025)))
  upper_fg_cox <- apply(fg_cox_boot, 2, function(x) ifelse(sum(is.na(x))>0,NA, quantile(x, prob = .975)))
  fg_cox <- c(fg_cox[,1], sd_cox_boot, lower_fg_cox, upper_fg_cox)
  
  res <- c(iter, naive, w_average, cox, fg_naive_est, fg_strat, fg_cox)
  
  write.table(res, file="results.txt",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  
} 