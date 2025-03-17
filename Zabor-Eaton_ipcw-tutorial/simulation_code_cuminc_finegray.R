library(survival)
library(ggplot2)
library(dplyr)
library(gtsummary)
library(cmprsk)

####### Function to simulate data

sim_data <- function(n = 100, censoring = 'none'){ 
  if (! (censoring %in% c('none','independent','baseline','time-varying'))) stop("censoring must be one of none, independent, baseline, or time-varying.")
  psa0 <- rnorm( n, mean = 1.86, sd = 0.64) ## From Kattan paper sims
  z1 <- cut(psa0, 
            breaks = c(-Inf, qnorm(.25, mean = 1.86, sd = 0.64) , 
                       qnorm(.5, mean = 1.86, sd = 0.64),
                           qnorm(.75, mean = 1.86, sd = 0.64), Inf),
            labels = 0:3)

  # For each person, time-var PSA is their baseline PSA + patient-specific slope*time
  slope <- rnorm(n, mean = 0.01, sd = 0.15) 
  psa1 <- psa0 + slope*1
  psa2 <- psa0 + slope*2
  psa3 <- psa0 + slope*3
  psa4 <- psa0 + slope*4
  psa5 <- psa0 + slope*5
  psa6 <- psa0 + slope*6
  psa7 <- psa0 + slope*7
  psa8 <- psa0 + slope*8
  psa9 <- psa0 + slope*9
  
  ## Back transform a uniform random variable to get the event time and type 
  U <- runif(n)
  lp <- ifelse(z1 == 0, 0, ifelse(z1 ==1,beta11_sh, ifelse(z1 == 2, beta12_sh,beta13_sh)))  
  sdf <- function(x, lp) 1-(1-p*(1-exp(-x)))^exp(lp)
  event_type <- ifelse(U < sapply(lp, function(lp) sdf(Inf, lp)), 1, 2) # 1 means the person had a type 1 event (mets), 2 means death
  time <-  sapply(1:n, function(i) ifelse(event_type[i] == 1, 
                    uniroot(function(x) sdf(x, lp[i]) - U[i], c(0, 100), extendInt = 'upX')$root, rexp(1, exp(-1 ))))
  ## Add censoring times 
  if (censoring == 'independent')  {
      cens <- rexp(n, rate = 1/8)
  }
    if (censoring == 'baseline')  { 
      # censoring depends on psa quartile at BL 
      cens <- rexp(n, rate = (1/4)*exp(-1*(as.numeric(as.character(z1))-1.5)))
    }
    if (censoring == 'time-varying')  { 
      # censoring depends on time-varying PSA, measured every year
      cens[cens>1] <- 1 + rexp(n, rate = (1/4)*exp(-2*(psa1 - 1.5)))[cens>1]
      cens[cens>2] <- 2 + rexp(n, rate = (1/4)*exp(-2*(psa2 - 1.5)))[cens>2]
      cens[cens>3] <- 3 + rexp(n, rate = (1/4)*exp(-2*(psa3 - 1.5)))[cens>3]
      cens[cens>4] <- 4 + rexp(n, rate = (1/4)*exp(-2*(psa4 - 1.5)))[cens>4]
      cens[cens>5] <- 5 + rexp(n, rate = (1/4)*exp(-2*(psa5 - 1.5)))[cens>5]
      cens[cens>6] <- 6 + rexp(n, rate = (1/4)*exp(-2*(psa6 - 1.5)))[cens>6]
      cens[cens>7] <- 7 + rexp(n, rate = (1/4)*exp(-2*(psa7 - 1.5)))[cens>7]
      cens[cens>8] <- 8 + rexp(n, rate = (1/4)*exp(-2*(psa8 - 1.5)))[cens>8]
      cens[cens>9] <- 9 + rexp(n, rate = (1/4)*exp(-2*(psa9 - 1.5)))[cens>9]
    }
  
  event_type <- ifelse(time<cens, event_type, 0)
  time <-  ifelse(time<cens, time, cens)
  
  # Treat time to death and mets as competing risks
  event_type <- factor(event_type, 0:2, labels=c("censor", "mets", "death"))
  
  time[time == 0] <- min(time[time > 0])/2 # Sometimes the simulated data has time exactly zero 
    
  dat<- data.frame(z1 = z1,  event_type=event_type, time = time, 
                   psa0, psa1, psa2, psa3, psa4, psa5, psa6, psa7, psa8, psa9)
  dat$id <- 1:nrow(dat)
  return(dat) 
}

####### Functions to get cumulative incidence estimates 

# Function - naive approach, ignore dependent censoring 
cuminc_naive <- function(dat, esttimes){
  fit <- survfit(Surv(time, event_type) ~ 1, data=dat)
  est <- summary(fit, times = esttimes, extend = T)
  est$pstate[,2]
}

# Function - estimate cumulative incidence separately in each stratum, take weighted average 
cuminc_waverage <- function(dat, esttimes = seq(from = 0 , to =10 , length.out = 100)){
  strat <- survfit(Surv(time, event_type) ~ strata(z1), data=dat)
  est <- summary(strat, times = esttimes, extend = T)
  pstate <- matrix(est$pstate[,2], nrow = length(esttimes), ncol = 4)
  weighted <- (pstate[,1]*sum(dat$z1 == 0) + pstate[,2]*sum(dat$z1 == 1)+
                 pstate[,3]*sum(dat$z1 == 2)+pstate[,4]*sum(dat$z1 == 3))/nrow(dat)
  # Do not provide an estimate unless there is at least one person still
  # being followed in each group. 
  weighted[esttimes >= min(max(dat$time[dat$z1 == 0]), 
                               max(dat$time[dat$z1 == 1]),
                                   max(dat$time[dat$z1 == 2]),
                                       max(dat$time[dat$z1 == 3]))] <- NA
  weighted
}

# Function - make wide data into long data in order to fit a Cox model for censoring
wide_to_long <- function(dat){
  dat$censored <-  (dat$event_type == 'censor')
  dat$event2_time <- ifelse(dat$event_type =='death', dat$time, NA)
  # Get weights for each person and time - code to calculate weights is from Willem
  dat$Tstart <- 0
  times <- sort(unique(dat$time[dat$censored == T]))
  data.long <- survSplit(dat,
                         cut = times,
                         end = "time",
                         start = "Tstart",
                         event = "event_type")
  data.long <- data.long[order(data.long$id,data.long$time),]
  data.long.cens <- survSplit(dat,
                              cut = times,
                              end = "time",
                              start = "Tstart",
                              event = "censored")
  data.long.cens <- data.long.cens[ order(data.long.cens$id,data.long.cens$time),]
  data.long$censored <- data.long.cens$censored # data.long.cens has the correct "censored" indicator for fitting the model to estimate the censoring distribution 
  data.long
}

# Function - add IPCW weights 
add_ipcw_weights <- function(data.long, strat = 'no', new_data = NULL){
  # If no new data is given, we want to estimate the probability of not
  # being censored between time zero and the start of an interval in data.long.
  no_new_data <- ifelse(is.null(new_data), T, F)
  if (no_new_data){
    new_data <- data.long
    new_data$time <- new_data$Tstart
    new_data$Tstart <- 0 
  }
  # strat = 'yes' means weights come from estimating the censoring distribution in each stratum 
  if (strat == 'yes'){
    # Fit a model for the censoring distribution 
    fit_cens <- list(survfit(Surv(Tstart, time, censored) ~ 1,
                      data = data.long[data.long$z1 == 0,], timefix = F),
                     survfit(Surv(Tstart, time, censored) ~ 1,
                      data = data.long[data.long$z1 == 1,], timefix = F),
                     survfit(Surv(Tstart, time, censored) ~ 1,
                      data = data.long[data.long$z1 == 2,], timefix = F),
                     survfit(Surv(Tstart, time, censored) ~ 1,
                      data = data.long[data.long$z1 == 3,], timefix = F))
    # Add weights
    KZti <- NULL
    for (j in 1:nrow(new_data)){
      if (new_data$time[j]==0) KZti <- c(KZti,1)
      else KZti <- c(KZti, summary(fit_cens[[as.numeric(as.character(new_data$z1[j]))+1]], new_data$time[j], extend = T)$surv/
        summary(fit_cens[[as.numeric(as.character(new_data$z1[j]))+1]], new_data$Tstart[j], extend = T)$surv)
    }
    new_data$p_notcens <- KZti
#    new_data <- new_data[new_data$p_notcens >0,]
  }
  # strat = 'no' means you estimate the censoring distribution with a Cox model 
  if (strat == 'no'){
    # Fit a model for the censoring distribution 
    fit_cens <- coxph(Surv(Tstart, time, censored) ~ z1,
                    data = data.long, timefix = F)
    # Add weights and get cuminc estimates
    new_data$p_notcens[new_data$time == 0] <- 1
    new_data$p_notcens[is.na(new_data$p_notcens)] <-
      exp(-predict(fit_cens, newdata =
                     new_data[is.na(new_data$p_notcens),], 
                 type = 'expected'))
  }
  if (no_new_data) {
    data.long$p_notcens <- new_data$p_notcens
    return(data.long)
  }
  return(new_data)
}

# Function - get ipcw weighted estimates of cumulative incidence using 
# probability of yet not being censored already stored in the dataset: p_notcens 
cuminc_ipcw <- function(data.long, esttimes = seq(from = 0 , to =10 , length.out = 100), extend = T){
    fit <- survfit(Surv(Tstart, time, event_type) ~ 1,
                 data=data.long,  weights=1/p_notcens, id =id, timefix = F)
    ipcw <- summary(fit, times = esttimes, extend = T)
    
    # Do not provide an estimate unless there is at least one person still
    # being followed in each group. 
    if (extend == F){
      ipcw$pstate[,2][esttimes >= 
                          min(max(data.long$time[data.long$z1 == 0]), 
                               max(data.long$time[data.long$z1 == 1]),
                               max(data.long$time[data.long$z1 == 2]),
                               max(data.long$time[data.long$z1 == 3]))] <- NA
    }
    ipcw$pstate[,2]
   }

######## Functions to do FG regression  

# Function - 
fg_split <- function(data.long){
  # Split the data at any time an event of type 2 happened, bc weights need to reflect the probability of being censored after event 2 time
 # times <- unique(data.long$time[data.long$event_type=='death'])
#  data.long <- survSplit(data.long,
#                         cut = times,
#                         end = "time",
#                         start = "Tstart",
#                         event = "censored")
  # Add rows AFTER an event of type 2, until the last event time in the data, 
  # split at every censoring time
#  times <- unique(data.long$time[data.long$event_type=='censor'])
  times <- sort(unique(data.long$time[data.long$censored == 1]))
  for (id in unique(data.long$id[data.long$event_type =='death'])){
    death_time <- data.long$time[data.long$event_type == 'death'&
                                   data.long$id == id]
    data.long <- rbind(data.long,
                  data.frame(
                    z1 = data.long$z1[data.long$id ==id][1],
                    psa0 = data.long$psa0[data.long$id ==id][1],
                    psa1 = data.long$psa1[data.long$id ==id][1],
                    psa2 = data.long$psa2[data.long$id ==id][1],
                    psa3 = data.long$psa3[data.long$id ==id][1],
                    psa4 = data.long$psa4[data.long$id ==id][1],
                    psa5 = data.long$psa5[data.long$id ==id][1],
                    psa6 = data.long$psa6[data.long$id ==id][1],
                    psa7 = data.long$psa7[data.long$id ==id][1],
                    psa8 = data.long$psa8[data.long$id ==id][1],
                    psa9 = data.long$psa9[data.long$id ==id][1],
                    id = id, 
                    censored = 0, 
                    event2_time = data.long$event2_time[data.long$id ==id][1],
                    Tstart = c(death_time,times[times>death_time]),
                    time = c(times[times>death_time], Inf), 
                    event_type = 'censor'))
  }
  data.long[order(data.long$id,data.long$Tstart),]
}

# Function - add fg weights based on a Cox model for censoring distribution (strat = 'no') 
# or KM curves of the censoring distribution in each stratum (strat = 'yes')
add_fg_weights <- function(data.long, strat = 'no'){
  # When we add the weights, we need the time variables to be named Tstart and time
  temp <- data.frame(Tstart = data.long$event2_time,
                   time = data.long$Tstart, 
                   censored = NA,
                   z1 = data.long$z1)
  
  temp$p_notcens_after_death[is.na(temp$Tstart)] <- 1
  temp$p_notcens_after_death[temp$Tstart >= temp$time] <- 1
  if (strat == 'yes'){
    temp$p_notcens_after_death[is.na(temp$p_notcens_after_death)] <- 
      add_ipcw_weights(data.long[data.long$Tstart<data.long$event2_time | 
                                   is.na(data.long$event2_time), ], strat = 'yes', 
                       new_data = temp[is.na(temp$p_notcens_after_death),])$p_notcens
  }
  if (strat == 'no'){
    temp$p_notcens_after_death[is.na(temp$p_notcens_after_death)] <- 
      add_ipcw_weights(data.long[data.long$Tstart<data.long$event2_time| 
                                   is.na(data.long$event2_time), ], strat = 'no', 
                       new_data = temp[is.na(temp$p_notcens_after_death),])$p_notcens
  }
  data.long$p_notcens_after_death <- temp$p_notcens_after_death
  data.long
}


### Function - naive FG SHR estimation - assume censoring is independent 
fg_naive <- function(dat){
  pdata <- finegray(Surv(time, event_type) ~ . , data=dat, timefix = F)
    m1 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ z1,
                     weight=fgwt, data=pdata, timefix = F)
    summary(m1)$coef[,c(1,3)]
}

# Function - weighted FG - weighted estimates with weights that are already in the dataset: p_notcens_after_death
fg_weighted <- function(data.long, extend = T){
  if (extend == F) {
    # Do not use data unless there is at least one person still
    # being followed in each group. 
    tau <- min(max(data.long$time[data.long$z1 == 0 &
                                    data.long$time<=data.long$event2_time]), 
               max(data.long$time[data.long$z1 == 1 &
                                    data.long$time<=data.long$event2_time]),
               max(data.long$time[data.long$z1 == 2 &
                                    data.long$time<=data.long$event2_time]),
               max(data.long$time[data.long$z1 == 3 &
                                    data.long$time<=data.long$event2_time]))
    data.long$event_type[data.long$time > tau] <- 'censor'
    data.long$time[data.long$time > tau] <- tau
    data.long <- data.long[data.long$time > data.long$Tstart,]
  }
  m1 <- coxph(Surv(Tstart, time, event_type=='mets') ~ z1 + cluster(id),
        weight = p_notcens_after_death, 
        data=data.long[data.long$p_notcens_after_death>0,], timefix = F)
  summary(m1)$coef[,c(1,4)]
}

########## Set true parameters for simulating data 

# SH for metastasis will follow a proportional SH model 
# Effects of psa on SH for event 1 (metastasis)
beta11_sh = log(1.5); beta12_sh = log(2.25); beta13_sh = log(3.4) 
# Proportion of type 1 events when z1 = 0
p <- 0.3

# True cumulative incidence of mets - overall, in all PSA groups combined - over time
sdf <- function(x, lp) 1-(1-p*(1-exp(-x)))^exp(lp)
true_cuminc <- function(t) .25*sdf(t, 0) + .25*sdf(t, beta11_sh) +
  .25*sdf(t, beta12_sh) + .25*sdf(t, beta13_sh) 
truth <- true_cuminc(seq(from = 0, to = 10, length.out=100))

##### Run simulations for cumulative incidence
nsim <- 500
res_naive <- matrix(NA, ncol = 100, nrow = nsim)
res_weighted_average <- matrix(NA, ncol = 100, nrow = nsim)
res_ipcw <- matrix(NA, ncol = 100, nrow = nsim) 
res_ipcw_cox <- matrix(NA, ncol = 100, nrow = nsim)
  
set.seed(1989)
for (i in 1:nsim){
  print(i)
  dat <- sim_data(500, censoring = 'baseline')
  data.long  <-  wide_to_long(dat)
  
  res_naive[i,] <- cuminc_naive(dat, seq(from = 0, to =10, length.out = 100))
  res_weighted_average[i,] <- cuminc_waverage(dat, seq(from = 0, to =10, length.out = 100))
  res_ipcw[i,] <- cuminc_ipcw(add_ipcw_weights(data.long, strat = 'yes'), seq(from = 0, to =10, length.out = 100), extend = F)
  res_ipcw_cox[i,] <- cuminc_ipcw(add_ipcw_weights(data.long, strat = 'no'), seq(from = 0, to =10, length.out = 100))
}

##### Run simulations for Fine and Gray regression 

nsim<- 500
set.seed(2345)
z11 <- matrix(NA, nrow = nsim, ncol = 10)
z12 <- matrix(NA, nrow = nsim, ncol = 10)
z13 <- matrix(NA, nrow = nsim, ncol = 10)
for (j in 1:nsim){
  print(j)
  dat <- sim_data(500, censoring = 'baseline')
  data.long  <-  wide_to_long(dat)
  data.long.fg <- fg_split(data.long)
  
  naive <- fg_naive(dat)
  z11[j,1:2] <- naive[1,]
  z12[j,1:2] <- naive[2,]
  z13[j,1:2] <- naive[3,]
  
  strat <- fg_weighted(add_fg_weights(data.long.fg, strat = 'yes'), extend = F)
  z11[j,3:4] <- strat[1,]
  z12[j,3:4] <- strat[2,]
  z13[j,3:4] <- strat[3,]
  
  cox <- fg_weighted(add_fg_weights(data.long.fg, strat = 'no'))
  z11[j,5:6] <- cox[1,]
  z12[j,5:6] <- cox[2,]
  z13[j,5:6] <- cox[3,]

  adata <- finegray(Surv(time, event_type) ~ . + strata(z1),
             data=dat, na.action=na.pass)
  mod <- coxph(Surv(fgstart, fgstop, fgstatus) ~ z1 + cluster(id),
                     weight=fgwt, data=adata)
  z11[j,7:8] <- summary(mod)$coef[1,c(1,4)]
  z12[j,7:8] <- summary(mod)$coef[2,c(1,4)]
  z13[j,7:8] <- summary(mod)$coef[3,c(1,4)]
  
  mod2 <- crr(dat$time, as.numeric(dat$event)-1, 
              cov1 = model.matrix(~dat$z1)[,-1], cengroup = dat$z1)
  z11[j,9:10] <- c(mod2$coef[1], sqrt(mod2$var[1,1]))
  z12[j,9:10] <- c(mod2$coef[2], sqrt(mod2$var[2,2]))
  z13[j,9:10] <- c(mod2$coef[3], sqrt(mod2$var[3,3]))
}