
library(survival)
library(dplyr)

# Function - simulate competing risks data
sim_data_CR <- function (n = 100, censoring = 'none'){
  if (! ( censoring %in% c('none','independent','baseline'))) stop(" censoring must
be one of none, independent, or baseline")
  z1 <- as.factor(sample(0:3, n, replace = T))
#  psa0 <- rnorm( n, mean = 1.86 , sd = 0.64) ## From Kattan paper sims
#  z1 <- cut(psa0 ,
#            breaks = c(-Inf , qnorm (.25 , mean = 1.86 , sd = 0.64) ,
#                       qnorm (.5, mean = 1.86 , sd = 0.64) ,
#                       qnorm (.75 , mean = 1.86 , sd = 0.64) , Inf),
#            labels = 0:3)
  ## Back transform a uniform random variable to get the event time and type
  U <- runif(n)
  lp <- ifelse (z1 == 0, 0, ifelse (z1 ==1, beta11_sh, 
                                    ifelse (z1 == 2, beta12_sh, beta13_sh)))
  sdf <- function (x, lp) 1-(1-p*(1-exp(-x)))^exp(lp)
  delta <- ifelse (U < sapply (lp , function (lp) sdf(Inf , lp)), 1, 2) # 1 means the person had a type 1 event (mets), 2 means death
  t <- sapply (1:n, function (i) ifelse (delta[i] == 1,
                                         uniroot ( function (x) sdf(x, lp[i]) - U[i],
                                                   c(0, 100) , extendInt = 'upX')$root , rexp (1, exp(-1 ))))
  ## Add censoring times
  if (censoring == 'independent') {
    cens <- rexp(n, rate = 1/8)
  }
  
  if (censoring == 'baseline') {
    # censoring depends on psa quartile at BL
    cens <- rexp(n, rate = (1/4)*exp(-1*(as.numeric(as.character (z1)) -1.5)))
  }
  delta <- ifelse (t<cens , delta , 0)
  t <- ifelse(t<cens , t, cens)
  # Treat time to death and mets as competing risks
  delta <- factor(delta, 0:2, labels =c("censor", "event_1", "event_2"))
  t[t == 0] <- min(t[t > 0])/2 # Sometimes the simulated data has time exactly zero
  dat <- data.frame(z1 = z1 , delta=delta , t = t)
 #  dat$id <- 1: nrow(dat)
 # dat$event2_time <- ifelse (dat$delta =='death ', dat$t, NA)
  return (dat)
}

# Function - make wide CR data into long data. 
# In the wide data, t should be event or censoring time, delta should be a factor variable
# with levels "censor", "event_1", "event_2" where "censor" is censored and "event_1" is the event of interest.
# The long data will include multiple rows per person, split at every event or censoring time 
wide_to_long_CR <- function(dat){
  dat$id <- 1: nrow(dat)
  dat$censor <-  (dat$delta == 'censor')
  dat$event2_time <- ifelse(dat$delta == 'event_2', dat$t, NA)
  dat$tstart <- 0
  times <- sort(unique(dat$t[dat$censor == T]))
  data_long1 <- survSplit(dat,
                         cut = times,
                         end = "t",
                         start = "tstart",
                         event = "delta")
  data_long1 <- data_long1[order(data_long1$id, data_long1$t),]
  data_long2 <- survSplit(dat,
                              cut = times,
                              end = "t",
                              start = "tstart",
                              event = "censor")
  data_long2 <- data_long2[ order(data_long2$id,data_long2$t),]
  data_long1$censor <- data_long2$censor # data_long2 has the correct "censor" indicator for fitting the model to estimate the censoring distribution 
  data_long1 <- data_long1 |> rename(tstop = t)
  data_long1
}

# Function - add IPCW weights 
add_ipcw_weights <- function(data_long, strat = 'no', new_data = NULL, by.start = T){
  # The variable you want to base the weights on needs to be called z1, 
  # and if you are doing strat = 'yes', it should be a factor variable. 
  # You can get the probability of not being censored between time 0 and tstart 
  # or the probability of not being censored in the interval from tstart to tstop. 
  # The former is most often needed to IPCW and it's the default (by.start = T).

  # If no new data is given, get P(not cens) for data.long
  if (is.null(new_data)) new_data <- data_long
  
  # If by.start = T, we will estimate the probability of not being censored between time zero and 
  # the start of an interval, tstart
  # If by.start = F we will estimate the probability of not being censored during the interval from time tstart and tstop
  temp <- new_data[,c('z1','tstart','tstop','censor')]
  if (by.start){
    temp$tstop <- temp$tstart
    temp$tstart <- 0 
  }
  temp$p_notcens <- NA
  # strat = 'yes' means weights come from estimating the censoring distribution in each stratum 
  if (strat == 'yes'){
    # Fit a model for the censoring distribution 
    cens_mod <- lapply(levels(data_long$z1), 
                       function(x) survfit(Surv(tstart, tstop, censor) ~ 1,
                                           data = data_long[data_long$z1 == x,], timefix = F))
    # Add weights
    for (i in 1:length(levels(data_long$z1))){
      temp$p_notcens[temp$z1==levels(data_long$z1)[i]] <- 
        summary(cens_mod[[i]], temp$tstop[temp$z1==levels(data_long$z1)[i]], extend = T)$surv/
        summary(cens_mod[[i]], temp$tstart[temp$z1==levels(data_long$z1)[i]], extend = T)$surv
    }
  }
  # strat = 'no' means you estimate the censoring distribution with a Cox model 
  if (strat == 'no'){
    # Fit a model for the censoring distribution 
    cens_mod <- coxph(Surv(tstart, tstop, censor) ~ z1,
                      data = data_long, timefix = F)
    # Add weights 
    temp$p_notcens[temp$tstop == 0] <- 1
    temp$p_notcens[is.na(temp$p_notcens)] <- exp(-predict(cens_mod, newdata =
                              temp[is.na(temp$p_notcens),], type = 'expected'))
  }
  new_data$p_notcens <- temp$p_notcens
  return(new_data)
}

cuminc_naive <- function(dat, esttimes){
  fit <- survfit(Surv(t, delta) ~ 1, data=dat)
  est <- summary(fit, times = esttimes, extend = T)
  est$pstate[,2]
}

cuminc_waverage <- function(dat, esttimes = seq(from = 0 , to =10 , length.out = 100)){
  strat <- survfit(Surv(t, delta) ~ strata(z1), data=dat)
  est <- summary(strat, times = esttimes, extend = T)
  pstate <- matrix(est$pstate[,2], nrow = length(esttimes), ncol = 4)
  weighted <- (pstate[,1]*sum(dat$z1 == 0) + pstate[,2]*sum(dat$z1 == 1)+
                 pstate[,3]*sum(dat$z1 == 2)+pstate[,4]*sum(dat$z1 == 3))/nrow(dat)
  # Do not provide an estimate unless there is at least one person still
  # being followed in each group. 
  weighted[esttimes >= min(max(dat$t[dat$z1 == 0]), 
                           max(dat$t[dat$z1 == 1]),
                           max(dat$t[dat$z1 == 2]),
                           max(dat$t[dat$z1 == 3]))] <- NA
  weighted
}

cuminc_ipcw <- function(data_long, esttimes = seq(from = 0 , to =10 , length.out = 100), extend = T){
  fit <- survfit(Surv(tstart, tstop, delta) ~ 1,
                 data=data_long,  weights=1/p_notcens, id = id, timefix = F)
  ipcw <- summary(fit, times = esttimes, extend = T)
  
  # Do not provide an estimate unless there is at least one person still
  # being followed in each group. 
  if (extend == F){
    ipcw$pstate[,2][esttimes >= 
                      min(max(data.long$tstop[data.long$z1 == 0]), 
                          max(data.long$tstop[data.long$z1 == 1]),
                          max(data.long$tstop[data.long$z1 == 2]),
                          max(data.long$tstop[data.long$z1 == 3]))] <- NA
  }
  ipcw$pstate[,2]
}

fg_split <- function(data_long){
  # Add rows AFTER an event of type 2, until the last event time in the data, 
  # split at every censoring time
  times <- sort(unique(data_long$tstop[data_long$censor == 1]))

  event2_dat <- data_long[data_long$delta=='event_2',] |> select(id, z1, event2_time)
  cens_dat <- data.frame(tstart = c(0, times), tstop = c(times, Inf), censor = 0, delta = 'censor')
  
  fulldat <- merge(event2_dat, cens_dat)
  fulldat <- fulldat[fulldat$event2_time < fulldat$tstop,]
  fulldat$tstart[fulldat$event2_time > fulldat$tstart] <- fulldat$event2_time[fulldat$event2_time > fulldat$tstart]
  fulldat <- rbind(data_long, fulldat[names(data_long)])
  fulldat[order(fulldat$id, fulldat$tstart),]
}

### Add fg weights based on a Cox model for cens distribution or stratified 
# estimates of the censoring distribution 
add_fg_weights <- function(data_long_fg, strat = 'no'){
  # When we add the weights, we need the time variables to be named tstart and tstop
  temp <- data.frame(tstart = data_long_fg$event2_time,
                     tstop = data_long_fg$tstart, 
                     censor = NA,
                     z1 = data_long_fg$z1)
  temp$p_notcens_after_death[is.na(temp$tstart)] <- 1
  temp$p_notcens_after_death[temp$tstart >= temp$tstop] <- 1
    temp$p_notcens_after_death[is.na(temp$p_notcens_after_death)] <- 
      add_ipcw_weights(data_long_fg[data_long_fg$tstart<data_long_fg$event2_time | 
                                   is.na(data_long_fg$event2_time), ], strat = strat, 
                       new_data = temp[is.na(temp$p_notcens_after_death),], by.start = F)$p_notcens
    
    data_long_fg$p_notcens_after_death <- temp$p_notcens_after_death
    data_long_fg
}

fg_naive <- function(dat){
  pdata <- finegray(Surv(t, delta) ~ . , data=dat, timefix = F)
  m1 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ z1,
              weight=fgwt, data=pdata, timefix = F)
  summary(m1)$coef[,c(1,3)]
}

fg_weighted <- function(data_long_fg, extend = T){
  if (extend == F) {
    # Do not use data unless there is at least one person still
    # being followed in each group. 
    tau <- min(max(data_long_fg$tstop[data_long_fg$z1 == 0 &
                    (is.na(data_long_fg$event2_time) | 
                       data_long_fg$tstop<=data_long_fg$event2_time)]), 
               max(data_long_fg$tstop[data_long_fg$z1 == 1 &
                     (is.na(data_long_fg$event2_time) | 
                        data_long_fg$tstop<=data_long_fg$event2_time)]),
               max(data_long_fg$tstop[data_long_fg$z1 == 2 &
                      (is.na(data_long_fg$event2_time) | 
                         data_long_fg$tstop<=data_long_fg$event2_time)]),
               max(data_long_fg$tstop[data_long_fg$z1 == 3 &
                      (is.na(data_long_fg$event2_time) | 
                         data_long_fg$tstop<=data_long_fg$event2_time)]))
    data_long_fg$delta[data_long_fg$tstop > tau] <- 'censor'
    data_long_fg$tstop[data_long_fg$tstop > tau] <- tau
    data_long_fg <- data_long_fg[data_long_fg$tstop > data_long_fg$tstart,]
  }
  m1 <- coxph(Surv(tstart, tstop, delta=='event_1') ~ z1 + cluster(id),
              weight = p_notcens_after_death, 
              data=data_long_fg[data_long_fg$p_notcens_after_death>0,], timefix = F)
  summary(m1)$coef[,c(1,4)]
}


