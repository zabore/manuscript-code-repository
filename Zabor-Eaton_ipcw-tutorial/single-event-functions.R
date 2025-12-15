# Function to get the IPCW ----
## In this function the variable names are all fixed - event, censor, tstart, time
## This is something that could be made more flexible later
## data = dataset obtained from gen_sim_data() function
## This version only returns the original unstabilized wgt:
### wgt = unstabilized weight based on pooled data
get_ipcw_wgt <- function(data) {
  
  times <- sort(unique(data$t[date$delta == 0]))
  
  dat_prep <- 
    data |> 
    mutate(
      censor = 1 - delta,
      tstart = 0,
      id = row_number()
    )
  
  dat_long1 <- 
    survSplit(
      dat_prep,
      cut = times,
      end = "t",
      start = "tstart",
      event = "delta"
    ) |> 
    arrange(id, t) 
  
  dat_long2 <- 
    survSplit(
      dat_prep,
      cut = times,
      end = "t",
      start = "tstart",
      event = "censor"
    ) |> 
    arrange(id, t)
  
  dat_long0 <-
    dat_long1 |> 
    select(-censor) |> 
    add_column(censor = dat_long2$censor) |> 
    rename(tstop = t)
  
  cens_mod <- coxph(Surv(tstart, tstop, censor) ~ W2,
                    data = dat_long0, timefix = FALSE)
  
  dat_long3 <- 
    dat_long0 |> 
    mutate(
      # to estimate the probability of not being censored between time zero and the start of an interval (tstart)
      tstop = tstart,
      tstart = 0,
      # the weight is one for intervals from 0 to 0
      inv_wgt = case_when(
        tstop == 0 ~ 1
      )
    )
  
  dat_long3$inv_wgt[is.na(dat_long3$inv_wgt)] <-
    exp(-predict(cens_mod, newdata =
                   dat_long3[is.na(dat_long3$inv_wgt),], type = 'expected'))
  
  dat_long3$wgt <- 1 / dat_long3$inv_wgt
  
  dat_long <- 
    dat_long0 |> 
    full_join(
      dat_long3 |> 
        select(id, tstop, wgt),
      by = c("id", "tstart" = "tstop")
    )
  
  return(dat_long)
}


# Function to get the IPCW Cox regression result ----
# Return the log hazard, SE, HR, and 95% CI
# Note that the confidence intervals returned are based on the robust SE
## data = dataset obtained from gen_sim_data() function
# Using truncated stabilized weights now
get_ipcw_cox_fit <- function(data, weight) {
  ipcw_cox_fit <- coxph(Surv(tstart, tstop, delta) ~ x + cluster(id), 
                        data = data, 
                        weights = data[[weight]], 
                        timefix = FALSE)
  full_join(
    tidy(ipcw_cox_fit)[, 1:4] |>
      rename(log_hr = estimate, log_hr_se = std.error, log_hr_rob_se = robust.se),
    tidy(ipcw_cox_fit, exponentiate = TRUE, conf.int = TRUE)[, c(1, 2, 7, 8)] |>
      rename(hr = estimate, hr_ci_low = conf.low, hr_ci_high = conf.high),
    by = "term"
  )
}


# Function to get the IPCW Kaplan-Meier estimated probabilities, by binary covariate ----
## data = dataset obtained from get_ipcw()function
# Using original wgt
get_ipcw_km_prob_x <- function(data, pre_times = seq(0, 50, 1)) {
  ipcw_km_surv_fit <- survfit(Surv(tstart, tstop, delta) ~ x, data = data,
                              weights = wgt, timefix = FALSE)
  tibble(
    time = summary(ipcw_km_surv_fit, times = pre_times)$time, 
    surv = summary(ipcw_km_surv_fit, times = pre_times)$surv, 
    x = summary(ipcw_km_surv_fit, times = pre_times)$strata
  )
}
