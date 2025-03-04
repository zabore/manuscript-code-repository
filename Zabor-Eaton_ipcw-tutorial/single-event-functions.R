# Function to get the IPCW ----
## In this function the variable names are all fixed - event, censor, tstart, time
## This is something that could be made more flexible later
## data = dataset obtained from gen_sim_data() function
## This version only returns the original unstabilized wgt:
### wgt = unstabilized weight based on pooled data
get_ipcw_wgt <- function(data) {
  
  times <- sort(unique(data$t[date$t == 0]))
  
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
    mutate(id_nest = id) |> 
    nest(.by = id_nest) |> 
    mutate(
      prob_num = map(data, ~ summary(
        survfit(Surv(tstart, tstop, censor) ~ 1, timefix = FALSE, data = .x, id = id),
        times = .x[["tstart"]])$surv
      ),
      prob_den = map(data, ~ summary(
        survfit(cens_mod, newdata = .x, id = id, timefix = FALSE),
        times = .x[["tstart"]])$surv
      )
    ) |> 
    unnest(c(data, prob_num, prob_den)) |> 
    mutate(
      wgt = 1 / prob_den,
      wgt_stab = prob_num / prob_den
    ) |> 
    arrange(id, tstart)
  
  return(dat_long3)
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