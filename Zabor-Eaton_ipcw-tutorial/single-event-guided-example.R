# This code is for the single dataset in the guided example section for single events

# Load needed R packages
library(tibble)
library(tidyr)
library(dplyr)
library(survival)
library(purrr)
library(broom)

# Generate single simulated dataset ----
# Fix parameters
n <- 500 
alpha <- 0.05
x_prop <- 0.5
a <- 2
sigma <- 500
beta <- log(0.25)
lambda <- 0.01
phi <- -5

# Set the seed
set.seed(20240429)

# Generate the data
Y1 <- rexp(n, rate = 1)
Y2 <- rexp(n, rate = 1)

Z <- rgamma(n, shape = alpha, rate = 1)

W1 <- (1 + Y1 / Z)^(-alpha)
W2 <- (1 + Y2 / Z)^(-alpha)

x <- rbinom(n, 1, prob = x_prop)

S <- sigma * (-log(1 - W1) / exp(beta * x))^(1 / a)

C <- rexp(n, rate = lambda * exp(phi * W2))

t <- pmin(S, C)
delta <- ifelse(S <= C, 1, 0)

dat <- tibble(S, t, delta, x, W2)

# Note that in the paper, x was renamed trt and W2 was renamed anc, in line with the motivating example
# Here we keep them in their original forms for use with the functions

# Save the single simulated dataset 
# Note you may need to include an appropriate path for where you want this saved
save(dat, file = "single_example_dat.rda")


# Prepare the data for analysis ----
# Ordered vector of observed times
times <- sort(unique(dat$t[dat$delta == 0]))

# Add censoring indicator, id, and tstart
dat_prep <- 
  dat |> 
  mutate(
    censor = 1 - delta,
    id = row_number(),
    tstart = 0
  )

# Apply survSplit() to event indicator
dat_long1 <- 
  survSplit(
    dat_prep,
    cut = times,
    end = "t",
    start = "tstart",
    event = "delta"
  ) |> 
  arrange(id, t) 

# Apply survSplit() to censoring indicator
dat_long2 <- 
  survSplit(
    dat_prep,
    cut = times,
    end = "t",
    start = "tstart",
    event = "censor"
  ) |> 
  arrange(id, t)

# Add censor from dat_long2 into dat_long1 
# Rename t as tstop
dat_long0 <-
  dat_long1 |> 
  select(-censor) |> 
  add_column(censor = dat_long2$censor) |> 
  rename(tstop = t)

# Fit Cox model for censoring
# Note that in the paper, W2 was renamed anc
cens_mod <- coxph(Surv(tstart, tstop, censor) ~ W2,
                  data = dat_long0, timefix = FALSE)

# Estimate IPC weights
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
  exp(-predict(cens_mod,  
               newdata = dat_long3[is.na(dat_long3$inv_wgt),], 
               type = 'expected'))

dat_long3$wgt <- 1 / dat_long3$inv_wgt

dat_long <- 
  dat_long0 |> 
  full_join(
    dat_long3 |> 
      select(id, tstop, wgt),
    by = c("id", "tstart" = "tstop")
  )

# Save the single simulated dataset in long format with weights
# Note you may need to include an appropriate path for where you want this saved
save(dat_long, file = "single_example_ipcw_dat.rda")


# Conduct the analyses ---- 
# Plot the IPCW K-M estimates
survfit2(Surv(tstart, tstop, delta) ~ x, 
         data = dat_long, weights = wgt, timefix = FALSE) |> 
  ggsurvfit() + 
  xlim(c(0, 1000)) +
  scale_color_hue(labels = c("Chemotherapy", "TKI")) +
  labs(
    x = "Days from treatment start",
    y = "Progression-free survival probability"
  )

# Fix the Cox regression model
# Note that in the paper, x was renamed trt
ipcw_cox_fit <- coxph(Surv(tstart, tstop, delta) ~ x + cluster(id), 
                      data = dat_long, weights = wgt, timefix = FALSE)


# Get the bootstrap percentile intervals ----
# Set the number of bootstrap samples
B <- 500

# Set a seed for replication
set.seed(20240917) 

# Allocate 100 cores to run in parallel
# Note that I raw the below two sections of code on a server with 100 cores
# Set up any parallel computing as appropriate for your own computing resources
# Below
# plan(multicore(workers = 100))

# Take B bootstrap samples and estimate the IPC weights for each
ipcw_boot_dat1 <- 
  map(
    1:B,
    ~ slice_sample(dat, prop = 1, replace = T)
  ) |> 
  furrr::future_map(
    ~ get_ipcw_wgt(.)
  )

# Fit the IPCW K-M curves to each of the 500 bootstrap samples
ipcw_km_boot_fit1 <- 
  ipcw_boot_dat1 |> 
  map(~get_ipcw_km_prob_x(., pre_times = seq(0, 2429, 1)))

# Fit the weighted Cox model for each of the 500 bootstrap samples
ipcw_cox_boot_fit1 <-
  ipcw_boot_dat1 |> 
  map(~ get_ipcw_cox_fit(., "wgt"))

# Extract the B log HRs 
ipcw_boot_fit1_lhr <- 
  map_dbl(ipcw_cox_boot_fit1, "log_hr")

# Then calculate the bootstrap variance
ipcw_boot_var1 <- sum((ipcw_boot_fit1_lhr - mean(ipcw_boot_fit1_lhr))^2) / B

# Calculate the percentile interval
ipcw_boot_pci1 <- quantile(ipcw_boot_fit1_lhr, c(0.025, 0.975))

# Make a table showing the results for both log(HR) and HR scale
ipcw_cox_res_tab <- 
  tibble(
    Scale = c("log(HR)", "HR"),
    Estimate = c(ipcw_cox_fit[["coefficients"]], 
                 exp(ipcw_cox_fit[["coefficients"]])),
    `Lower 95% PI` = c(ipcw_boot_pci1[[1]], exp(ipcw_boot_pci1[[1]])),
    `Upper 95% PI` = c(ipcw_boot_pci1[[2]], exp(ipcw_boot_pci1[[2]]))
  )

gt(ipcw_cox_res_tab) |> 
  fmt_number(decimals = 2) |>
  tab_footnote(
    footnote = "PI = bootstrap percentile interval",
    locations = cells_column_labels(columns = c(`Lower 95% PI`, `Upper 95% PI`))
  ) |> 
  opt_footnote_marks("standard") |> 
  tab_header("IPCW Cox hazard ratio estimates, with bootstrap percentile intervals")
