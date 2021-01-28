# Load/install needed packages ------------------------------------------------
# remotes::install_github("anneae/cwcens") # run once initially to install the cwcens packages
library(cwcens)
library(tidyverse)
library(survival)
library(survminer)
library(gtsummary)
library(icenReg)
library(flexsurv)


# Load functions --------------------------------------------------------------
# You will need to set these file paths appropriately to wherever you have saved the files
load("fn-create-single-dataset.R") 
load("fn-fit-mods.R")
load("fn-tidy-ic_sp.R")


# Generate many datasets for a given scenario ---------------------------------
# Note that for 1000 this took 18 seconds on my machine
S <- 1000 # number of simulated datasets
set.seed(01202021) # set seed

datalist <- map(1:S, 
                ~makesimdat(hr12 = 1.5, 
                            hr13 = 1.5,
                            hr23 = 1, 
                            mos_z0 = 6,
                            mos_z1 = 6,
                            n = 100,
                            dv = FALSE)
                )


# Then fit the models to the datasets -----------------------------------------
# Note that this took 815 seconds (~13.5 minutes) on my machine
modfitlist <- 
  datalist %>% 
  # First we need to code the various endpoints
  map(~ .x %>% 
        mutate(
          rfs_event = ifelse(dstatus == 1 | state2obs != Inf, 1, 0),
          rfs_months = case_when(
            rfs_event == 1 ~ pmin(state2obs, dtime), 
            rfs_event == 0 ~ laststate1),
          rfs_months_2 = case_when(
            rfs_event == 1 ~ pmin(state2obs, dtime), 
            rfs_event == 0  ~ dtime), # for this RFS endpoint, patients are censored at dtime
          left_time = case_when(
            state2obs != Inf ~ laststate1, # last known disease-free for recurrence
            state2obs == Inf & dstatus == 1 ~ dtime, # death is uncensored
            state2obs == Inf & dstatus == 0 ~ laststate1 # last known disease-free for no event
          ),
          left_time_2 = case_when(
            state2obs != Inf ~ laststate1, # last known disease-free for recurrence
            state2obs == Inf & dstatus == 1 ~ dtime, # death is uncensored
            state2obs == Inf & dstatus == 0 ~ dtime # dtime for no event
          ),
          right_time = case_when(
            state2obs != Inf ~ state2obs, # observed recurrence time for recurrence
            state2obs == Inf & dstatus == 1 ~ dtime, # death is uncensored
            state2obs == Inf & dstatus == 0 ~ Inf # no event is right censored
          ),
          # CS Recurrence
          rcr_event = ifelse(state2obs != Inf, 1, 0),
          rcr_months = case_when(
            rcr_event == 1 ~ state2obs,
            rcr_event == 0 & dstatus == 1 ~ dtime,
            rcr_event == 0 & dstatus == 0 ~ laststate1
          ),
          # CS Cox death
          death_event_1 = ifelse(dstatus == 1 & state2obs == Inf, 1, 0),
          death_time_1 = case_when(
            death_event_1 == 1 ~ dtime,
            death_event_1 == 0 & state2obs != Inf ~ state2obs,
            death_event_1 == 0 & state2obs == Inf ~ dtime
          ),
          # CS SC Cox recurrence
          left_time2 = laststate1, # last known disease free 
          right_time2 = case_when(
            state2obs != Inf ~ state2obs, # observed recurrence time for recurrence
            TRUE ~ Inf
          )
        )
  ) %>% 
  # Then fit the models to each dataset
  map(~fitsimmod(df = .x))


# Prepare data to summarize --------------------------------------------------
# This makes the list into a dataframe
res_df <- 
  do.call("rbind", modfitlist) %>% 
  rename(log_hr = estimate) 


# Check whether there were any errors or warnings
table(res_df$issue, useNA = 'ifany')


# Boxplot result -------------------------------------------------------------
res_box <- 
  res_df %>% 
  select(-issue) %>% 
  mutate(
    truth = log(1.5),
    bias = log_hr - truth) %>%
  mutate(
    model = fct_relevel(model, "Cox RFS - lastobs1", "Cox RFS - dtime",
                        "Interval RFS - lastobs1", "Interval RFS - dtime",
                        "CS Cox death", 
                        "CS Cox recurrence", "CS Interval recurrence")
  )

ggplot(res_box, aes(x = model, y = bias)) + 
  geom_boxplot() +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  labs(
    x = "",
    y = "Bias (log(HR) scale)"
  ) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# Table result ---------------------------------------------------------------
res_tbl <- 
  res_df %>% 
  mutate(
    truth = log(1.5),
    bias = log_hr - truth) %>% 
  group_by(model) %>% 
  summarize(avg_log_hr = mean(log_hr),
            sd_log_hr = sd(log_hr),
            avg_bias = mean(bias),
            sd_bias = sd(bias)) %>% 
  ungroup() %>%
  mutate(
    model = fct_relevel(model, "Cox RFS - lastobs1", "Cox RFS - dtime",
                        "Interval RFS - lastobs1", "Interval RFS - dtime",
                        "CS Cox death", 
                        "CS Cox recurrence", "CS Interval recurrence")
  )

knitr::kable(res_tbl, digits = 3)
