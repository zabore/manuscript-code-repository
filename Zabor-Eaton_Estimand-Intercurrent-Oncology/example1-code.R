# Load the needed libraries ---
library(survival)
library(cmprsk)
library(tibble)
library(dplyr)
library(patchwork)
library(purrr)
library(survminer)
library(ggsurvfit)

# Generate the simulated data ---
n <- 500 # sample size

set.seed(22) # seed for randomly generated data, for reproducibility

enroll_day <- runif(n, 0, 365 * 4) # enrollment is uniformly distributed over 4 years (Jan 1, 2018-Dec 31, 2021)
age <- rnorm(n, 65, 10) # age is normally distributed with mean 65 and standard deviation of 10
trt <- rbinom(n, 1, p = ifelse(enroll_day < 790, .4, .9))  # treatment is generated from a binomial distribution with 40% probability of reduced intensity pre-Covid, and 90% probability of reduced intensity post-Covid (Covid = March 1, 2020)

# Cancer death day is exponentially distributed, and the rate increases with age
cancer_death <- enroll_day + rexp(n , rate = exp(-10 + .0012 * age ^ 2))

# Covid death is exponentially distributed, and the rate increases with age - risk of Covid death starts after day 790
covid_death <- pmax(enroll_day, 790) + rexp(n, rate =  exp(-10 + .0015 * age ^ 2))

cancer_death <- cancer_death - enroll_day
covid_death <- covid_death - enroll_day
cens <- 365 * 4 - enroll_day # censored on December 31, 2021

time <- pmin(cancer_death, covid_death, cens)
event <- 1*(time == cancer_death) + 2 * (time == covid_death)

dat <- tibble(
  time = time,
  event = event,
  trt = trt,
  age = age,
  agesq = age^2
)


# K-M plot of overall survival with COVID-19 deaths included -------------------
# This is Figure 1A

# First, get the location of crosses for covid deaths
# first get model fit
fit_onep <- survfit(Surv(time/365, event>0)~trt, data = dat)

# vector of covid event times
covid_event_times <- dat$time[dat$event == 2]

# then make dataframe of resulting times and survival probs
# only keep times where there were events
# then match the times to the covid event times with rounding to find times where covid deaths happened
# only keep the times where a covid event happened
# 196 because there are 196 unique covid event times
dat_onep <- tibble(
  time = fit_onep[["time"]],
  surv = fit_onep[["surv"]],
  event = ifelse(fit_onep[["n.event"]] > 0, 1, 0)
  ) %>%
  filter(event == 1) %>%
  mutate(
    covid = round(time, 10) %in% round(covid_event_times/365, 10)
  ) %>%
  filter(covid) %>%
  mutate(type = "COVID-19 death")

one <- ggsurvplot(
  fit_onep,
  censor = F,
  xlim = c(0, 2),
  legend.labs = c('Standard schedule','Reduced frequency schedule'),
  xlab = "Years",
  legend.title = "",
  palette = "Set1",
  # palette = ezfun::ccf_palette("contrast")[c(1, 5)],
  font.y = 12,
  font.x = 12,
  font.legend = 12,
  break.time.by = 0.5,
  risk.table = T,
  tables.theme = theme_cleantable(),
  risk.table.y.text = F
)

# Modify the plot but not the risktable
one$plot <- one$plot +
  geom_point(data = dat_onep, aes(x = time, y = surv, shape = type)) +
  scale_shape_manual(values = 3) +
  theme(legend.title = element_blank())


# Kaplan-Meier plot with COVID-19 deaths censored ------------------------------
# This is Figure 1B
two <- ggsurvplot(
  survfit(Surv(time/365, event==1)~trt, data = dat),
  censor = F,
  xlim = c(0, 2),
  legend.labs = c('Standard schedule','Reduced frequency schedule'),
  xlab = "Years",
  legend.title = "",
  legend = "none",
  palette = "Set1",
  # palette = ezfun::ccf_palette("contrast")[c(1, 5)],
  font.y = 12,
  font.x = 12,
  font.legend = 12,
  break.time.by = 0.5,
  risk.table = T,
  tables.theme = theme_cleantable(),
  risk.table.y.text = F
)


# K-M plot with IPCW -----------------------------------------------------------
# This is Figure 1C

### Plot KM curves dealing with death due to COVID-19 using IPCW
### First, make a long dataset (counting process format)  with separate rows for
### before and after COVID, so you can fit a model for censoring stratified on covid "era"
### as a time-varying stratification factor
before_covid <- data.frame(id = 1:n, enroll_day, agesq = age^2, time, event, trt,
                           covid = 0,
                           start = 0, end = pmin(pmax(0, 790-enroll_day), time),
                           cens_cov = (event==2 & time <= pmax(0, 790-enroll_day)),
                           death_cancer= (event == 1 & time <= pmax(0, 790-enroll_day)))

after_covid <- data.frame(id = 1:n, enroll_day, agesq = age^2, time, event, trt,
                          covid = 1,
                          start = pmin(pmax(0, 790-enroll_day), time), end =  time,
                          cens_cov = (event==2 & time > pmax(0, 790-enroll_day)),
                          death_cancer= (event == 1 & time > pmax(0, 790-enroll_day)))
dat_long<- rbind(before_covid, after_covid)
dat_long<-dat_long[dat_long$end>dat_long$start, ]
dat_long<-dat_long[order(dat_long$id, dat_long$start),]

# Since there is no censoring due to COVID before COVID, you don't actually need to fit
# a stratified model. Before COVID, there are no events, so you can just fit a model
# to the after-COVID data.
# mod <- coxph(Surv(start, end, cens_cov)~agesq + strata(covid), data=dat_long)
mod <- coxph(Surv(start, end, cens_cov)~agesq, data=dat_long[dat_long$covid==1,])

# Next, make a super-long dataset with a row for each person, for each event time.
# This is needed because
times <- sort(unique(dat_long$end))
data.long <- survSplit(dat_long, cut = times, end = "end",start = "start",
                       event = "death_cancer")

data.long$p_nocens <- NA
data.long$p_nocens_tot <- NA
for(i in 1:nrow(data.long))
{
  # print(i)
  if  (data.long$covid[i] == 0) data.long$p_nocens[i] <- 1
  else {
    sfiCZ <- survfit(mod, newdata = data.long[i,])
    data.long$p_nocens[i] <- summary(sfiCZ,times = data.long$end[i])$surv/summary(sfiCZ,times = data.long$start[i])$surv
    #   KZti <- c(KZti, summary(sfiCZ,times = dat_long$end[i])$surv/summary(sfiCZ,times = dat_long$start[i])$surv)
  }
}

# Multiply the chance of not being censored due to COVID during each time interval
# to get the probability of not being censored from time zero to the end of the interval
IDs <- unique(data.long$id)
for (i in IDs){
  whichIDs <- which(data.long$id == i)
  datai <- data.long[whichIDs,]
  data.long$p_nocens_tot[whichIDs] <- cumprod(datai$p_nocens)
}

# Weight is 1/(probability still being followed)
data.long$wc <- 1/data.long$p_nocens_tot
summary(data.long$wc)

# Fit the weighted model
fit_four <- survfit(Surv(start/365, end/365, death_cancer)~ trt,
                    weight = wc, data=data.long)

four <- ggsurvplot(
  fit_four,
  censor = F,
  xlim = c(0, 2),
  legend.labs = c('Standard schedule','Reduced frequency schedule'),
  xlab = "Years",
  legend.title = "",
  legend = "none",
  palette = "Set1",
  # palette = ezfun::ccf_palette("contrast")[c(1, 5)],
  font.y = 12,
  font.x = 12,
  font.legend = 12,
  break.time.by = 0.5,
  risk.table = T,
  tables.theme = theme_cleantable(),
  risk.table.y.text = F
)



# Combined plot for manuscript ---
((one$plot | two$plot) / (four$plot | plot_spacer())) / guide_area() + 
  plot_layout(guides = "collect", heights = c(4, 4, 1)) +
  plot_annotation(tag_levels = 'A', tag_suffix = ')')
