library(ggsurvfit)
library(purrr)

beta11_sh <- log (1.5)
beta12_sh <- log (2.25)
beta13_sh <- log (3.4)
p <- 0.3

set.seed(9843)
dat <- sim_data_CR(500, censoring = 'baseline')
dat_long  <-  wide_to_long_CR(dat)

# Plot cumulative incidence by PSA quartile 
plot(survfit(Surv(t, delta) ~ z1, data=dat), col = 1:4, lty = rep(1:2, each = 4))

# Plot censoring distribution 
plot(survfit(Surv(t, delta=='censor') ~ z1, data=dat))

# Get bootstrap samples to use for SE estimation 
set.seed (20240917)
ipcw_boot_dat <- map(1:500 , ~ slice_sample(dat , prop = 1, replace = T))
ipcw_boot_dat <- map(ipcw_boot_dat, function(x) {x[,'id'] <- 1:nrow(x); return(x)})
ipcw_boot_dat_long <- map(ipcw_boot_dat, wide_to_long_CR)

# Estimate cumulative incidence with weighted average, standard error 
w_average <- cuminc_waverage(dat, sort(dat$t))
w_average_boot <- ipcw_boot_dat |> map(~ cuminc_waverage(.,  sort(dat$t)))
w_average_boot <- matrix(unlist(w_average_boot), ncol = length(dat$t), byrow = T)
sd_w_average <- apply(w_average_boot, 2, sd)
lower_w_average <- apply(w_average_boot, 2, function(x) ifelse(sum(is.na(x))>0,NA, quantile(x, prob = .025)))
upper_w_average <- apply(w_average_boot, 2, function(x) ifelse(sum(is.na(x))>0,NA, quantile(x, prob = .975)))

# Estimate cumulative incidence with weights from a Cox model, standard error 
cox <- cuminc_ipcw(add_ipcw_weights(dat_long, strat = 'no'), sort(dat$t))
cox_boot <- ipcw_boot_dat_long |> map(~ cuminc_ipcw(add_ipcw_weights(., strat = 'no'), sort(dat$t)))
cox_boot <- matrix(unlist(cox_boot), ncol = length(dat$t), byrow = T)
sd_cox <- apply(cox_boot, 2, sd)
lower_cox <- apply(cox_boot, 2, function(x) ifelse(sum(is.na(x))>0,NA, quantile(x, prob = .025)))
upper_cox <- apply(cox_boot, 2, function(x) ifelse(sum(is.na(x))>0,NA, quantile(x, prob = .975)))

to_plot <- rbind(data.frame(time = sort(dat$t), est = w_average, method = 'w_average', lty='a', gp=1),
                 data.frame(time = sort(dat$t), est = lower_w_average, method = 'w_average', lty='b', gp=2),
                 data.frame(time = sort(dat$t), est = upper_w_average, method = 'w_average', lty='b', gp=3),
                 data.frame(time = sort(dat$t), est = cox, method = 'cox', lty='a', gp=4),
                 data.frame(time = sort(dat$t), est = lower_cox, method = 'cox', lty='b', gp=5),
                 data.frame(time = sort(dat$t), est = upper_cox, method = 'cox', lty='b', gp=6))
ggplot(to_plot) + geom_step( mapping=aes(x=time, y=est, colour = method, linetype = lty, group = gp))+
  xlim(0,5)+xlab('Time')+ ylab('Cumulative incidence of mets') +ylim(0,.65)+ 
  theme_bw(base_size = 14) + theme(legend.position="bottom") + 
  scale_color_hue(labels = c("Cox model IPCW","Simple weighted/\n non-parametric IPCW"))+
  theme(legend.title=element_blank())+guides(color = guide_legend(reverse = TRUE), linetype = 'none')

# Histogram of non-parametric weights
data.strat <- add_ipcw_weights(dat_long, strat = 'yes')
ggplot(data.strat, aes(x=1/p_notcens)) + geom_histogram()

# Weight over time by PSA group
fit <- summary(survfit(Surv(t, delta=='censor') ~ z1, data=dat))
fit <- data.frame(time = fit$time[fit$surv>0], weight = 1/fit$surv[fit$surv>0], z1 = fit$strata[fit$surv>0])
ggplot(fit) + geom_step( mapping=aes(x=time, y=weight, group = z1, colour = z1))+ xlab('Time')+ ylab('Weight') + theme_bw(base_size = 14) + theme(legend.position="bottom") + scale_color_hue(labels = c("PSA Q1", "PSA Q2","PSA Q3","PSA Q4"))+theme(legend.title=element_blank())

# Histogram of Cox model weights
data.cox <- add_ipcw_weights(dat_long, strat = 'no')
ggplot(data.cox, aes(x=1/p_notcens)) + geom_histogram()

# Fine gray estimates 
dat_long_fg <- fg_split(dat_long)
ipcw_boot_dat_fg <- map(ipcw_boot_dat_long, fg_split)

## Stratified
fg_strat <- fg_weighted(add_fg_weights(dat_long_fg, strat = 'yes'), extend = F)
fg_strat_boot <- map(ipcw_boot_dat_fg, ~fg_weighted(add_fg_weights(., strat = 'yes'), extend = F)[,1])
fg_strat_boot <- matrix(unlist(fg_strat_boot), ncol = 3, byrow = T)
sd_strat_boot <- apply(fg_strat_boot, 2, sd)
lower_fg_strat <- apply(fg_strat_boot, 2, function(x) ifelse(sum(is.na(x))>0,NA, quantile(x, prob = .025)))
upper_fg_strat <- apply(fg_strat_boot, 2, function(x) ifelse(sum(is.na(x))>0,NA, quantile(x, prob = .975)))

exp(fg_strat[,1])
exp(lower_fg_strat)
exp(upper_fg_strat)

# Cox
fg_cox <- fg_weighted(add_fg_weights(dat_long_fg, strat = 'no'))
fg_cox_boot <- map(ipcw_boot_dat_fg, ~fg_weighted(add_fg_weights(., strat = 'no'))[,1])
fg_cox_boot <- matrix(unlist(fg_cox_boot), ncol = 3, byrow = T)
sd_cox_boot <- apply(fg_cox_boot, 2, sd)
lower_fg_cox <- apply(fg_cox_boot, 2, function(x) ifelse(sum(is.na(x))>0,NA, quantile(x, prob = .025)))
upper_fg_cox <- apply(fg_cox_boot, 2, function(x) ifelse(sum(is.na(x))>0,NA, quantile(x, prob = .975)))

exp(fg_cox[,1])
exp(lower_fg_cox)
exp(upper_fg_cox)

# Plot the weights for two random individuals in each stratum 
data.strat <- add_fg_weights(dat_long_fg, strat = 'yes')

set.seed(454)
ids <- c(sample(unique(data.strat$id[data.strat$p_notcens_after_death!=1 & data.strat$z1 == 0]), 2),
         sample(unique(data.strat$id[data.strat$p_notcens_after_death!=1 & data.strat$z1 == 1]), 2),
         sample(unique(data.strat$id[data.strat$p_notcens_after_death!=1 & data.strat$z1 == 2]), 2),
         sample(unique(data.strat$id[data.strat$p_notcens_after_death!=1 & data.strat$z1 == 3]), 2))
ggplot(data.strat[data.strat$id %in% ids,]) + 
  geom_step( mapping=aes(x=tstart, y=p_notcens_after_death, group = as.factor(id), 
                         colour = as.factor(z1))) + xlab('Time')+ ylab('Weight')  + 
  theme_bw(base_size = 14) + theme(legend.position="bottom") + 
  scale_color_hue(labels = c("PSA Q1", "PSA Q2","PSA Q3","PSA Q4"))+theme(legend.title=element_blank())

