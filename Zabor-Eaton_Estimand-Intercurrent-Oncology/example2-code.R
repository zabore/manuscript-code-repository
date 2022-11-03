### Updating so that people enroll at a constant rate over the calendar year 2019 
### After MArch 2020, 25% of people continue to attend visits every 8 weeks, 
### 75% start to attend every 24 weeks. 

# Load needed libraries ---
library(cwcens)
library(survival)
library(survminer)
library(patchwork)


# Generate the simulated data ---

# Generate visit schedule every 8 weeks, on days scale
visit.schedule.8 <- 7 * seq(8, 104, 8)

# Set the seed for reproducibility
set.seed(108)

dat75 <- NULL

#### Treatment group 
for (i in 1:500){
  if (rbinom(1,1, .25)==1) visit.schedule <- visit.schedule.8
  else {
    enroll<- runif(1, min=0, max = 7*48) # When in calendar year 2019 they enrolled
    ttcovid <- (7*48-enroll) + 7*10 # followup in 2019 + 10 weeks of 2020  before covid
    visit.schedule <- c(visit.schedule.8[visit.schedule.8<ttcovid]) # pre covid visits
    last.normal <- max(0, tail(visit.schedule, 1))
    if (last.normal+7*24 < 7*104) visit.schedule <- 
      c(visit.schedule, seq(from = last.normal+7*24, to = 7*104, by = 7*24))
  }
  # Follow everyone until 104 weeks
  dat <- simdat(1, visit.schedule = visit.schedule, vital.lfu = c(7*105), 
                scale12 = (1/0.5)*1/16e-04, scale13 = (1/0.5)*1/2e-04)
  dat$grp <- 'trt'
  dat$lastnormal<-last.normal
  maxvis <- max(dat$nvisits, dat75$nvisits)
  if (i > 1){
    if (dat$nvisits<maxvis){
      dat[,paste('t',(dat$nvisits+1):maxvis, sep = '')]<-NA
    }
    if (max(dat75$nvisits)<maxvis){
      dat75[,paste('t',(dat75$nvisits+1):maxvis, sep = '')]<-NA
    }
  }
  dat75<-rbind(dat75, dat[,c('grp','dtime','lastnormal', 'state2obs','laststate1','dstatus', 'nvisits', 
                             paste('t', 1:maxvis, sep = ''))])
}

#### Add control group
for (i in 1:500){
  if (rbinom(1,1, .25)==1) visit.schedule <- visit.schedule.8
  else {
    enroll<- runif(1, min=0, max = 7*48) # When in calendar year 2019 they enrolled
    ttcovid <- (7*48-enroll) + 7*10 # followup in 2019 + 10 weeks of 2020  before covid
    visit.schedule <- c(visit.schedule.8[visit.schedule.8<ttcovid]) # pre covid visits
    last.normal <- max(0, tail(visit.schedule, 1))
    if (last.normal+7*24 < 7*104) visit.schedule <- 
      c(visit.schedule, seq(from = last.normal+7*24, to = 7*104, by = 7*24))
  }
  dat <- simdat(1, visit.schedule = visit.schedule, vital.lfu = c(7*105), 
                scale12 = 1/16e-04)
  dat$grp <- 'con'
  dat$lastnormal<-last.normal
  maxvis <- max(dat$nvisits, dat75$nvisits)
  #  if (i > 1){
  if (dat$nvisits<maxvis){
    dat[,paste('t',(dat$nvisits+1):maxvis, sep = '')]<-NA
  }
  if (max(dat75$nvisits)<maxvis){
    dat75[,paste('t',(dat75$nvisits+1):maxvis, sep = '')]<-NA
  }
  # }
  dat75<-rbind(dat75, dat[,c('grp','dtime','lastnormal','state2obs','laststate1','dstatus', 'nvisits', 
                             paste('t', 1:maxvis, sep = ''))])
}
dat75$nvisits<-max(dat75$nvisits)

# KM curves ,  censor at last disease free
dat75$time <- pmin(dat75$dtime, dat75$state2obs)
dat75$event <- ifelse(dat75$state2obs<Inf, 1, dat75$dstatus)
dat75$time[dat75$event==0]<-dat75$laststate1[dat75$event == 0]
dat75$pctdelay<-75


# Make Figure 2 ---
mycolors <- RColorBrewer::brewer.pal(4, 'Set1')[3:4]

# Figure 2A
fit1 <- survfit(Surv(time/7, event)~grp, data=dat75)

ex2p1 <- 
  ggsurvplot(
    fit1, 
    censor = F, 
    xlim = c(0, 80),
    legend.labs = c('Standard of care','New treatment'),
    xlab = "Weeks from randomization",
    ylab = "Progression-free survival probability",
    legend.title = "",
    palette = mycolors,
    font.y = 12,
    font.x = 12,
    font.legend = 12, 
    break.x.by = 8
  )

# Figure 2B
ex2p2v2 <- 
  dat75 %>% 
  select(grp, t1:t13) %>% 
  pivot_longer(
    cols = t1:t13,
    names_prefix = "t",
    names_to = "visit_number",
    values_to = "visit_time"
  ) %>% 
  filter(!is.na(visit_time)) %>% 
  mutate(
    visit_weeks = visit_time / 7
  ) %>% 
  filter(visit_weeks <= 80) %>% 
  ggplot(aes(x = visit_weeks, fill = grp)) +
  geom_histogram(bins = 80, show.legend = FALSE) +
  facet_grid(rows = vars(grp)) +
  scale_fill_manual(values = mycolors) +
  xlab("Weeks from randomization") + 
  ylab("Count") +
  scale_x_continuous(breaks = seq(8, 80, 8)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank()
  )

# Figure 2C
# the variables starting with "t" are the actual times of the visits for each participant, up to 13 visits
# I think the times are in days
library(dplyr)
library(tidyr)

for (i in 1:nrow(dat75)){
  dat75$eight[i] <- (sum(dat75[i,paste('t',1:dat75$nvisits[1],sep = '')] >= 6*7 & 
                           dat75[i,paste('t',1:dat75$nvisits[1],sep = '')] <= 10*7  , na.rm = T) >=1)
  if (dat75$state2obs[i] < 7*6 | dat75$dtime[i] < 7*6) dat75$eight[i] <- NA
  
  dat75$sixteen[i] <- (sum(dat75[i,paste('t',1:dat75$nvisits[1],sep = '')] >= 14*7 & 
                             dat75[i,paste('t',1:dat75$nvisits[1],sep = '')] <= 18*7  , na.rm = T) >=1)
  if (dat75$state2obs[i] < 14*6 | dat75$dtime[i] < 14*6) dat75$sixteen[i] <- NA
  
  dat75$twentyfour[i] <- (sum(dat75[i,paste('t',1:dat75$nvisits[1],sep = '')] >= 22*7 & 
                                dat75[i,paste('t',1:dat75$nvisits[1],sep = '')] <= 26*7  , na.rm = T) >=1)
  if (dat75$state2obs[i] < 7*22 | dat75$dtime[i] < 7*22) dat75$twentyfour[i] <- NA
  
  dat75$thirtytwo[i] <- (sum(dat75[i,paste('t',1:dat75$nvisits[1],sep = '')] >=30*7 & 
                               dat75[i,paste('t',1:dat75$nvisits[1],sep = '')] <= 34*7  , na.rm = T) >=1)
  if (dat75$state2obs[i] < 7*30 | dat75$dtime[i] < 7*30) dat75$thirtytwo[i] <- NA
  
  dat75$forty[i] <- (sum(dat75[i,paste('t',1:dat75$nvisits[1],sep = '')] >=38*7 & 
                           dat75[i,paste('t',1:dat75$nvisits[1],sep = '')] <= 42*7  , na.rm = T) >=1)
  if (dat75$state2obs[i] < 7*38 | dat75$dtime[i] < 7*38) dat75$forty[i] <- NA
  
  dat75$fortyeight[i] <- (sum(dat75[i,paste('t',1:dat75$nvisits[1],sep = '')] >=46*7 & 
                                dat75[i,paste('t',1:dat75$nvisits[1],sep = '')] <= 50*7  , na.rm = T) >=1)
  if (dat75$state2obs[i] < 7*46 | dat75$dtime[i] < 7*46) dat75$fortyeight[i] <- NA
  
  dat75$fiftysix[i] <- (sum(dat75[i,paste('t',1:dat75$nvisits[1],sep = '')] >=54*7 & 
                              dat75[i,paste('t',1:dat75$nvisits[1],sep = '')] <= 58*7  , na.rm = T) >=1)
  if (dat75$state2obs[i] < 7*54 | dat75$dtime[i] < 7*54) dat75$fiftysix[i] <- NA
  
  dat75$sixtyfour[i] <- (sum(dat75[i,paste('t',1:dat75$nvisits[1],sep = '')] >=62*7 & 
                               dat75[i,paste('t',1:dat75$nvisits[1],sep = '')] <= 66*7  , na.rm = T) >=1)
  if (dat75$state2obs[i] < 7*62 | dat75$dtime[i] < 7*62) dat75$sixtyfour[i] <- NA
  
  dat75$seventytwo[i] <- (sum(dat75[i,paste('t',1:dat75$nvisits[1],sep = '')] >=70*7 & 
                                dat75[i,paste('t',1:dat75$nvisits[1],sep = '')] <= 74*7  , na.rm = T) >=1)
  if (dat75$state2obs[i] < 7*70 | dat75$dtime[i] < 7*70) dat75$seventytwo[i] <- NA
  
  dat75$eighty[i] <- (sum(dat75[i,paste('t',1:dat75$nvisits[1],sep = '')] >=78*7 & 
                            dat75[i,paste('t',1:dat75$nvisits[1],sep = '')] <= 82*7  , na.rm = T) >=1)
  if (dat75$state2obs[i] < 7*78 | dat75$dtime[i] < 7*78) dat75$eighty[i] <- NA
  
}

dat2 <- 
  dat75 %>% 
  select(grp, eight, sixteen, twentyfour, thirtytwo, forty, fortyeight, 
         fiftysix, sixtyfour, seventytwo, eighty) %>% 
  pivot_longer(
    cols = -grp,
    names_to = "visit_time",
    values_to = "visited"
  ) %>% 
  mutate(
    grp = recode(grp, con = "Standard of care", trt = "New treatment")
  ) %>% 
  filter(!is.na(visited)) %>% 
  group_by(grp, visit_time) %>% 
  summarize(
    n_visits = sum(visited),
    n_at_risk = n()
  ) %>% 
  mutate(
    pct_visits = n_visits / n_at_risk,
    visit_weeks = case_when(
      visit_time == "eight" ~ 8,
      visit_time == "sixteen" ~ 16,
      visit_time == "twentyfour" ~ 24,
      visit_time == "thirtytwo" ~ 32,
      visit_time == "forty" ~ 40,
      visit_time == "fortyeight" ~ 48,
      visit_time == "fiftysix" ~ 56,
      visit_time == "sixtyfour" ~ 64,
      visit_time == "seventytwo" ~ 72,
      visit_time == "eighty" ~ 80
    ))


ex2p2 <- 
  ggplot(dat2, aes(x = visit_weeks, y = pct_visits, fill = grp)) + 
  geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
  xlab("Weeks from randomization") + 
  ylab("Visit completion rate") + 
  scale_fill_manual(values = mycolors) + 
  theme_classic() + 
  ylim(c(0, 1)) + 
  scale_x_continuous(breaks = seq(8, 80, 8))


# Combine the plots ---
(ex2p1$plot + ex2p2v2 / ex2p2) / guide_area() + 
  plot_layout(guides = "collect", heights = c(4, 1)) +
  plot_annotation(tag_levels = 'A', tag_suffix = ')')