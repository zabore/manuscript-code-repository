###############################################################################
# This code generates a single simulated dataset
# The dataset was used in the Guided Example section
###############################################################################

# load needed libraries --------------------------------------------------------
# first run install.packages(".") to install package "." from CRAN if needed
library(MASS)
library(tibble)
library(dplyr)


# Set all of the fixed parameters ----------------------------------------------

# sample size
n <- 3500 

# Weibull shape and scale parameters
event_shape <- 1.6 
event_scale <- 122 
cens_shape <- 2.6 
cens_scale <- 8.2 

# log(HR) for each covariate 
beta_x1 <- log(0.8) 
beta_x2 <- log(1.05) 
beta_x3 <- log(1.25)
beta_x4 <- log(1.54)
beta_x5 <- log(1.18)
beta_x6 <- log(1.45)
beta_x7 <- log(1.10)
beta_x8 <- log(0.76)
beta_x9 <- log(0.64)
beta_x10 <- log(1.25)
beta_x11 <- log(0.48)

# proportion of missingness for each covariate
p_m_x1 <- 0.05
p_m_x3 <- 0.15
p_m_x4 <- 0.30

# log(OR) for joint missingness
gamma_11 <- 0
gamma_13 <- log(2.59)
gamma_14 <- log(0.31)

# log(OR) for covariates related to missingness
gamma_21 <- log(1.05)
gamma_23 <- log(0.8)
gamma_24 <- log(0.7)

# Covariate mean vector 
mymu <- c(0.6145, 57.6495, 2.3665, 0.4225, 0.1538, 0.2421, 0.4142, 0.8680, 
          0.1695, 0.1636, 0.8891)

# Covariate covariance matrix
mysigma <- matrix(
  c(0.2370, -1.3349, 0.0812, 0.0247, 0.0298, 0.0134, 0.0280, -0.0069, 
    0.0039, 0.0002, 0.0223, -1.3349, 196.2990, 0.6471, -0.6314, -0.0759, 
    0.3167, -0.7944, 0.0156, -0.1261, 0.0294, -0.8057, 0.0812, 0.6471, 
    1.0967, 0.0680, 0.0452, 0.0063, 0.0740, -0.0190, -0.0080, -0.0147, 
    -0.0075, 0.0247, -0.6314, 0.0680, 0.2441, 0.0063, -0.0089, 0.0468, 
    -0.0581, 0.0484, 0.0003, -0.0139, 0.0298, -0.0759, 0.0452, 0.0063, 
    0.1302, 0.0070, 0.0214, -0.0026, 0.0016, -0.0001, 0.0043, 0.0134, 
    0.3167, 0.0063, -0.0089, 0.0070, 0.0558, -0.0015, 0.0050, -0.0026, 
    0.0003, 0.0000, 0.0280, -0.7944, 0.0740, 0.0468, 0.0214, -0.0015, 
    0.2427, -0.0119, 0.0131, 0.0017, 0.0012, -0.0069, 0.0156, -0.0190, 
    -0.0581, -0.0026, 0.0050, -0.0119, 0.1146, -0.0318, 0.0001, 0.0086,
    0.0039, -0.1261, -0.0080, 0.0484, 0.0016, -0.0026, 0.0131, -0.0318, 
    0.1408, 0.0050, -0.0267, 0.0002, 0.0294, -0.0147, 0.0003, -0.0001, 
    0.0003, 0.0017, 0.0001, 0.0050, 0.1369, 0.0014, 0.0223, -0.8057, 
    -0.0075, -0.0139, 0.0043, 0.0000, 0.0012, 0.0086,  -0.0267, 0.0014, 0.0986),
  nrow = 11, ncol = 11, byrow = TRUE)


# Generate complete data ------------------------------------------------------

# set seed for replication
set.seed(20240815) 

# Generate covariates
covs <- 
  mvrnorm(n = n, mu = mymu, Sigma = mysigma) |> 
  as_tibble(.name_repair = ~ paste0("x", 1:11)) |> 
  mutate(
    x1 = ifelse(x1 <= 0.5, 0, 1),
    x4 = ifelse(x4 <= 0.5, 0, 1),
    x5 = ifelse(x5 <= 0.5, 0, 1),
    x7 = ifelse(x7 <= 0.5, 0, 1),
    x8 = ifelse(x8 <= 0.5, 0, 1),
    x9 = ifelse(x9 <= 0.5, 0, 1),
    x10 = ifelse(x10 <= 0.5, 0, 1),
    x11 = ifelse(x11 <= 0.5, 0, 1))

# Generate event and censoring times
u <- runif(n) 
t <- event_scale * 
  ((-log(u)) / 
     exp(beta_x1 * covs$x1 + 
           beta_x2 * covs$x2 +
           beta_x3 * covs$x3 + 
           beta_x4 * covs$x4 +
           beta_x5 * covs$x5 + 
           beta_x6 * covs$x6 + 
           beta_x7 * covs$x7 + 
           beta_x8 * covs$x8 + 
           beta_x9 * covs$x9 + 
           beta_x10 * covs$x10 + 
           beta_x11 * covs$x11))^(1/event_shape) 

c <- rweibull(n, cens_shape, cens_scale) # weibull censoring times

# Combine into single time/event indicator
s <- pmin(t, c)
delta <- ifelse(t <= c, 1, 0)

# Generate missingness
gamma_01 <- log(p_m_x1 / (1 - p_m_x1)) - gamma_21 * mean(covs$x2)
xg1 <- gamma_01 + gamma_21 * covs$x2
mp1 <- exp(xg1) / (1 + exp(xg1))
m1 <- rbinom(n, 1, prob = mp1) 
x1_miss = ifelse(m1 == 1, NA, covs$x1) 

gamma_03 <- log(p_m_x3 / (1 - p_m_x3)) - gamma_13 * p_m_x1 - 
  gamma_23 * mean(covs$x5)
xg3 <- gamma_03 + gamma_13 * ifelse(is.na(covs$x1), 1, 0) + gamma_23 * covs$x5
mp3 <- exp(xg3) / (1 + exp(xg3))
m3 <- rbinom(n, 1, prob = mp3) 
x3_miss = ifelse(m3 == 1, NA, covs$x3) 

gamma_04 <- log(p_m_x4 / (1 - p_m_x4)) - gamma_14 * p_m_x3 - 
  gamma_24 * mean(covs$x6)
xg4 <- gamma_04 + gamma_14 * ifelse(is.na(covs$x3), 1, 0) + gamma_24 * covs$x6
mp4 <- exp(xg4) / (1 + exp(xg4))
m4 <- rbinom(n, 1, prob = mp4) 
x4_miss = ifelse(m4 == 1, NA, covs$x4) 

# Complete and missing covariate data
dat0 <- tibble(
  s,
  delta,
  x1_miss,
  x3_miss,
  x4_miss
  ) |> 
  bind_cols(covs)


# Induce covariate missingness -------------------------------------------------

# Imputation models
imp_x1 <- glm(x1_miss ~ x2 + x5 + x6 + x7 + x8 + x9 + x10 + x11,
              data = dat0[!is.na(dat0$x1_miss), ],
              family = "binomial")

imp_x3 <- glm(x3_miss ~ x2 + x5 + x6 + x7 + x8 + x9 + x10 + x11,
              data = dat0[!is.na(dat0$x3_miss), ],
              family = "gaussian")

imp_x4 <- glm(x4_miss ~ x2 + x5 + x6 + x7 + x8 + x9 + x10 + x11,
              data = dat0[!is.na(dat0$x4_miss), ],
              family = "binomial")

# Predicted value
x1_pred <- ifelse(
  predict(imp_x1, newdata = dat0, type = "response") > 0.5, 1, 0)

x3_pred <- predict(imp_x3, newdata = dat0)

x4_pred <- ifelse(
  predict(imp_x4, newdata = dat0, type = "response") > 0.5, 1, 0)

# Combine observed and predicted
dat <-
  dat0 |>
  mutate(
    x1_imp = case_when(
      is.na(x1_miss) ~ x1_pred,
      TRUE ~ x1
    ),
    x3_imp = case_when(
      is.na(x3_miss) ~ x3_pred,
      TRUE ~ x3
    ),
    x4_imp = case_when(
      is.na(x4_miss) ~ x4_pred,
      TRUE ~ x4
    )
  )