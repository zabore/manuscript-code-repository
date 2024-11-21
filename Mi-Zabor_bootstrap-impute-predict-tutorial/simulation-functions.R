#############################################################################
# Functions for data simulation: 
# 0. Supporting functions - missing_or, cloglog, invcloglog
# 1. Function to simulate samples using mean and covariance data
# 2. Function to generate missing data
# 3. Function to impute missing data
# 4. Function to calculate performance and prediction for full, cc, and bi
# 5. Function to summarize bias of performance and prediction results
# 6. Function to check imputation accuracy
# 7. Function to check imputation success
# 99. Function to do internal validation for a single dataset
#############################################################################

# 0. Supporting functions ---------------------------------------------------

# Function to calculate OR of missingness
# Example: missing_or(x1_missing = 0.05, x2_missing = 0.05, 
#                     joint_missing = 0.015)
missing_or <- function(x1_missing, x2_missing, joint_missing) {
  
  x1_missing_only = x1_missing - joint_missing
  x2_missing_only = x2_missing - joint_missing
  no_missing = 1 - x1_missing - x2_missing_only
  
  missing_or = (no_missing * joint_missing) /
    (x1_missing_only * x2_missing_only)
  
  return(missing_or)
  
}

# Functions to do the complementary-log-log transformation and 
# inverse-complementary-log-log transformation
cloglog <- function(x) log(-log(1 - x))
inv_cloglog <- function(x) 1 - exp(-exp(x))


# 1. Function to simulate samples using mean and covariance data ------------
# Example: sim_sample(sample_size = 750)
sim_sample <- function(sample_size) {
  
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
  
  covs <- rmvnorm(n = sample_size, mean = mymu, sigma = mysigma) |>
    as_tibble(.name_repair = ~ paste0("x", 1:11)) |> 
    mutate(
      x1 = ifelse(x1 <= 0.5, 0, 1),
      x4 = ifelse(x4 <= 0.5, 0, 1),
      x5 = ifelse(x5 <= 0.5, 0, 1),
      x7 = ifelse(x7 <= 0.5, 0, 1),
      x8 = ifelse(x8 <= 0.5, 0, 1),
      x9 = ifelse(x9 <= 0.5, 0, 1),
      x10 = ifelse(x10 <= 0.5, 0, 1),
      x11 = ifelse(x11 <= 0.5, 0, 1)
    )
  
  n <- sample_size # sample size
  event_shape <- 1.6 # Weibull shape no PMRT
  event_scale <- 122 # Weibull scale no PMRT
  cens_shape <- 2.6 # Weibull shape censor
  cens_scale <- 8.2 # Weibull scale censor
  beta_x1 <- log(0.8) # log HR for PMRT
  beta_x2 <- log(1.05) # log HR for age
  beta_x3 <- log(1.25) # log HR for tumor size
  beta_x4 <- log(1.54) # log HR for grade
  beta_x5 <- log(1.18) # log HR for positive ln
  beta_x6 <- log(1.45) # log HR for nodal positive ratio
  beta_x7 <- log(1.10) # log HR for lymph vascular space invasion
  beta_x8 <- log(0.76) # log HR for er or pr status
  beta_x9 <- log(0.64) # log HR for HER2
  beta_x10 <- log(1.25) # log HR for tumor quadrant
  beta_x11 <- log(0.48) # log HR for optimal systemic therapy
  
  # Generate event and censoring times
  u <- runif(n) # uniform random variable
  t <- event_scale *
    ((-log(u)) / exp(beta_x1 * covs$x1 + beta_x2 * covs$x2 +
                       beta_x3 * covs$x3 + beta_x4 * covs$x4 + 
                       beta_x5 * covs$x5 + beta_x6 * covs$x6 +
                       beta_x7 * covs$x7 + beta_x8 * covs$x8 + 
                       beta_x9 * covs$x9 + beta_x10 * covs$x10 + 
                       beta_x11 * covs$x11))^(1/event_shape) # weibull event times
  c <- rweibull(n, cens_shape, cens_scale) # weibull censoring times
  
  # Combine into single time/event indicator
  os_yrs <- pmin(t, c)
  os_event <- ifelse(t <= c, 1, 0)
  
  sample <- tibble(
    true_event_time = t,
    os_yrs = os_yrs,
    os_event = os_event) |> 
      bind_cols(covs) |> 
      mutate(id = row_number())
  
  return(sample)
  
}


# 2. Function to generate missing data -----------------------------------
#
# Nine missing patterns:
# 1. x1 missing in 5%
# 2. x1 missing in 15%
# 3. x1 missing in 60%
# 4. x1, x3, x4 missing in 5%, 5%, 5%
# 5. x1, x3, x4 missing in 5%, 15%, 30%
# 6. x1, x3, x4 missing in 15%, 30%, 60%
# 7. x1, x3, x4, x7, x10, x11 missing in 5%, 5%, 5%, 5%, 5%, 5%
# 8. x1, x3, x4, x7, x10, x11 missing in 5%, 5%, 15%, 15%, 30%, 30%
# 9. x1, x3, x4, x7, x10, x11 missing in 15%, 15%, 30%, 30%, 60%, 60%
#
# Example: missing_dat(dataset = dat_750, missing_var = c('x1','x3','x4'), 
#                      missing_pct = c(0.05, 0.15, 0.3))
missing_dat <- function(dataset, missing_var, missing_pct) {
  
  # log(OR) for joint missingness
  gamma_0505 <- log(11.2) # based on 1.5% missing both values
  gamma_0515 <- log(4.21) # based on 2% missing both values
  gamma_1515 <- log(3.75) # based on 5% missing both values
  gamma_1530 <- log(2.78) # based on 7.5% missing both values
  gamma_3030 <- log(1.25) # based on 10% missing both values
  gamma_3060 <- log(1.50) # based on 20% missing both values
  gamma_6060 <- log(2.00) # based on 40% missing both values
  
  # log(OR) for covariates related to missingness
  gamma_x2 <- log(1.05)
  gamma_x5 <- log(0.8)
  gamma_x6 <- log(0.7)
  gamma_x8 <- log(0.9)
  gamma_x9 <- log(0.6)
  
  data_missing <- dataset
  n <- nrow(data_missing)
  
  # generate missing variables
  gamma_01 <- log(missing_pct[1] / (1 - missing_pct[1])) - 
    gamma_x2 * mean(data_missing$x2)
  xg1 <- gamma_01 + gamma_x2 * data_missing$x2
  mp1 <- exp(xg1) / (1 + exp(xg1))
  
  data_missing$m1 <- rbinom(n, 1, prob = mp1) # randomly generate missing indicator
  data_missing <- data_missing %>% mutate(x1_miss = ifelse(m1 == 1, NA, x1))
  
  if (length(missing_var) > 1) {
    
    if (length(missing_var) == 3) {
      
      if (all(missing_pct == c(0.05, 0.05, 0.05))) {
        
        gamma_03 <- log(missing_pct[2] / (1 - missing_pct[2])) - 
          gamma_0505 * missing_pct[1] - gamma_x5 * mean(data_missing$x5)
        xg3 <- gamma_03 + gamma_0505 * ifelse(is.na(data_missing$x1_miss), 1, 0) + 
          gamma_x5 * data_missing$x5
        mp3 <- exp(xg3) / (1 + exp(xg3))
        
        data_missing$m3 <- rbinom(n, 1, prob = mp3)
        data_missing <- data_missing %>% 
          mutate(x3_miss = ifelse(m3 == 1, NA, x3))
        
        gamma_04 <- log(missing_pct[3] / (1 - missing_pct[3])) - 
          gamma_0505 * missing_pct[2] - gamma_x6 * mean(data_missing$x6)
        xg4 <- gamma_04 + gamma_0505 * ifelse(is.na(data_missing$x3_miss), 1, 0) + 
          gamma_x6 * data_missing$x6
        mp4 <- exp(xg4) / (1 + exp(xg4))
        
        data_missing$m4 <- rbinom(n, 1, prob = mp4)
        data_missing <- data_missing %>% 
          mutate(x4_miss = ifelse(m4 == 1, NA, x4))
        
      }
      
      if (all(missing_pct == c(0.05, 0.15, 0.3))) {
        
        gamma_03 <- log(missing_pct[2] / (1 - missing_pct[2])) - 
          gamma_0515 * missing_pct[1] - gamma_x5 * mean(data_missing$x5)
        xg3 <- gamma_03 + gamma_0515 * ifelse(is.na(data_missing$x1_miss), 1, 0) + 
          gamma_x5 * data_missing$x5
        mp3 <- exp(xg3) / (1 + exp(xg3))
        
        data_missing$m3 <- rbinom(n, 1, prob = mp3)
        data_missing <- data_missing %>% 
          mutate(x3_miss = ifelse(m3 == 1, NA, x3))
        
        gamma_04 <- log(missing_pct[3] / (1 - missing_pct[3])) - 
          gamma_1530 * missing_pct[2] - gamma_x6 * mean(data_missing$x6)
        xg4 <- gamma_04 + gamma_1530 * ifelse(is.na(data_missing$x3_miss), 1, 0) + 
          gamma_x6 * data_missing$x6
        mp4 <- exp(xg4) / (1 + exp(xg4))
        
        data_missing$m4 <- rbinom(n, 1, prob = mp4)
        data_missing <- data_missing %>% 
          mutate(x4_miss = ifelse(m4 == 1, NA, x4))
        
      }
      
      if (all(missing_pct == c(0.15, 0.3, 0.6))) {
        
        gamma_03 <- log(missing_pct[2] / (1 - missing_pct[2])) - 
          gamma_1530 * missing_pct[1] - gamma_x5 * mean(data_missing$x5)
        xg3 <- gamma_03 + gamma_1530 * ifelse(is.na(data_missing$x1_miss), 1, 0) + 
          gamma_x5 * data_missing$x5
        mp3 <- exp(xg3) / (1 + exp(xg3))
        
        data_missing$m3 <- rbinom(n, 1, prob = mp3)
        data_missing <- data_missing %>% 
          mutate(x3_miss = ifelse(m3 == 1, NA, x3))
        
        gamma_04 <- log(missing_pct[3] / (1 - missing_pct[3])) - 
          gamma_3060 * missing_pct[2] - gamma_x6 * mean(data_missing$x6)
        xg4 <- gamma_04 + gamma_3060 * ifelse(is.na(data_missing$x3_miss), 1, 0) + 
          gamma_x6 * data_missing$x6
        mp4 <- exp(xg4) / (1 + exp(xg4))
        
        data_missing$m4 <- rbinom(n, 1, prob = mp4)
        data_missing <- data_missing %>% 
          mutate(x4_miss = ifelse(m4 == 1, NA, x4))
        
      }      
      
    }
    
    if (length(missing_var) == 6) {
      
      if (all(missing_pct == c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05))) {
        
        gamma_03 <- log(missing_pct[2] / (1 - missing_pct[2])) - 
          gamma_0505 * missing_pct[1] - gamma_x5 * mean(data_missing$x5)
        xg3 <- gamma_03 + gamma_0505 * ifelse(is.na(data_missing$x1_miss), 1, 0) + 
          gamma_x5 * data_missing$x5
        mp3 <- exp(xg3) / (1 + exp(xg3))
        
        data_missing$m3 <- rbinom(n, 1, prob = mp3)
        data_missing <- data_missing %>% 
          mutate(x3_miss = ifelse(m3 == 1, NA, x3))
        
        gamma_04 <- log(missing_pct[3] / (1 - missing_pct[3])) - 
          gamma_0505 * missing_pct[2] - gamma_x6 * mean(data_missing$x6)
        xg4 <- gamma_04 + gamma_0505 * ifelse(is.na(data_missing$x3_miss), 1, 0) + 
          gamma_x6 * data_missing$x6
        mp4 <- exp(xg4) / (1 + exp(xg4))
        
        data_missing$m4 <- rbinom(n, 1, prob = mp4)
        data_missing <- data_missing %>% 
          mutate(x4_miss = ifelse(m4 == 1, NA, x4))
        
        gamma_07 <- log(missing_pct[4] / (1 - missing_pct[4])) - 
          gamma_0505 * missing_pct[3] - gamma_x8 * mean(data_missing$x8)
        xg7 <- gamma_07 + gamma_0505 * ifelse(is.na(data_missing$x4_miss), 1, 0) + 
          gamma_x8 * data_missing$x8
        mp7 <- exp(xg7) / (1 + exp(xg7))
        
        data_missing$m7 <- rbinom(n, 1, prob = mp7)
        data_missing <- data_missing %>% 
          mutate(x7_miss = ifelse(m7 == 1, NA, x7))
        
        gamma_010 <- log(missing_pct[5] / (1 - missing_pct[5])) - 
          gamma_0505 * missing_pct[4] - gamma_x9 * mean(data_missing$x9)
        xg10 <- gamma_010 + gamma_0505 * ifelse(is.na(data_missing$x7_miss), 1, 0) + 
          gamma_x9 * data_missing$x9
        mp10 <- exp(xg10) / (1 + exp(xg10))
        
        data_missing$m10 <- rbinom(n, 1, prob = mp10)
        data_missing <- data_missing %>% 
          mutate(x10_miss = ifelse(m10 == 1, NA, x10))
        
        gamma_011 <- log(missing_pct[6] / (1 - missing_pct[6])) - 
          gamma_0505 * 
          missing_pct[5] - gamma_x2 * mean(data_missing$x2)
        xg11 <- gamma_011 + gamma_0505 * ifelse(is.na(data_missing$x10_miss), 1, 0) + 
          gamma_x2 * data_missing$x2
        mp11 <- exp(xg11) / (1 + exp(xg11))
        
        data_missing$m11 <- rbinom(n, 1, prob = mp11)
        data_missing <- data_missing %>% 
          mutate(x11_miss = ifelse(m11 == 1, NA, x11))
        
      }
      
      if (all(missing_pct == c(0.05, 0.05, 0.15, 0.15, 0.3, 0.3))) {
        
        gamma_03 <- log(missing_pct[2] / (1 - missing_pct[2])) - 
          gamma_0505 * missing_pct[1] - gamma_x5 * mean(data_missing$x5)
        xg3 <- gamma_03 + gamma_0505 * ifelse(is.na(data_missing$x1_miss), 1, 0) + 
          gamma_x5 * data_missing$x5
        mp3 <- exp(xg3) / (1 + exp(xg3))
        
        data_missing$m3 <- rbinom(n, 1, prob = mp3)
        data_missing <- data_missing %>% 
          mutate(x3_miss = ifelse(m3 == 1, NA, x3))
        
        gamma_04 <- log(missing_pct[3] / (1 - missing_pct[3])) - 
          gamma_0515 * 
          missing_pct[2] - gamma_x6 * mean(data_missing$x6)
        xg4 <- gamma_04 + gamma_0515 * ifelse(is.na(data_missing$x3_miss), 1, 0) + 
          gamma_x6 * data_missing$x6
        mp4 <- exp(xg4) / (1 + exp(xg4))
        
        data_missing$m4 <- rbinom(n, 1, prob = mp4)
        data_missing <- data_missing %>% 
          mutate(x4_miss = ifelse(m4 == 1, NA, x4))
        
        gamma_07 <- log(missing_pct[4] / (1 - missing_pct[4])) - 
          gamma_1515 * missing_pct[3] - gamma_x8 * mean(data_missing$x8)
        xg7 <- gamma_07 + gamma_1515 * ifelse(is.na(data_missing$x4_miss), 1, 0) + 
          gamma_x8 * data_missing$x8
        mp7 <- exp(xg7) / (1 + exp(xg7))
        
        data_missing$m7 <- rbinom(n, 1, prob = mp7)
        data_missing <- data_missing %>% 
          mutate(x7_miss = ifelse(m7 == 1, NA, x7))
        
        gamma_010 <- log(missing_pct[5] / (1 - missing_pct[5])) - 
          gamma_1530 * missing_pct[4] - gamma_x9 * mean(data_missing$x9)
        xg10 <- gamma_010 + gamma_1530 * ifelse(is.na(data_missing$x7_miss), 1, 0) + 
          gamma_x9 * data_missing$x9
        mp10 <- exp(xg10) / (1 + exp(xg10))
        
        data_missing$m10 <- rbinom(n, 1, prob = mp10)
        data_missing <- data_missing %>% 
          mutate(x10_miss = ifelse(m10 == 1, NA, x10))
        
        gamma_011 <- log(missing_pct[6] / (1 - missing_pct[6])) - 
          gamma_3030 * missing_pct[5] - gamma_x2 * mean(data_missing$x2)
        xg11 <- gamma_011 + gamma_3030 * ifelse(is.na(data_missing$x10_miss), 1, 0) + 
          gamma_x2 * data_missing$x2
        mp11 <- exp(xg11) / (1 + exp(xg11))
        
        data_missing$m11 <- rbinom(n, 1, prob = mp11)
        data_missing <- data_missing %>% 
          mutate(x11_miss = ifelse(m11 == 1, NA, x11)) 
        
      }
      
      if (all(missing_pct == c(0.15, 0.15, 0.3, 0.3, 0.6, 0.6))) {
        
        gamma_03 <- log(missing_pct[2] / (1 - missing_pct[2])) - 
          gamma_1515 * missing_pct[1] - gamma_x5 * mean(data_missing$x5)
        xg3 <- gamma_03 + gamma_1515 * ifelse(is.na(data_missing$x1_miss), 1, 0) + 
          gamma_x5 * data_missing$x5
        mp3 <- exp(xg3) / (1 + exp(xg3))
        
        data_missing$m3 <- rbinom(n, 1, prob = mp3)
        data_missing <- data_missing %>% 
          mutate(x3_miss = ifelse(m3 == 1, NA, x3))
        
        gamma_04 <- log(missing_pct[3] / (1 - missing_pct[3])) - 
          gamma_1530 * missing_pct[2] - gamma_x6 * mean(data_missing$x6)
        xg4 <- gamma_04 + gamma_1530 * ifelse(is.na(data_missing$x3_miss), 1, 0) + 
          gamma_x6 * data_missing$x6
        mp4 <- exp(xg4) / (1 + exp(xg4))
        
        data_missing$m4 <- rbinom(n, 1, prob = mp4)
        data_missing <- data_missing %>% 
          mutate(x4_miss = ifelse(m4 == 1, NA, x4))
        
        gamma_07 <- log(missing_pct[4] / (1 - missing_pct[4])) - 
          gamma_3030 * missing_pct[3] - gamma_x8 * mean(data_missing$x8)
        xg7 <- gamma_07 + gamma_3030 * ifelse(is.na(data_missing$x4_miss), 1, 0) + 
          gamma_x8 * data_missing$x8
        mp7 <- exp(xg7) / (1 + exp(xg7))
        
        data_missing$m7 <- rbinom(n, 1, prob = mp7)
        data_missing <- data_missing %>% 
          mutate(x7_miss = ifelse(m7 == 1, NA, x7))
        
        gamma_010 <- log(missing_pct[5] / (1 - missing_pct[5])) - 
          gamma_3060 * missing_pct[4] - gamma_x9 * mean(data_missing$x9)
        xg10 <- gamma_010 + gamma_3060 * ifelse(is.na(data_missing$x7_miss), 1, 0) + 
          gamma_x9 * data_missing$x9
        mp10 <- exp(xg10) / (1 + exp(xg10))
        
        data_missing$m10 <- rbinom(n, 1, prob = mp10)
        data_missing <- data_missing %>% 
          mutate(x10_miss = ifelse(m10 == 1, NA, x10))
        
        gamma_011 <- log(missing_pct[6] / (1 - missing_pct[6])) - 
          gamma_6060 * missing_pct[5] - gamma_x2 * mean(data_missing$x2)
        xg11 <- gamma_011 + gamma_6060 * ifelse(is.na(data_missing$x10_miss), 1, 0) + 
          gamma_x2 * data_missing$x2
        mp11 <- exp(xg11) / (1 + exp(xg11))
        
        data_missing$m11 <- rbinom(n, 1, prob = mp11)
        data_missing <- data_missing %>% 
          mutate(x11_miss = ifelse(m11 == 1, NA, x11))
        
      }      
      
    }        
    
  }    
  
  return(data_missing)
  
}


# 3. Function to impute missing data ---------------------------------------------
#
# Three imputation methods:
# 1 - imupte all
# 2 - only impute covs missing > 10%,
# 3 - only impute subjects missing 1-2 covs
#
# Example: imp_dat(dataset = data_missing, impute_method = 1)
imp_dat <- function(dataset, impute_method) {
  
  data_impute <- dataset
  
  # summarize missing N by row
  
  data_impute$n_miss <- apply(data_impute, 1, function(x) sum(is.na(x)))
  
  # summarize missing covs and missing percent
  
  missing_col <- colnames(data_impute)[colSums(is.na(data_impute)) > 0]
  
  missing_var <- sub("_.*", "", missing_col)
  
  missing_pct <- as.numeric(colSums(is.na(data_impute))[missing_col] / 
                              nrow(data_impute))
  
  complete_var <- colnames(data_impute)[
    names(data_impute) %in% paste("x", seq(1, 11), sep = "") & 
      !(names(data_impute) %in% missing_var)]
  
  for (i in 1:length(missing_var)) {
    
    x_ite <- missing_var[i]
    
    x_binary <- ifelse(length(unique(data_impute[[x_ite]])) == 2, T, F)
    
    mod_family <- ifelse(x_binary,'binomial','gaussian')
    
    if(!(impute_method == 2 & missing_pct[i] < 0.1)) {
      
      imp_mod <- tryCatch({
        glm(as.formula(
          paste0(x_ite, " ~ ", paste0(complete_var, collapse = " + "))), 
          data = data_impute[!is.na(data_impute[[paste0(x_ite,'_miss')]]), ], 
          family = mod_family)}, 
        warning = function(war){
          NA}, 
        error = function(err){
          NA
        })
      
      if(x_binary){
        
        if(any(class(imp_mod) == "logical")) {
          x_ite_pred <- NA
          } else {
          x_ite_pred <- ifelse(
            predict(imp_mod, newdata = data_impute, type = "response") > 0.5, 
            1, 0)
          }
        
      }else{
        
        if(any(class(imp_mod) == "logical")) {
          x_ite_pred <- NA
          } else {
          x_ite_pred <- predict(imp_mod, newdata = data_impute)
          }
        
      }
      
      if(impute_method == 3) {
        
        # 3 - only impute subjects missing 1-2 covs
        if(any(class(imp_mod) == "logical")){
          data_impute[[paste0(x_ite,'_pred')]] <- rep(NA, nrow(data_impute))
        } else {
          data_impute[[paste0(x_ite,'_pred')]] <- 
            ifelse((is.na(data_impute[[paste0(x_ite,'_miss')]]) & 
                      data_impute$n_miss < 3),
                   x_ite_pred, 
                   data_impute[[paste0(x_ite,'_miss')]])
        }
        
      }else{
        
        if(any(class(imp_mod) == "logical")){
          data_impute[[paste0(x_ite,'_pred')]] <- rep(NA, nrow(data_impute))
        } else {
        
        data_impute[[paste0(x_ite,'_pred')]] <- 
          ifelse(is.na(data_impute[[paste0(x_ite,'_miss')]]),
                 x_ite_pred, 
                 data_impute[[paste0(x_ite,'_miss')]])
        }
      }
      
    }else{
      
      data_impute[[paste0(x_ite,'_pred')]] <- 
        data_impute[[paste0(x_ite,'_miss')]]
      
    }
  }
  
  return(data_impute)
  
}



# 4. Function to calculate performance and prediction by full, cc, and bi -------
#
# imp_dat: single imputed dataset
# boot_imp_dat: B bootstrap samples of imp_dat 
#
# Example: calc_perf_pred(impute_list[[49]])
calc_perf_pred <- function(imp_dat, boot_imp_dat, predtime = 5) {
  
  # Fix "no information performance" for each metric
  gamma_auc <- 0.5
  gamma_brier <- 0.25
  
  # identify missing cols, complete cols, and predicted cols

  missing_col <- colnames(imp_dat)[colSums(is.na(imp_dat)) > 0 & 
                                     grepl("_miss$", colnames(imp_dat))]
  
  missing_var <- sub("_.*", "", missing_col)
  
  pred_col <- gsub("_miss", "_pred", missing_col)
  
  missing_pct <- as.numeric(colSums(is.na(imp_dat))[missing_col] / nrow(imp_dat))
  
  complete_var <- colnames(imp_dat)[
    names(imp_dat) %in% paste("x", seq(1, 11), sep = "") & 
      !(names(imp_dat) %in% missing_var)]
  
  # Full data
  
  mod_full <- tryCatch({
    coxph(Surv(os_yrs, os_event) ~ x1 + x2 + x3 + x4 + x5 + x6 + 
            x7 + x8 + x9 + x10 + x11, 
          data = imp_dat, x = TRUE, y = TRUE)},
    warning = function(war){
      NA}, 
    error = function(err){
      NA
    })
  
  mod_full_success <- ifelse(isTRUE(is.na(mod_full)[[1]]), 0, 1)
  
  # Calculate performance in original data - full data
  app_auc_brier_full <- 
    tryCatch({
      Score(list("fit" = mod_full),
            formula = as.formula("Surv(os_yrs, os_event) ~ x1 + x2 + x3 + x4 + 
            x5 + x6 + x7 + x8 + x9 + x10 + x11"),
            data = imp_dat,
            se.fit = FALSE,
            conf.int = FALSE,
            times = predtime)},
      error = function(err){
        NA
      })
  
  app_auc_full <- 
    tryCatch({
      app_auc_brier_full[["AUC"]][["score"]][["AUC"]]},
      error = function(err){
        NA
      })
  
  app_brier_full <- tryCatch({
    app_auc_brier_full[["Brier"]][["score"]][["Brier"]][[2]]},
    error = function(err){
      NA
    })
  
  # Fit a bootstrap model to each bootstrap sample - full data
  boot_mod_full <- 
    map(
      boot_imp_dat,
      function(z){
        tryCatch({
          coxph(as.formula("Surv(os_yrs, os_event) ~ x1 + x2 + x3 + x4 + 
                           x5 + x6 + x7 + x8 + x9 + x10 + x11"), data = z, 
                x = TRUE, y = TRUE)},
          warning = function(war){
            NA}, 
          error = function(err){
            NA
          })
      }
     ) 
  
  
  boot_mod_full_success = as.numeric(sum(!is.na(boot_mod_full)))
   
   # Bootstrap AUC and Brier
   b_auc_brier_full <-
     map2(
       boot_imp_dat,
       boot_mod_full,
       function(x, y) {
         tryCatch({
           Score(list("b_fit" = y),
                 formula = as.formula("Surv(os_yrs, os_event) ~ x1 + x2 + x3 + 
                 x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11"),
                 data = x,
                 se.fit = FALSE,
                 conf.int = FALSE,
                 times = predtime)},
           error = function(err){
             NA
           })
       }
       )
   
   boot_auc_brier_full_success = as.numeric(sum(!is.na(b_auc_brier_full)))
   
   b_auc_full <-
     map_dbl(
       b_auc_brier_full,
       function(x){
         tryCatch({
           x[["AUC"]][["score"]][["AUC"]]},
           error = function(err){
             NA
           })
       }
     )
   
   b_brier_full <-
     map_dbl(
       b_auc_brier_full,
       function(x){
         tryCatch({
           x[["Brier"]][["score"]][["Brier"]][[2]]},
           error = function(err){
             NA
           })
       }
     )
   
   # Calculate performance in original data
   o_auc_brier_full <-
     map(
       boot_mod_full,
       function(x) {
         tryCatch({
           Score(list("ofit" = x),
                 formula = as.formula("Surv(os_yrs, os_event) ~ x1 + x2 + x3 + 
                 x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11"),
                 data = imp_dat,
                 se.fit = FALSE,
                 conf.int = FALSE,
                 times = predtime)},
           error = function(err){
             NA
           })
       }
     )
   
   o_auc_full <-
     map_dbl(
       o_auc_brier_full,
       function(x){
         tryCatch({
           x[["AUC"]][["score"]][["AUC"]]},
           error = function(err){
             NA
           })
       }
     )
   
   o_brier_full <-
     map_dbl(
       o_auc_brier_full,
       function(x){
         tryCatch({
           x[["Brier"]][["score"]][["Brier"]][[2]]},
           error = function(err){
             NA
           })
       }
     )
   
   # Calculate bootstrap optimism corrected AUC and Brier
   boot_auc_full <- app_auc_full - mean(b_auc_full - o_auc_full, na.rm = T)
   boot_brier_full <- app_brier_full - 
     mean(b_brier_full - o_brier_full, na.rm = T)
   
   # .632 performance
   
   # Obtain test datasets for each bootstrap sample
   test_dat_full <- 
     boot_imp_dat |> 
     map(
       function(x){
         if(nrow(x) == 0) {
           x
         } else{
           imp_dat |> filter(!(id %in% unique(x[["id"]])))
         }
       }
     )
   
   # Calculate performance in each test dataset
   test_auc_brier_full <-
     map2(
       test_dat_full,
       boot_mod_full,
       function(x, y){
         tryCatch({
           Score(list("b_fit" = y),
                 formula = as.formula("Surv(os_yrs, os_event) ~ x1 + x2 + x3 + 
                 x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11"),
                 data = x,
                 se.fit = FALSE,
                 conf.int = FALSE,
                 times = predtime)},
           error = function(err){
             NA
           })
       }
     )
   
   test_auc_full <-
     map_dbl(
       test_auc_brier_full,
       function(x){
         tryCatch({
           x[["AUC"]][["score"]][["AUC"]]},
           error = function(err){
             NA
           })
       }
     )
   
   test_brier_full <-
     map_dbl(
       test_auc_brier_full,
       function(x){
         tryCatch({
           x[["Brier"]][["score"]][["Brier"]][[2]]},
           error = function(err){
             NA
           })
       }
     )
   
   # Calculate .632 AUC and Brier
   boot632_auc_full <- .368 * app_auc_full + 
     .632 * mean(test_auc_full, na.rm = T)
   boot632_brier_full <- .368 * app_brier_full + 
     .632 * mean(test_brier_full, na.rm = T)
   
   # .632+ performance
   
   # Define relative overfitting rate
   r_auc_full <- (mean(test_auc_full, na.rm = T) - app_auc_full) / 
     (gamma_auc - app_auc_full)
   r_brier_full <- (mean(test_brier_full, na.rm = T) - app_brier_full) / 
     (gamma_brier - app_brier_full)
   
   # Define weights
   w_auc_full <- .632 / (1 - .368 * r_auc_full)
   w_brier_full <- .632 / (1 - .368 * r_brier_full)
   
   # Calculate the .632+ AUC and Brier
   boot632plus_auc_full <- (1 - w_auc_full) * app_auc_full + 
     w_auc_full * mean(test_auc_full, na.rm = T)
   boot632plus_brier_full <- (1 - w_brier_full) * app_brier_full + 
     w_brier_full * mean(test_brier_full, na.rm = T)
   
   # Get predictions from the full model across the bootstrap samples
   preds_full <- 
     boot_imp_dat |> 
     map(
       ~predict(mod_full, 
                newdata = .x |> mutate(time = predtime), 
                type = "survival")
     )
   
   preds_avg_full <- inv_cloglog(rowMeans(cloglog(
     as.data.frame(do.call(cbind, preds_full))
   ), na.rm = T)) 
   
  
  # Complete case data
  
  # Make a complete case dataset
  imp_dat_cc <- 
    imp_dat |> 
    filter(if_all(all_of(missing_col), ~ !is.na(.)))
  
  # Make complete case bootstrap datasets
  boot_imp_dat_cc <- 
    boot_imp_dat |> 
    map(
      ~ .x |> filter(if_all(all_of(missing_col), ~ !is.na(.)))
    )

  # We have two different models here 
  # mod_cc0 is on the full dataset that includes missing values, and will be used to get predictions with missing values
  mod_cc0 <- tryCatch({
    coxph(as.formula(paste0("Surv(os_yrs, os_event) ~ ", 
                            paste0(complete_var, collapse = " + "), "+", 
                            paste0(missing_col, collapse = " + "))),
          data = imp_dat, x = TRUE, y = TRUE)}, 
    warning = function(war){
      NA
    }, 
    error = function(err){
      NA
    })
  
  # check if mod_cc0 is success - return 0 if failed
  mod_cc0_success <- ifelse(isTRUE(is.na(mod_cc0)[[1]]), 0, 1)
  
  # Fit the complete case model on the original data - complete cases only
  # Need this mod_cc model for Score to work below - it doesn't play well with NAs
  mod_cc <- tryCatch({
    coxph(Surv(os_yrs, os_event) ~ x1 + x2 + x3 + x4 + x5 + x6 + 
            x7 + x8 + x9 + x10 + x11, 
          data = imp_dat_cc, x = TRUE, y = TRUE)}, 
    warning = function(war){
      NA
    }, 
    error = function(err){
      NA
    })
  
  # Get AUC and Brier using the complete case model on the original data
  app_auc_brier_cc <- 
    tryCatch({
      Score(list("fit" = mod_cc),
            formula = as.formula("Surv(os_yrs, os_event) ~ x1 + x2 + x3 +
            x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11"),
            data = imp_dat_cc,
            se.fit = FALSE,
            conf.int = FALSE,
            times = predtime)},
      error = function(err){
        NA
      })
  
  app_auc_cc <- 
    tryCatch({
      app_auc_brier_cc[["AUC"]][["score"]][["AUC"]]},
      error = function(err){
        NA
      })
  app_brier_cc <- tryCatch({
    app_auc_brier_cc[["Brier"]][["score"]][["Brier"]][[2]]},
    error = function(err){
      NA
    })
  
  # Fit a bootstrap model to each bootstrap sample
  boot_mod_cc <- 
    map(
      boot_imp_dat_cc,
      function(z){
        tryCatch({
          coxph(as.formula("Surv(os_yrs, os_event) ~ x1 + x2 + x3 +
                x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11"), data = z, 
                x = TRUE, y = TRUE)},
          warning = function(war){
            NA}, 
          error = function(err){
            NA
          })
      }
    ) 
  
  boot_mod_cc_success = as.numeric(sum(!is.na(boot_mod_cc)))
  
  # Get AUC and Brier using each bootstrap model on the bootstrap sample
  boot_auc_brier_cc <-
    map2(
      boot_imp_dat_cc,
      boot_mod_cc,
      function(x, y) {
        tryCatch({
            Score(list("b_fit" = y),
                formula = as.formula("Surv(os_yrs, os_event) ~ x1 + x2 + x3 +
                x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11"),
                data = x,
                se.fit = FALSE,
                conf.int = FALSE,
                times = predtime)
          },
          error = function(err){
            NA
          })
      }
    )
  
  boot_auc_brier_cc_success = as.numeric(sum(!is.na(boot_auc_brier_cc)))
  
  b_auc_cc <-
    map_dbl(
      boot_auc_brier_cc,
      function(x){
        tryCatch({
          x[["AUC"]][["score"]][["AUC"]]},
          error = function(err){
            NA
          })
      }
    )
  
  b_brier_cc <-
    map_dbl(
      boot_auc_brier_cc,
      function(x){
        tryCatch({
          x[["Brier"]][["score"]][["Brier"]][[2]]},
          error = function(err){
            NA
          })
      }
    )
  
  # Get AUC and Brier using each bootstrap model on the original data
  o_auc_brier_cc <-
    map(
      boot_mod_cc,
      function(x) {
        tryCatch({
          Score(list("ofit" = x),
                formula = as.formula("Surv(os_yrs, os_event) ~ x1 + x2 + x3 +
                x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11"),
                data = imp_dat_cc,
                se.fit = FALSE,
                conf.int = FALSE,
                times = predtime)},
          error = function(err){
            NA
          })
      }
    )
  
  o_auc_cc <-
    map_dbl(
      o_auc_brier_cc,
      function(x){
        tryCatch({
          x[["AUC"]][["score"]][["AUC"]]},
          error = function(err){
            NA
          })
      }
    )
  
  o_brier_cc <-
    map_dbl(
      o_auc_brier_cc,
      function(x){
        tryCatch({
          x[["Brier"]][["score"]][["Brier"]][[2]]},
          error = function(err){
            NA
          })
      }
    )
  
  # Calculate bootstrap optimism corrected AUC and Brier
  boot_auc_cc <- app_auc_cc - mean(b_auc_cc - o_auc_cc, na.rm = T)
  boot_brier_cc <- app_brier_cc - mean(b_brier_cc - o_brier_cc, na.rm = T)
  
  
  # .632 performance
  
  # Obtain test datasets for each bootstrap sample
  test_dat_cc <- 
    boot_imp_dat_cc |> 
    map(
      function(x){
        if(nrow(x) == 0) {
          x
        } else{
          imp_dat_cc |> filter(!(id %in% unique(x[["id"]])))
        }
      }
    )
  
  # Calculate performance in each test dataset
  test_auc_brier_cc <-
    map2(
      test_dat_cc,
      boot_mod_cc,
      function(x, y){
        tryCatch({
          Score(list("b_fit" = y),
                formula = as.formula("Surv(os_yrs, os_event) ~ x1 + x2 + x3 + 
                 x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11"),
                data = x,
                se.fit = FALSE,
                conf.int = FALSE,
                times = predtime)},
          error = function(err){
            NA
          })
      }
    )
  
  test_auc_cc <-
    map_dbl(
      test_auc_brier_cc,
      function(x){
        tryCatch({
          x[["AUC"]][["score"]][["AUC"]]},
          error = function(err){
            NA
          })
      }
    )
  
  test_brier_cc <-
    map_dbl(
      test_auc_brier_cc,
      function(x){
        tryCatch({
          x[["Brier"]][["score"]][["Brier"]][[2]]},
          error = function(err){
            NA
          })
      }
    )
  
  # Calculate .632 AUC and Brier
  boot632_auc_cc <- .368 * app_auc_cc + .632 * mean(test_auc_cc, na.rm = T)
  boot632_brier_cc <- .368 * app_brier_cc + .632 * mean(test_brier_cc, na.rm = T)
  
  # .632+ performance

  # Define relative overfitting rate
  r_auc_cc <- (mean(test_auc_cc, na.rm = T) - app_auc_cc) / 
    (gamma_auc - app_auc_cc)
  r_brier_cc <- (mean(test_brier_cc, na.rm = T) - app_brier_cc) / 
    (gamma_brier - app_brier_cc)
  
  # Define weights
  w_auc_cc <- .632 / (1 - .368 * r_auc_cc)
  w_brier_cc <- .632 / (1 - .368 * r_brier_cc)
  
  # Calculate the .632+ AUC and Brier
  boot632plus_auc_cc <- (1 - w_auc_cc) * app_auc_cc + 
    w_auc_cc * mean(test_auc_cc, na.rm = T)
  boot632plus_brier_cc <- (1 - w_brier_cc) * app_brier_cc + 
    w_brier_cc * mean(test_brier_cc, na.rm = T)
  
  # Get predictions from the complete case model above
  preds_cc <- 
    boot_imp_dat |> 
    map(
      function(x){
        tryCatch({
          predict(mod_cc0, 
                  newdata = x |> mutate(time = predtime), 
                  type = "survival")
        },
        error = function(err){
          NA
        })
      }
    )
  
  preds_avg_cc <- inv_cloglog(rowMeans(cloglog(
    as.data.frame(do.call(cbind, preds_cc))
  ), na.rm = TRUE)) 
  
    
  # Bootstrap-impute data
  
  # For impute methods 2 and 3, need complete case data
  imp_dat_bi <- 
    imp_dat |> 
    filter(if_all(all_of(pred_col), ~ !is.na(.)))
  
  boot_imp_dat_bi <- 
    boot_imp_dat |> 
    map(
      ~ .x |> filter(if_all(all_of(pred_col), ~ !is.na(.)))
    )
  
  
  # mod_bi0 will be used to get a full vector of predictions, with missing values
  # fit this to the dataset with all patients included
  if(nrow(imp_dat_bi) == 0) {
    mod_bi0 <- NA
  } else {
    mod_bi0 <- coxph(as.formula(
      paste0("Surv(os_yrs, os_event) ~ ", 
             paste0(complete_var, collapse = " + "), "+", 
             paste0(pred_col, collapse = " + "))), 
      data = imp_dat, x = TRUE, y = TRUE)
  }
  
  mod_bi0_success <- ifelse(isTRUE(is.na(mod_bi0)[[1]]), 0, 1)
  
  # mod_bi will be used for the Score function, which doesn't play well with NAs
  # fit this to the bi complete case dataset
  if(nrow(imp_dat_bi) == 0){
    mod_bi <- NA
  } else {
    mod_bi <- coxph(as.formula(
      paste0("Surv(os_yrs, os_event) ~ ", 
             paste0(complete_var, collapse = " + "), "+", 
             paste0(pred_col, collapse = " + "))), 
      data = imp_dat_bi, x = TRUE, y = TRUE)
  }
  
  # Get AUC and Brier using the complete case model on the original data
  # return NA if there was imputation model failure for any of the predicted covariates
  if(nrow(imp_dat_bi) == 0) {
    app_auc_brier_bi <- NA
    app_auc_bi <- NA
    app_brier_bi <- NA
  } else{
    app_auc_brier_bi <- 
      Score(list("fit" = mod_bi),
            formula = mod_bi[["formula"]],
            data = imp_dat_bi,
            se.fit = FALSE,
            conf.int = FALSE,
            times = predtime)
    
    app_auc_bi <- app_auc_brier_bi[["AUC"]][["score"]][["AUC"]]
    app_brier_bi <- app_auc_brier_bi[["Brier"]][["score"]][["Brier"]][[2]]
  }
  
  # Fit a bootstrap model to each bootstrap sample
  # This returns NA for any bootstrap samples that had an imputation model failure for any predicted covariate
  boot_mod_bi <- 
    boot_imp_dat_bi |> 
    map(
      function(x){
        tryCatch({
          coxph(as.formula(
            paste0("Surv(os_yrs, os_event) ~ ", 
                   paste0(complete_var, collapse = " + "), "+", 
                   paste0(pred_col, collapse = " + "))), 
            data = x, x = TRUE, y = TRUE)}, 
          warning = function(war){
            NA}, 
          error = function(err){
            NA
          })
      }
    ) 
  
  boot_mod_bi_success = as.numeric(sum(!is.na(boot_mod_bi)))
  
  # Get AUC and Brier using each bootstrap model on the bootstrap sample
  # The error part catches the times when there was a imputation model failure in one of the bootstrap samples, and returns NA
  boot_auc_brier_bi <-
    map2(
      boot_imp_dat_bi,
      boot_mod_bi,
      function(x, y){
        tryCatch({
          Score(list("bfit" = y),
                formula = y[["formula"]],
                data = x,
                se.fit = FALSE,
                conf.int = FALSE,
                times = predtime)}, 
          error = function(err){
            NA
          })
      }
    )
  
  boot_auc_brier_bi_success = as.numeric(sum(!is.na(boot_auc_brier_bi)))
  
  b_auc_bi <-
    map_dbl(
      boot_auc_brier_bi,
      function(x){
        tryCatch({
          x[["AUC"]][["score"]][["AUC"]]},
          error = function(err){
            NA
          })
      }
    )
  
  b_brier_bi <-
    map_dbl(
      boot_auc_brier_bi,
      function(x){
        tryCatch({
          x[["Brier"]][["score"]][["Brier"]][[2]]},
          error = function(err){
            NA
          })
      }
    )
  
  # Get AUC and Brier using each bootstrap model on the original data
  # Need to handle case when the original imputed data had an imputation model failure too, in addition to in the bootstrap samples
  if(nrow(imp_dat_bi) == 0){
    o_auc_brier_bi <- NA
    o_auc_bi <- NA
    o_brier_bi <- NA
  } else {
    o_auc_brier_bi <-
      map(
        boot_mod_bi,
        function(x){
          tryCatch({
            Score(list("ofit" = x),
                  formula = x[["formula"]],
                  data = imp_dat_bi,
                  se.fit = FALSE,
                  conf.int = FALSE,
                  times = predtime)}, 
            error = function(err){
              NA
            })
        }
      )
    
    o_auc_bi <-
      map_dbl(
        o_auc_brier_bi,
        function(x){
          tryCatch({
            x[["AUC"]][["score"]][["AUC"]]},
            error = function(err){
              NA
            })
        }
      )
    
    o_brier_bi <-
      map_dbl(
        o_auc_brier_bi,
        function(x){
          tryCatch({
            x[["Brier"]][["score"]][["Brier"]][[2]]},
            error = function(err){
              NA
            })
        }
      )
  }
    
  # Calculate bootstrap optimism corrected AUC and Brier
  # This is only based on bootstrap samples where there wasn't imputation model failure for any of the predictors
  # If the original imputed dataset had imputation model failure, then this will be NA
  boot_auc_bi <- app_auc_bi - mean(b_auc_bi - o_auc_bi, na.rm = TRUE)
  boot_brier_bi <- app_brier_bi - mean(b_brier_bi - o_brier_bi, na.rm = TRUE)
  
  # .632 performance
  
  # Obtain test datasets for each bootstrap sample
  test_dat_bi <- 
    boot_imp_dat_bi |> 
    map(
      function(x){
        if(nrow(x) == 0) {
          x
        } else{
          imp_dat_bi |> filter(!(id %in% unique(x[["id"]])))
        }
      }
    )
  
  # Calculate performance in each test dataset
  test_auc_brier_bi <-
    map2(
      test_dat_bi,
      boot_mod_bi,
      function(x, y){
        tryCatch({
          Score(list("b_fit" = y),
                formula = mod_bi[["formula"]],
                data = x,
                se.fit = FALSE,
                conf.int = FALSE,
                times = predtime)}, 
          error = function(err){
            NA
          })
      }
    )
  
  test_auc_bi <-
    map_dbl(
      test_auc_brier_bi,
      function(x){
        tryCatch({
          x[["AUC"]][["score"]][["AUC"]]},
          error = function(err){
            NA
          })
      }
    )
  
  test_brier_bi <-
    map_dbl(
      test_auc_brier_bi,
      function(x){
        tryCatch({
          x[["Brier"]][["score"]][["Brier"]][[2]]},
          error = function(err){
            NA
          })
      }
    )
  
  # Calculate .632 AUC and Brier
  boot632_auc_bi <- .368 * app_auc_bi + .632 * mean(test_auc_bi, na.rm = T)
  boot632_brier_bi <- .368 * app_brier_bi + .632 * 
    mean(test_brier_bi, na.rm = T) 
  
  # .632+ performance

  # Define relative overfitting rate
  # only based on bootstrap samples where there wasn't imputation model failure for any of the predictors
  r_auc_bi <- (mean(test_auc_bi, na.rm = T) - app_auc_bi) / 
    (gamma_auc - app_auc_bi)
  r_brier_bi <- (mean(test_brier_bi, na.rm = T) - app_brier_bi) / 
    (gamma_brier - app_brier_bi)
  
  # Define weights
  w_auc_bi <- .632 / (1 - .368 * r_auc_bi)
  w_brier_bi <- .632 / (1 - .368 * r_brier_bi)
  
  # Calculate the .632+ AUC and Brier
  boot632plus_auc_bi <- (1 - w_auc_bi) * app_auc_bi + 
    w_auc_bi * mean(test_auc_bi, na.rm = T)
  boot632plus_brier_bi <- (1 - w_brier_bi) * app_brier_bi + 
    w_brier_bi * mean(test_brier_bi, na.rm = T)
  
  # Get predictions from the complete case model above
  # Gives NA for any bootstrap samples where there was imputation model failure
  preds_bi <- 
    boot_imp_dat |> 
    map(
      function(x){
        tryCatch({
          predict(mod_bi0, 
                  newdata = x |> mutate(time = predtime), 
                  type = "survival")
        },
      error = function(err){
        NA
      })
      }
    )

  preds_avg_bi <- inv_cloglog(rowMeans(cloglog(
    as.data.frame(do.call(cbind, preds_bi))
  ), na.rm = TRUE))  

  # Return results
  return(
    list(
      # Full
      app_auc_full = app_auc_full,
      app_brier_full = app_brier_full,
      boot_auc_full = boot_auc_full,
      boot_brier_full = boot_brier_full,
      boot632_auc_full = boot632_auc_full, 
      boot632_brier_full = boot632_brier_full,
      boot632plus_auc_full = boot632plus_auc_full,
      boot632plus_brier_full = boot632plus_brier_full,
      preds_avg_full = preds_avg_full,
      mod_full_success = mod_full_success,
      boot_mod_full_success = boot_mod_full_success,
      boot_auc_brier_full_success = boot_auc_brier_full_success,
      # cc
      app_auc_cc = app_auc_cc,
      app_brier_cc = app_brier_cc,
      boot_auc_cc = boot_auc_cc,
      boot_brier_cc = boot_brier_cc,
      boot632_auc_cc = boot632_auc_cc, 
      boot632_brier_cc = boot632_brier_cc,
      boot632plus_auc_cc = boot632plus_auc_cc,
      boot632plus_brier_cc = boot632plus_brier_cc,
      preds_avg_cc = preds_avg_cc,
      mod_cc0_success = mod_cc0_success, 
      boot_mod_cc_success = boot_mod_cc_success,
      boot_auc_brier_cc_success = boot_auc_brier_cc_success,
      #bi
      app_auc_bi = app_auc_bi,
      app_brier_bi = app_brier_bi,
      boot_auc_bi = boot_auc_bi,
      boot_brier_bi = boot_brier_bi,
      boot632_auc_bi = boot632_auc_bi, 
      boot632_brier_bi = boot632_brier_bi,
      boot632plus_auc_bi = boot632plus_auc_bi,
      boot632plus_brier_bi = boot632plus_brier_bi,
      preds_avg_bi = preds_avg_bi,
      mod_bi0_success = mod_bi0_success, 
      boot_mod_bi_success = boot_mod_bi_success,
      boot_auc_brier_bi_success = boot_auc_brier_bi_success
      )
    )
  
}  


# 5. Function to summarize bias of performance and prediction results -----------
# Input the list returned from res_summary to summary bias
# Example: bias_summary(res)
bias_summary <- function(res) {
  
  cc_bias_res <- 
    tibble(
      auc_orig = res[["app_auc_full"]] - res[["app_auc_cc"]],
      auc_boot = res[["boot_auc_full"]] - res[["boot_auc_cc"]],
      auc_632 = res[["boot632_auc_full"]] - res[["boot632_auc_cc"]],
      auc_632_plus = res[["boot632plus_auc_full"]] - res[["boot632plus_auc_cc"]], 
      brier_orig = res[["app_brier_full"]] - res[["app_brier_cc"]], 
      brier_boot = res[["boot_brier_full"]] - res[["boot_brier_cc"]], 
      brier_632 = res[["boot632_brier_full"]] - res[["boot632_brier_cc"]], 
      brier_632_plus = res[["boot632plus_brier_full"]] - res[["boot632plus_brier_cc"]]
    )
  
  bi_bias_res <-
    tibble(
      auc_orig = res[["app_auc_full"]] - res[["app_auc_bi"]],
      auc_boot = res[["boot_auc_full"]] - res[["boot_auc_bi"]],
      auc_632 = res[["boot632_auc_full"]] - res[["boot632_auc_bi"]],
      auc_632_plus = res[["boot632plus_auc_full"]] - res[["boot632plus_auc_bi"]], 
      brier_orig = res[["app_brier_full"]] - res[["app_brier_bi"]], 
      brier_boot = res[["boot_brier_full"]] - res[["boot_brier_bi"]], 
      brier_632 = res[["boot632_brier_full"]] - res[["boot632_brier_bi"]], 
      brier_632_plus = res[["boot632plus_brier_full"]] - res[["boot632plus_brier_bi"]]
    )
    
  # Results table for prediction
  preds_bias_res <- 
    cbind(res[["preds_avg_full"]], res[["preds_avg_cc"]], 
          res[["preds_avg_bi"]]) |> 
    as_tibble(.name_repair = ~ c("preds_avg_full", "preds_avg_cc", 
                                 "preds_avg_bi")) |> 
    mutate(
      bias_cc = preds_avg_full - preds_avg_cc,
      bias_bi = preds_avg_full - preds_avg_bi
    ) |> 
    summarize(
      avg_pred_bias_cc = mean(bias_cc, na.rm = TRUE),
      avg_pred_bias_bi = mean(bias_bi, na.rm = TRUE)
    )
  
  return(list(cc_bias_res = cc_bias_res,
              bi_bias_res = bi_bias_res,
              preds_bias_res = preds_bias_res))
  
}


# 6. Function to check imputation accuracy -----------------------------------
#
# Binary x: Calculate accuracy
# Continuous x: Calculate mean and sd for x and x_pred
#
# Example: imp_accuracy(impute_list[[49]])
imp_accuracy <- function(imp,  deci = 3) {
  
  results <- list()
  missing_col <- colnames(imp)[colSums(is.na(imp)) > 0]
  origin_col <- sub("_.*", "", missing_col)
  
  for (x in origin_col) {
    x_missing_col <- paste0(x, "_miss")
    x_pred_col <- paste0(x, "_pred")
    
    # Filter rows where x_miss is NA
    filtered_data <- imp[is.na(imp[[x_missing_col]]), ]
    
    # Compare x and x_pred
    if (length(unique(filtered_data[[x]])) == 2) {
      
      # Binary x: Calculate accuracy
      accuracy <- as.character(round(mean(filtered_data[[x]] == filtered_data[[x_pred_col]]),deci))
      results[[x]] <- list(
        accuracy = accuracy
      )
    } else {
      
      # Continuous x: Calculate mean and sd for x and x_pred
      mean_x <- round(mean(filtered_data[[x]], na.rm = TRUE),deci)
      sd_x <- round(sd(filtered_data[[x]], na.rm = TRUE),deci)
      mean_x_pred <- round(mean(filtered_data[[x_pred_col]], na.rm = TRUE),deci)
      sd_x_pred <- round(sd(filtered_data[[x_pred_col]], na.rm = TRUE),deci)
      mean_sd = paste0("original ", mean_x, " (", sd_x, "), imputed ", mean_x_pred, " (", sd_x_pred, ")")
      
      results[[x]] <- list(
        mean_sd = mean_sd
      )
    }
  }
  
  return(results)
}


# 7. Function to check imputation success---------------------------------------------
#
# For each variable ending with '_pred', check imputation success by whether it is fully NA
# Summarize imputation success for a single dataset or across a list of datasets 
#
# Example: count_all_na_pred_unified(impute_list_nsim)
count_all_na_pred_unified <- function(data_or_list) {
  # Check if input is a list of datasets or a single dataset
  is_list_of_datasets <- is.list(data_or_list) && all(map_lgl(data_or_list, is.data.frame))
  
  # Find '_pred' variables
  pred_vars <- grep("_pred$", colnames(if (is_list_of_datasets) data_or_list[[1]] else data_or_list), value = TRUE)
  
  # Define the function to check if each '_pred' variable is fully NA
  pred_na_counts <- map(pred_vars, function(var) {
    if (is_list_of_datasets) {
      # If we have a list of datasets, count datasets where variable is fully NA
      all_na_counts <- map_lgl(data_or_list, ~ all(is.na(.x[[var]])))
      sum(all_na_counts)  # Sum TRUE values representing datasets with all NA for that variable
    } else {
      # If a single dataset, check if the variable is fully NA and return 1 (all NA) or 0 (not all NA)
      as.integer(all(is.na(data_or_list[[var]])))
    }
  })
  
  # Set variable names as names of result and return
  set_names(pred_na_counts, pred_vars)
}


# 99. Function to get the bootstrap-corrected AUC and Brier score, -------------------
#    the .632 AUC, the .632+ AUC and predictions
#
# fit = full cox model fit on original data
# dataset = name of original dataframe use to generate fit
# boot_data_list = list of bootstrapped datasets
# predtime = timepoint for time-dependent AUC and predictions, 
#            defaults to 5 (years)
# Call to coxph must include options x = TRUE and y = TRUE
#
# Example: do_internal_val(mod_full, dat, boot_dat)
do_internal_val <- function(fit, dataset, boot_data_list, predtime = 5) {
  
  # Calculate performance in original data
  app_auc_brier <- 
    Score(list("fit" = fit),
          formula = fit[["formula"]],
          data = dataset,
          se.fit = FALSE,
          conf.int = FALSE,
          times = predtime)
  
  app_auc <- app_auc_brier[["AUC"]][["score"]][["AUC"]]
  app_brier <- app_auc_brier[["Brier"]][["score"]][["Brier"]][[2]]
  
  # Fit a bootstrap model to each bootstrap sample
  boot_mod <- 
    map(
      boot_data_list,
      ~ coxph(fit[["formula"]], data = .x, x = TRUE, y = TRUE)
    ) 
  
  
  # Bootstrap AUC and Brier
  b_auc_brier <-
    map2(
      boot_data_list,
      boot_mod,
      ~ Score(list("b_fit" = .y),
              formula = fit[["formula"]],
              data = .x,
              se.fit = FALSE,
              conf.int = FALSE,
              times = predtime)
    )
  
  b_auc <-
    map_dbl(
      b_auc_brier,
      ~ .x[["AUC"]][["score"]][["AUC"]]
    )
  
  b_brier <-
    map_dbl(
      b_auc_brier,
      ~ .x[["Brier"]][["score"]][["Brier"]][[2]]
    )
  
  # Calculate performance in original data
  o_auc_brier <-
    map(
      boot_mod,
      ~ Score(list("ofit" = .x),
              formula = fit[["formula"]],
              data = dataset,
              se.fit = FALSE,
              conf.int = FALSE,
              times = predtime)
    )
  
  o_auc <-
    map_dbl(
      o_auc_brier,
      ~ .x[["AUC"]][["score"]][["AUC"]]
    )
  
  o_brier <-
    map_dbl(
      o_auc_brier,
      ~ .x[["Brier"]][["score"]][["Brier"]][[2]]
    )
  
  
  # Calculate bootstrap optimism corrected AUC and Brier
  boot_auc <- app_auc - mean(b_auc - o_auc)
  boot_brier <- app_brier - mean(b_brier - o_brier)
  
  # .632 performance
  
  # Obtain test datasets for each bootstrap sample
  test_dat <- 
    boot_data_list |> 
    map(
      ~ dataset |> filter(!(id %in% unique(.x[["id"]])))
    )
  
  # Calculate performance in each test dataset
  test_auc_brier <-
    map2(
      test_dat,
      boot_mod,
      ~ Score(list("b_fit" = .y),
              formula = fit[["formula"]],
              data = .x,
              se.fit = FALSE,
              conf.int = FALSE,
              times = predtime)
    )
  
  test_auc <-
    map_dbl(
      test_auc_brier,
      ~ .x[["AUC"]][["score"]][["AUC"]]
    )
  
  test_brier <-
    map_dbl(
      test_auc_brier,
      ~ .x[["Brier"]][["score"]][["Brier"]][[2]]
    )
  
  # Calculate .632 AUC and Brier
  boot632_auc <- .368 * app_auc + .632 * mean(test_auc)
  boot632_brier <- .368 * app_brier + .632 * mean(test_brier) 
  
  
  # .632+ performance
  
  # Fix "no information performance" for each metric
  gamma_auc <- 0.5
  gamma_brier <- 0.25
  
  # Define relative overfitting rate
  r_auc <- (mean(test_auc) - app_auc) / (gamma_auc - app_auc)
  r_brier <- (mean(test_brier) - app_brier) / (gamma_brier - app_brier)
  
  # Define weights
  w_auc <- .632 / (1 - .368 * r_auc)
  w_brier <- .632 / (1 - .368 * r_brier)
  
  # Calculate the .632+ AUC and Brier
  boot632plus_auc <- (1 - w_auc) * app_auc + w_auc * mean(test_auc)
  boot632plus_brier <- (1 - w_brier) * app_brier + w_brier * mean(test_brier)
  
  
  # Get predictions from the full model across the bootstrap samples
  preds <- 
    boot_data_list |> 
    map(
      ~predict(fit, 
               newdata = .x |> mutate(time = predtime), 
               type = "survival")
    )
  
  preds_avg <- inv_cloglog(rowMeans(cloglog(
    as.data.frame(do.call(cbind, preds))
  )))   
  
  return(
    list(
      app_auc = app_auc,
      app_brier = app_brier,
      boot_auc = boot_auc,
      boot_brier = boot_brier,
      boot632_auc = boot632_auc, 
      boot632_brier = boot632_brier,
      boot632plus_auc = boot632plus_auc,
      boot632plus_brier = boot632plus_brier,
      preds = preds_avg
    )
  )
  
}
