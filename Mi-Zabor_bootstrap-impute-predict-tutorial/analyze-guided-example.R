###############################################################################
# This code was used in the Guided Example
# It bootstraps, imputes, then predicts for a single simulated dataset
###############################################################################

# load needed libraries --------------------------------------------------------
# first run install.packages(".") to install package "." from CRAN if needed
library(purrr)
library(dplyr)
library(riskRegression)
library(survival)


# load the data ----------------------------------------------------------------
# Either run the code file "sim-guided-example-data.R" or load the file "dat0.rda"
# I am loading it from my local GitHub clone
load("D:/manuscript-code-repository/Mi-Zabor_bootstrap-impute-predict-tutorial/dat0.rda")


# Function to impute the three missing variables -------------------------------
# Allows for easier imputation into each bootstrap sample
# dataset = full original dataset
do_imp <- function(dataset) {
  # Imputation models
  imp_x1 <- glm(
    x1 ~ x2 + x5 + x6 + x7 + x8 + x9 + x10 + x11,
    data = dataset,
    family = "binomial")
  
  imp_x3 <- glm(
    x3 ~ x2 + x5 + x6 + x7 + x8 + x9 + x10 + x11,
    data = dataset,
    family = "gaussian")
  
  imp_x4 <- glm(
    x4 ~ x2 + x5 + x6 + x7 + x8 + x9 + x10 + x11,
    data = dataset,
    family = "binomial")
  
  # Predicted value
  x1_pred <- ifelse(
    predict(imp_x1, newdata = dataset, 
            type = "response") > 0.5, 1, 0)
  
  x3_pred <- predict(imp_x3, newdata = dataset)
  
  x4_pred <- ifelse(
    predict(imp_x4, newdata = dataset, 
            type = "response") > 0.5, 1, 0)
  
  # Combine observed and predicted
  dat <-
    dataset |>
    mutate(
      x1_imp = ifelse(is.na(x1), x1_pred, x1),
      x3_imp = ifelse(is.na(x3), x3_pred, x3),
      x4_imp = ifelse(is.na(x4), x4_pred, x4)
    )
  
  return(dat)
}


# impute into the original data ------------------------------------------------
dat <- do_imp(dat0)


# generate bootstrap data and impute into each bootstrap sample ----------------

# set seed for replication
set.seed(20240819)

# Take 500 bootstrap samples of the data
boot_dat <-
  map(1:500,
      ~ slice_sample(dat0, prop = 1, replace = TRUE))

# Then impute into the bootstrap samples
boot_imp_dat <- 
  boot_dat |> 
  map(~ do_imp(.x))


# Apparent performance ---------------------------------------------------------

# Fit a model to the original imputed data
app_mod <- 
  coxph(Surv(s, delta) ~ 
          x1_imp + x2 + x3_imp + x4_imp + x5 + 
          x6 + x7 + x8 + x9 + x10 + x11, 
        data = dat,  x = TRUE, y = TRUE)

# Calculate performance in original imputed data
app_auc_brier <- 
  Score(list("fit" = app_mod),
        formula = app_mod[["formula"]],
        data = dat,
        se.fit = FALSE,
        conf.int = FALSE,
        times = 5)

app_auc <- app_auc_brier[["AUC"]][["score"]][["AUC"]]
app_brier <- 
  app_auc_brier[["Brier"]][["score"]][["Brier"]][[2]]


# Bootstrap-corrected performance ----------------------------------------------

# Fit a model to each bootstrap sample
b_mod <- 
  map(
    boot_imp_dat,
    ~ coxph(app_mod[["formula"]], data = .x, 
            x = TRUE, y = TRUE)) 

# Calculate performance in each bootstrap sample
b_auc_brier <-
  map2(
    boot_imp_dat,
    b_mod,
    ~ Score(list("b_fit" = .y),
            formula = app_mod[["formula"]],
            data = .x,
            se.fit = FALSE,
            conf.int = FALSE,
            times = 5))

b_auc <-
  map_dbl(
    b_auc_brier,
    ~ .x[["AUC"]][["score"]][["AUC"]])

b_brier <-
  map_dbl(
    b_auc_brier,
    ~ .x[["Brier"]][["score"]][["Brier"]][[2]])

# Calculate performance in original data
o_auc_brier <-
  map(
    b_mod,
    ~ Score(list("ofit" = .x),
            formula = app_mod[["formula"]],
            data = dat,
            se.fit = FALSE,
            conf.int = FALSE,
            times = 5))

o_auc <-
  map_dbl(
    o_auc_brier,
    ~ .x[["AUC"]][["score"]][["AUC"]]
  )

o_brier <-
  map_dbl(
    o_auc_brier,
    ~ .x[["Brier"]][["score"]][["Brier"]][[2]])

# Calculate bootstrap-corrected AUC and Brier
boot_auc <- app_auc - mean(b_auc - o_auc)
boot_brier <- app_brier - mean(b_brier - o_brier)


# .632 performance -------------------------------------------------------------

# Obtain test datasets for each bootstrap sample
test_dat <- 
  boot_imp_dat |> 
  map(
    ~ dat |> filter(!(id %in% unique(.x[["id"]]))))

# Calculate performance in each test dataset
test_auc_brier <-
  map2(
    test_dat,
    b_mod,
    ~ Score(list("b_fit" = .y),
            formula = app_mod[["formula"]],
            data = .x,
            se.fit = FALSE,
            conf.int = FALSE,
            times = 5))

test_auc <-
  map_dbl(
    test_auc_brier,
    ~ .x[["AUC"]][["score"]][["AUC"]])

test_brier <-
  map_dbl(
    test_auc_brier,
    ~ .x[["Brier"]][["score"]][["Brier"]][[2]])

# Calculate .632 AUC and Brier
boot632_auc <- .368 * app_auc + .632 * mean(test_auc)
boot632_brier <- .368 * app_brier + .632 * mean(test_brier) 


# .632+ performance ------------------------------------------------------------

# Fix "no information performance" for each metric
gamma_auc <- 0.5
gamma_brier <- 0.25

# Define relative overfitting rate
r_auc <- (mean(test_auc) - app_auc) / (gamma_auc - app_auc)
r_brier <- (mean(test_brier) - app_brier) / 
  (gamma_brier - app_brier)

# Define weights
w_auc <- .632 / (1 - .368 * r_auc)
w_brier <- .632 / (1 - .368 * r_brier)

# Calculate the .632+ AUC and Brier
boot632plus_auc <- (1 - w_auc) * app_auc + 
  w_auc * mean(test_auc)
boot632plus_brier <- (1 - w_brier) * app_brier + 
  w_brier * mean(test_brier)