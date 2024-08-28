# Function to impute the three missing variables -------------------------------
# Allows for easier imputation into each bootstrap sample
# dataset = full original dataset
do_imp <- function(dataset) {
  # Imputation models
  imp_x1 <- glm(x1_miss ~ x2 + x5 + x6 + x7 + x8 + 
                  x9 + x10 + x11,
                data = dataset[!is.na(dataset$x1_miss), ],
                family = "binomial")
  
  imp_x3 <- glm(x3_miss ~ x2 + x5 + x6 + x7 + x8 +
                  x9 + x10 + x11,
                data = dataset[!is.na(dataset$x3_miss), ],
                family = "gaussian")
  
  imp_x4 <- glm(x4_miss ~ x2 + x5 + x6 + x7 + x8 + 
                  x9 + x10 + x11,
                data = dataset[!is.na(dataset$x4_miss), ],
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
      x1_imp = ifelse(is.na(x1_miss), x1_pred, x1),
      x3_imp = ifelse(is.na(x3_miss), x3_pred, x3),
      x4_imp = ifelse(is.na(x4_miss), x4_pred, x4)
    )
  
  return(dat)
}