# Function to tidy the results of a model fit with icenReg::ic_sp
tidy_ic_sp <- function(x, exponentiate =  FALSE, conf.level = 0.95, ...) {
  tidy <-
    tibble(
      term = names(x[["coefficients"]]),
      estimate = x[["coefficients"]],
      std.error = sqrt(diag(x[["var"]])),
      statistic = summary(x)$summaryParameters[, "z-value"],
      p.value = summary(x)$summaryParameters[, "p"],
      conf.low = confint(x, level = conf.level)[, 1],
      conf.high = confint(x, level = conf.level)[, 2]
    )
  
  if (exponentiate == TRUE)
    tidy <- dplyr::mutate_at(tidy, vars(estimate, conf.low, conf.high), exp)
  
  tidy
}