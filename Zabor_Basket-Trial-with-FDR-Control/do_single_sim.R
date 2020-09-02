# Function to generate observerd response vector and run the MEM
do_single_sim <- function(ss, trr, subs) {
  # Generate the observed response count
  or <- 
    map2(
      ss,
      trr,
      ~ rbinom(n = 1, size = .x, prob = .y)
    ) %>% 
    flatten_dbl()
  
  # Run the MEM
  mod <- mem_mcmc(
    responses = or,
    size = ss,
    name = subs, 
    p0 = 0.1,
    prior = diag(10)/(1/.75) + matrix(0.25, nrow = 10, ncol = 10)
  )
  
  # organize the results
  basket_res <- tibble(
    basket = subs,
    pp = mod$basket[["post.prob"]],
    hpd_l = mod$basket[["HPD"]][1, ],
    hpd_u = mod$basket[["HPD"]][2, ],
    mean_est = mod$basket[["mean_est"]],
    median_est = mod$basket[["median_est"]],
    ess = mod$basket[["ESS"]]
  )
  
  cluster_res <- tibble(
    cluster = mod$cluster[["name"]],
    pp = mod$cluster[["post.prob"]],
    hpd_l = mod$cluster[["HPD"]][1, ],
    hpd_u = mod$cluster[["HPD"]][2, ],
    mean_est = mod$cluster[["mean_est"]],
    median_est = mod$cluster[["median_est"]]
  )
  
  # Do the frequentist per-basket tests
  p_null <- NULL
  p_alt <- NULL
  exact_mean <- NULL
  exact_lci <- NULL
  exact_uci <- NULL
  for(i in 1:length(or)) {
    p_null[i] <- 
      binom.test(or[i], ss[i], p = 0.1, alternative = "greater")$p.value
    p_alt[i] <- 
      binom.test(or[i], ss[i], p = 0.3, alternative = "greater")$p.value
    ci <- binom.confint(or[i], ss[i], method = "exact")
    exact_mean[i] <- ci$mean
    exact_lci[i] <- ci$lower
    exact_uci[i] <- ci$upper
  }
  
  # organize the results
  frequentist_res <- tibble(
    basket = subs,
    p_null = p_null,
    p_alt = p_alt,
    exact_mean = exact_mean,
    exact_lci = exact_lci,
    exact_uci = exact_uci
  )
  
  # return the resulting dataframes
  return(list(
    basket = basket_res,
    cluster = cluster_res,
    frequentist_res = frequentist_res
  ))
}



# # Test with just 3 baskets for speed
# 
# # Fix the sample size vector
# ss <- c(26, 25, 18)
# 
# # Fix the true response rate vector
# trr <- rep(0.1, 3)
# 
# # Fix the basket names
# subs <- paste0("Basket", seq(1:3))
# 
# set.seed(20190925)
# 
# or <-
#   map2(
#     ss,
#     trr,
#     ~ rbinom(n = 1, size = .x, prob = .y)
#   ) %>%
#   flatten_dbl()
# 
# # Run the MEM
# mod <- mem_mcmc(
#   responses = or,
#   size = ss,
#   name = subs,
#   p0 = 0.1,
#   prior = diag(3)/(1/.75) + matrix(0.25, nrow = 3, ncol = 3)
# )
# 
# # Run the frequentist
# p_null <- NULL
# p_alt <- NULL
# exact_mean <- NULL
# exact_lci <- NULL
# exact_uci <- NULL
# for(i in 1:length(or)) {
#   p_null[i] <-
#     binom.test(or[i], ss[i], p = 0.1, alternative = "greater")$p.value
#   p_alt[i] <-
#     binom.test(or[i], ss[i], p = 0.3, alternative = "greater")$p.value
#   ci <- binom::binom.confint(or[i], ss[i], method = "exact")
#   exact_mean[i] <- ci$mean
#   exact_lci[i] <- ci$lower
#   exact_uci[i] <- ci$upper
# }
# 
# i <- 1
# prop.test(or[i], ss[i], p = 0.1, correct = FALSE)
# 
# 
# do_single_sim(ss, trr, subs)
