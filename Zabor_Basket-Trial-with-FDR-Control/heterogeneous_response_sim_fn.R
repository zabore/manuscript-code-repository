exact_3basket <- function(pi1, pi2, pi3, sampsize, prior_pep) {
  
  # Create the prior PEP based on the input 
  # Note that this is currently hardcoded to only allow equal entries for all pairwise prior exchangeability probabilities
  prior_pep_mat <-
    matrix(
      c(1, prior_pep, prior_pep, 
        prior_pep, 1, prior_pep, 
        prior_pep, prior_pep, 1), 
    3, 
    3
  )
  
  # Randomly generate subtyple labels from unif
  sub_name <-rep(c("A", "B", "C"), each = sampsize)
  
  # Generate a vector of responses based on obs_rsp_rate and n_sub
  obs_rsp <- map(
    c(pi1, pi2, pi3), 
    ~ rbinom(n = 1, size = sampsize, prob = .x)) %>% 
    flatten_dbl
  
  # Now create a patient-level dataframe
  dat <- 
    tibble::tibble(
      sub_name = c("A", "B", "C"), 
      enrolled = sampsize,
      responders = obs_rsp
    )
  
  # Run the exact basket design
  res <- mem_exact(
    responses = dat[["responders"]],
    size = dat[["enrolled"]],
    name = dat[["sub_name"]], 
    prior = prior_pep_mat)
  
  # Function returns
  return(
    list(
      N = sampsize,
      pi = c(pi1, pi2, pi3),
      prior = prior_pep,
      map = basket_map(res),
      pep = basket_pep(res),
      ess = res$basket$ESS,
      clusters = cluster_baskets(res),
      mean_rr = res$basket$mean_est,
      median_rr = res$basket$median_est
      )
    )
}