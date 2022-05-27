library(dplyr)
library(purrr)
library(tibble)
library(future)
library(furrr)
library(ppseq)


# load the simulated datasets
load(here::here("code", "pooled-control-group", "p-sim-dat.rda"))


# function to evaluate a combination of thresholds 
# from ppseq::calibrate_thresholds(), but customized for this setting of a pooled control group
# returns a dataframe with each interim look of the trial, the thresholds, the ppp, and whether or not the trial was positive
eval_thresh_pooled <- function(df_ctl, 
                               df_trt0, 
                               df_trt1,
                               df_trt2,
                               pp_threshold, 
                               ppp_threshold,
                               p0,
                               N,
                               direction,
                               delta,
                               prior,
                               S) {
  
  decision0 <- decision1 <- decision2 <- FALSE
  ppp0 <- ppp1 <- ppp2 <- NULL
  
  for (i in 1:nrow(df_ctl)) {
    # print(i)
    
    # calculate the predictive probability for a given posterior threshold
    # do this separately for each treatment group
    if(tail(decision0, n = 1) == FALSE){
      ppp0[i] <- calc_predictive(
        y = c(df_ctl$y1[i], df_trt0$y1[i]),
        n = c(df_ctl$n1[i], df_trt0$n1[i]),
        p0 = p0,
        N = N,
        direction = direction, 
        delta = delta,
        prior = prior, 
        S = S, 
        theta = pp_threshold
      )
      decision0[i] <- ppp0[i] < ppp_threshold
      if(i == nrow(df_ctl)) decision0[i] <- TRUE
    } 
    
    if(tail(decision1, n = 1) == FALSE){
      ppp1[i] <- calc_predictive(
        y = c(df_ctl$y1[i], df_trt1$y1[i]),
        n = c(df_ctl$n1[i], df_trt1$n1[i]),
        p0 = p0,
        N = N,
        direction = direction, 
        delta = delta,
        prior = prior, 
        S = S, 
        theta = pp_threshold
      )
      decision1[i] <- ppp1[i] < ppp_threshold
      if(i == nrow(df_ctl)) decision1[i] <- TRUE
    } 
    
    if(tail(decision2, n = 1) == FALSE){
      ppp2[i] <- calc_predictive(
        y = c(df_ctl$y1[i], df_trt2$y1[i]),
        n = c(df_ctl$n1[i], df_trt2$n1[i]),
        p0 = p0,
        N = N,
        direction = direction, 
        delta = delta,
        prior = prior, 
        S = S, 
        theta = pp_threshold
      )
      decision2[i] <- ppp2[i] < ppp_threshold
      if(i == nrow(df_ctl)) decision2[i] <- TRUE
    } 
    
    # stop the whole thing if all arms have stopped
    if(tail(decision0, n = 1) == TRUE & 
       tail(decision1, n = 1) == TRUE & 
       tail(decision2, n = 1) == TRUE) break
  }
  
  # only keep the rows of the data before the trial stopped
  # add in the thresholds and the ppp for each comparison
  res0 <- cbind(
    df_ctl[max(length(decision0), length(decision1), length(decision2)), ] %>% 
      rename(n0 = n1, y0 = y1),
    df_trt0[length(decision0), ] %>% rename(n1_0 = n1, y1_0 = y1),
    df_trt1[length(decision1), ] %>% rename(n1_1 = n1, y1_1 = y1),
    df_trt2[length(decision2), ] %>% rename(n1_2 = n1, y1_2 = y1)) %>% 
    add_column(
      pp_threshold = pp_threshold,
      ppp_threshold = ppp_threshold,
      ppp0 = ppp0[length(decision0)],
      ppp1 = ppp1[length(decision1)],
      ppp2 = ppp2[length(decision2)]
    )
  
  # add a column recording whether the trial was positive or not
  # if we went to the full sample size and the ppp is > given posterior threshold
  res <- mutate(
    res0,
    positive0 = case_when(
      sum(n0, n1_0) == sum(N) & ppp0 > pp_threshold ~ TRUE,
      sum(n0, n1_0) == sum(N) & ppp0 <= pp_threshold ~ FALSE,
      sum(n0, n1_0) != sum(N) ~ FALSE
    ),
    positive1 = case_when(
      sum(n0, n1_1) == sum(N) & ppp1 > pp_threshold ~ TRUE,
      sum(n0, n1_1) == sum(N) & ppp1 <= pp_threshold ~ FALSE,
      sum(n0, n1_1) != sum(N) ~ FALSE
    ),
    positive2 = case_when(
      sum(n0, n1_2) == sum(N) & ppp2 > pp_threshold ~ TRUE,
      sum(n0, n1_2) == sum(N) & ppp2 <= pp_threshold ~ FALSE,
      sum(n0, n1_2) != sum(N) ~ FALSE
    )
  )
  
  return(res)
}


# function to generate the result dataframes
# this returns a dataframe of results for all the different threshold combinations for each simulated dataset
gen_ppseq_res <- function(df_ctl, df_trt0, df_trt1, df_trt2) {
  
  res <-
    future_pmap(
      list(df_ctl, df_trt0, df_trt1, df_trt2),
      function(w, x, y, z) 
        future_pmap_dfr(
          cross_threshold,
          function(pp_threshold, ppp_threshold) 
            eval_thresh_pooled(
              w, x, y, z,
              pp_threshold, 
              ppp_threshold,
              p0 = p0, 
              N = N,
              direction = direction, 
              delta = delta, 
              prior = prior, 
              S = S), 
          .options = furrr_options(seed = TRUE)
        ), 
      .options = furrr_options(seed = TRUE)
    )
  
  dplyr::bind_rows(
    res,
    .id = "sim_num"
  )
}


# simulation setup -------------------------------------------------------------
plan(multicore(workers = 100))
set.seed(20220404)
N <- c(50, 50)
p0 <- NULL 
direction <- "greater"
delta <- 0 
prior <- c(0.5, 0.5) 
S <- 5000


# make a tibble of the cross-tabulated posterior and predictive thresholds -----
pp_threshold <- c(0, 0.7, 0.74, 0.78, 0.82, 0.86, 
                  0.9, 0.92, 0.93, 0.94, 0.95,
                  0.96, 0.97, 0.98, 0.99, 0.999, 
                  0.9999, 0.99999, 1)
ppp_threshold <- seq(0.05, 0.2, 0.05)

cross_threshold <- 
  purrr::cross_df(list(pp_threshold = pp_threshold, 
                       ppp_threshold = ppp_threshold))


# generate null results --------------------------------------------------------
p_res_ppseq_null <- gen_ppseq_res(p_sim_dat_ctl_null, 
                                  p_sim_dat_trt0_null, 
                                  p_sim_dat_trt1_null, 
                                  p_sim_dat_trt2_null)

# generate alternative results -------------------------------------------------
p_res_ppseq_alt <- gen_ppseq_res(p_sim_dat_ctl_alt, 
                                 p_sim_dat_trt0_alt, 
                                 p_sim_dat_trt1_alt, 
                                 p_sim_dat_trt2_alt)


# save the raw results ---------------------------------------------------------
# will still need to summarize and format these later
save(
  p_res_ppseq_null, 
  p_res_ppseq_alt, 
  file = here::here("code", "pooled-control-group", "p-ppseq-res.rda")
  )