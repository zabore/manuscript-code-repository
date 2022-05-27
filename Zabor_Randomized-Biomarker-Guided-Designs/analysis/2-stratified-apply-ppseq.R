library(dplyr)
library(purrr)
library(tibble)
library(future)
library(furrr)
library(ppseq)


# load the simulated datasets
load(here::here("code", "stratified-control-group", "s-sim-dat.rda"))


# function to evaluate a combination of thresholds - from ppseq::calibrate_thresholds()
eval_thresh <- function(data, 
                        pp_threshold, 
                        ppp_threshold, 
                        p0, 
                        N, 
                        direction = "greater", 
                        delta = NULL,
                        prior = c(0.5, 0.5), 
                        S = 5000) {
  decision <- NULL
  ppp <- NULL
  for (i in 1:nrow(data)) {
    if (ncol(data) == 4) {
      ppp[i] <- calc_predictive(
        y = c(data$y0[i], data$y1[i]),
        n = c(data$n0[i], data$n1[i]),
        direction = direction, 
        p0 = p0, 
        delta = delta,
        prior = prior, 
        S = S, 
        N = N,
        theta = pp_threshold
      )
    } else if (ncol(data) == 2) {
      ppp[i] <- calc_predictive(
        y = data$y1[i], 
        n = data$n1[i],
        direction = direction, 
        p0 = p0, 
        delta = delta,
        prior = prior, 
        S = S, 
        N = N,
        theta = pp_threshold
      )
    }
    decision[i] <- ppp[i] < ppp_threshold
    if (decision[i] == TRUE) break
  }
  res0 <- add_column(
    data[ifelse(any(decision == TRUE),
                which(decision == TRUE),
                length(decision)
    ), ],
    pp_threshold = pp_threshold,
    ppp_threshold = ppp_threshold,
    ppp = ppp[ifelse(any(decision == TRUE),
                     which(decision == TRUE),
                     length(decision)
    )]
  )
  
  if (ncol(data) == 4) {
    res <- mutate(
      res0,
      positive = case_when(
        sum(n0, n1) == sum(N) & ppp > pp_threshold ~ TRUE,
        sum(n0, n1) == sum(N) & ppp <= pp_threshold ~ FALSE,
        sum(n0, n1) != sum(N) ~ FALSE
      )
    )
  } else if (ncol(data) == 2) {
    res <- mutate(
      res0,
      positive = case_when(
        n1 == N & ppp > pp_threshold ~ TRUE,
        n1 == N & ppp <= pp_threshold ~ FALSE,
        n1 != N ~ FALSE
      )
    )
  }
  return(res)
}


# function to generate the result dataframes
gen_ppseq_res <- function(df) {
  
  res <-
    future_map(
      df,
      function(x) 
        future_pmap_dfr(
          cross_threshold,
          function(pp_threshold, ppp_threshold) 
            eval_thresh(
              x, 
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
set.seed(20220405)
p0 <- NULL 
N <- c(50, 50)
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
s_res_ppseq_trt0_null <- gen_ppseq_res(s_sim_dat_trt0_null)
s_res_ppseq_trt1_null <- gen_ppseq_res(s_sim_dat_trt1_null)
s_res_ppseq_trt2_null <- gen_ppseq_res(s_sim_dat_trt2_null)


# generate alternative results -------------------------------------------------
s_res_ppseq_trt0_alt <- gen_ppseq_res(s_sim_dat_trt0_alt)
s_res_ppseq_trt1_alt <- gen_ppseq_res(s_sim_dat_trt1_alt)
s_res_ppseq_trt2_alt <- gen_ppseq_res(s_sim_dat_trt2_alt)


# save the raw results ---------------------------------------------------------
# will still need to summarize and format these later
save(
  s_res_ppseq_trt0_null, s_res_ppseq_trt1_null, s_res_ppseq_trt2_null,
  s_res_ppseq_trt0_alt, s_res_ppseq_trt1_alt, s_res_ppseq_trt2_alt,
  file = here::here("code", "stratified-control-group", "s-ppseq-res.rda")
  )