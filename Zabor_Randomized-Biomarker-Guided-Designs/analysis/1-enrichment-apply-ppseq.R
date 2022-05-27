library(dplyr)
library(purrr)
library(tidyr)
library(ppseq)
library(tibble)
library(future)
library(furrr)


# Let's consider the results of the pooled data analysis (stage 1) and the stratified datasets as pairs of data
# Paired by row/simulation number basically

# load the pooled simulation results
load(here::here("code", "pooled-control-group", "p-ppseq-res.rda"))


# 1. remove the thresholds that we aren't actually considering from the pooled results
p_res_ppseq_null_sub <-
  p_res_ppseq_null %>% 
  filter(!pp_threshold %in% c(0, 0.999, 0.9999, 0.99999, 1))

p_res_ppseq_alt_sub <-
  p_res_ppseq_alt %>% 
  filter(!pp_threshold %in% c(0, 0.999, 0.9999, 0.99999, 1))


# 2. Now arrange and filter the stage 1 results to only keep positive subgroups
# We don't care about the control group from stage 1 as it won't be carried forward
# Only care about subgroups that reach full enrollment, i.e. were positive in stage 1
# Want to know the number of patients (will always be 50) and the number of responses

# null
e_res_null_long <- 
  full_join(
    p_res_ppseq_null_sub %>% 
      select(sim_num, pp_threshold, ppp_threshold,
             starts_with("n")) %>% 
      pivot_longer(
        cols = starts_with("n1"),
        names_prefix = "n1_",
        names_to = "Subgroup",
        values_to = "n1"
      ),
    p_res_ppseq_null_sub %>% 
      select(sim_num, pp_threshold, ppp_threshold,
             starts_with("y")) %>% 
      pivot_longer(
        cols = starts_with("y1"),
        names_prefix = "y1_",
        names_to = "Subgroup",
        values_to = "y1"
      )
  ) %>% 
  filter(n1 == 50) %>% 
  select(sim_num, pp_threshold, ppp_threshold, Subgroup, n1, y1) %>% 
  rename(n1_stage1 = n1, y1_stage1 = y1)

# alt
e_res_alt_long <- 
  full_join(
    p_res_ppseq_alt_sub %>% 
      select(sim_num, pp_threshold, ppp_threshold,
             starts_with("n")) %>% 
      pivot_longer(
        cols = starts_with("n1"),
        names_prefix = "n1_",
        names_to = "Subgroup",
        values_to = "n1"
      ),
    p_res_ppseq_alt_sub %>% 
      select(sim_num, pp_threshold, ppp_threshold,
             starts_with("y")) %>% 
      pivot_longer(
        cols = starts_with("y1"),
        names_prefix = "y1_",
        names_to = "Subgroup",
        values_to = "y1"
      )
  ) %>% 
  filter(n1 == 50) %>% 
  select(sim_num, pp_threshold, ppp_threshold, Subgroup, n1, y1) %>% 
  rename(n1_stage1 = n1, y1_stage1 = y1)


# 3. Unlist the stage 2 data - stack with a simnum indicator
# And merge with the information from above to sum stage 1 results into stage 2 data
# These resulting datasets will all have different numbers of rows
# It depends on which stage 1 datasets for each combo of pp and ppp were positive

# load the stratified simulated datasets
load(here::here("code", "stratified-control-group", "s-sim-dat.rda"))

prep_stage2_dat <- function(stage2_df, stage1_res) {
  stage2_df %>% 
    bind_rows(
      .id = "sim_num"
    ) %>% 
    inner_join(stage1_res) %>% 
    arrange(sim_num, pp_threshold, ppp_threshold) %>% 
    mutate(n1 = n1 + n1_stage1,
           y1 = y1 + y1_stage1) %>% 
    select(sim_num, pp_threshold, ppp_threshold, n0, n1, y0, y1)
}

s_sim_dat_trt0_null_sub <- 
  prep_stage2_dat(
    s_sim_dat_trt0_null, 
    e_res_null_long %>% filter(Subgroup == 0)
    )

s_sim_dat_trt1_null_sub <- 
  prep_stage2_dat(
    s_sim_dat_trt1_null, 
    e_res_null_long %>% filter(Subgroup == 1)
  )

s_sim_dat_trt2_null_sub <- 
  prep_stage2_dat(
    s_sim_dat_trt2_null, 
    e_res_null_long %>% filter(Subgroup == 2)
  )

s_sim_dat_trt0_alt_sub <- 
  prep_stage2_dat(
    s_sim_dat_trt0_alt, 
    e_res_alt_long %>% filter(Subgroup == 0)
  )

s_sim_dat_trt1_alt_sub <- 
  prep_stage2_dat(
    s_sim_dat_trt1_alt, 
    e_res_alt_long %>% filter(Subgroup == 1)
  )

s_sim_dat_trt2_alt_sub <- 
  prep_stage2_dat(
    s_sim_dat_trt2_alt, 
    e_res_alt_long %>% filter(Subgroup == 2)
  )


# 4. Now evaluate the thresholds
eval_thresh <- function(data, 
                        pp_threshold, 
                        ppp_threshold, 
                        N = c(50, 100)) {
  decision <- NULL
  ppp <- NULL
  for (i in 1:nrow(data)) {
    ppp[i] <- calc_predictive(
      y = c(data$y0[i], data$y1[i]),
      n = c(data$n0[i], data$n1[i]),
      p0 = NULL, 
      N = N, 
      direction = "greater", 
      delta = 0,
      prior = c(0.5, 0.5), 
      S = 5000,
      theta = pp_threshold
    )
    decision[i] <- ppp[i] < ppp_threshold
    if (decision[i] == TRUE) break
  }
  
  res <- add_column(
    data[ifelse(any(decision == TRUE),
                which(decision == TRUE),
                length(decision)
    ), ],
    ppp = ppp[ifelse(any(decision == TRUE),
                     which(decision == TRUE),
                     length(decision)
    )]) %>% 
    mutate(
      positive = case_when(
        sum(n0, n1) == sum(N) & ppp > pp_threshold ~ TRUE,
        sum(n0, n1) == sum(N) & ppp <= pp_threshold ~ FALSE,
        sum(n0, n1) != sum(N) ~ FALSE
    )
  )
  
  return(res)
}

plan(multicore(workers = 100))
set.seed(20220517)

e_res_ppseq_trt0_null <- 
  s_sim_dat_trt0_null_sub %>% 
  group_by(sim_num, pp_threshold, ppp_threshold) %>% 
  group_split() %>% 
  future_map_dfr(
    ~ eval_thresh(.x, unique(.x$pp_threshold), unique(.x$ppp_threshold), c(50, 100)), 
    .options = furrr_options(seed = TRUE)
    )

e_res_ppseq_trt1_null <- 
  s_sim_dat_trt1_null_sub %>% 
  group_by(sim_num, pp_threshold, ppp_threshold) %>% 
  group_split() %>% 
  future_map_dfr(
    ~ eval_thresh(.x, unique(.x$pp_threshold), unique(.x$ppp_threshold), c(50, 100)), 
    .options = furrr_options(seed = TRUE)
  )

e_res_ppseq_trt2_null <- 
  s_sim_dat_trt2_null_sub %>% 
  group_by(sim_num, pp_threshold, ppp_threshold) %>% 
  group_split() %>% 
  future_map_dfr(
    ~ eval_thresh(.x, unique(.x$pp_threshold), unique(.x$ppp_threshold), c(50, 100)), 
    .options = furrr_options(seed = TRUE)
  )

e_res_ppseq_trt0_alt <- 
  s_sim_dat_trt0_alt_sub %>% 
  group_by(sim_num, pp_threshold, ppp_threshold) %>% 
  group_split() %>% 
  future_map_dfr(
    ~ eval_thresh(.x, unique(.x$pp_threshold), unique(.x$ppp_threshold), c(50, 100)), 
    .options = furrr_options(seed = TRUE)
  )

e_res_ppseq_trt1_alt <- 
  s_sim_dat_trt1_alt_sub %>% 
  group_by(sim_num, pp_threshold, ppp_threshold) %>% 
  group_split() %>% 
  future_map_dfr(
    ~ eval_thresh(.x, unique(.x$pp_threshold), unique(.x$ppp_threshold), c(50, 100)), 
    .options = furrr_options(seed = TRUE)
  )

e_res_ppseq_trt2_alt <- 
  s_sim_dat_trt2_alt_sub %>% 
  group_by(sim_num, pp_threshold, ppp_threshold) %>% 
  group_split() %>% 
  future_map_dfr(
    ~ eval_thresh(.x, unique(.x$pp_threshold), unique(.x$ppp_threshold), c(50, 100)), 
    .options = furrr_options(seed = TRUE)
  )
  

# 5. save the results
save(
  e_res_ppseq_trt0_null, 
  e_res_ppseq_trt1_null,
  e_res_ppseq_trt2_null,
  e_res_ppseq_trt0_alt, 
  e_res_ppseq_trt1_alt,
  e_res_ppseq_trt2_alt,
  file = here::here("code", "enrichment-design", "e-ppseq-stage2-comb-res.rda")
  )
