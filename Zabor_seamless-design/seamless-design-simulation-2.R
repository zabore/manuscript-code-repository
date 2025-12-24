# load packages
library(BOIN)
library(dplyr)
library(janitor)
library(purrr)
library(tidyr)
library(ppseq)

# In this second setting, we have only 3 dose levels and the highest dose is the MTD

# source functions
source(here::here("code", "seamless-design-functions.R"))

# Make a vectorized version of the calc_posterior function to use inside mutate
calc_posterior_v <- Vectorize(calc_posterior)

# fixed parameters
j <- 2 # number of biomarker subgroups
k <- 3 # number of doses
p01 <- 0.1 # null response rate subgroup 1
p11 <- 0.1 # alt response rate subgroup 1
p02 <- 0.1 # null response rate subgroup 2
p12 <- 0.3 # alt response rate subgroup 2
d <- c(0.10, 0.15, 0.25) # true toxicity per dose
d0 <- 0.3 # target toxicity
cohort_size <- 3 # cohort size for BOIN
max_cohorts <- 10 # max cohorts for BOIN
start_dose <- 1 # starting dose level for BOIN
max_cohort_size <- 9 # max cohort size for BOIN

# load dose escalation design
load(here::here("data", "cal-thresh-1.rda"))
load(here::here("data", "opt-eff-design-1.rda"))
load(here::here("data", "eff-rules-1.rda"))
load(here::here("data", "tox-rules-1.rda"))

# This isn't known until the above thresholds are determined
dose_expansion_n <- max(eff_rules_1$n) # sample size per dosing arm in dose expansion phase
n_interim <- seq(5, dose_expansion_n, 5) # Interim look after every 5 subjects

# load efficacy design
load(here::here("data", "cal-thresh-2.rda"))
load(here::here("data", "opt-eff-design-2.rda"))
load(here::here("data", "eff-rules-2.rda"))
load(here::here("data", "tox-rules-2.rda"))

# This isn't known until the above thresholds are determined
eff_expansion_n <- max(eff_rules_2$n0) # sample size per dosing arm in dose expansion phase
n_interim_eff <- seq(25, eff_expansion_n, 5) # Interim look after every 5 subjects
n_for_control <- eff_expansion_n # sample size for control groups is whatever the efficacy_n is
n_for_trt <- eff_expansion_n - dose_expansion_n # sample size for treatment groups is efficacy_n - dose_expansion_n

# Get the BOIN escalation and de-escalation boundaries
# Only need to do this once per d0, cohort_size, max_cohorts, d combination
boin_bounds <-
  get.boundary(
    target = d0,
    ncohort = max_cohorts,
    cohortsize = cohort_size
  )

# Format the boundaries as a tibble
boin_bounds_tab <- 
  as_tibble(
    t(boin_bounds[["boundary_tab"]])
  ) |> 
  clean_names()

# Fix the number of simulations
S <- 5000

# List of response probability matrices
# This is the only thing we need to map over
p <- 
  list(
    matrix(c(0.1, 0.1, 0.1,
             0.1, 0.1, 0.1),
           nrow = 2, byrow = T),
    matrix(c(0.1, 0.1, 0.1,
             0.3, 0.3, 0.3),
           nrow = 2, byrow = T),
    matrix(c(0.1, 0.1, 0.1,
             0.1, 0.3, 0.3),
           nrow = 2, byrow = T),
    matrix(c(0.1, 0.1, 0.1,
             0.1, 0.2, 0.3),
           nrow = 2, byrow = T),
    matrix(c(0.3, 0.3, 0.3,
             0.3, 0.3, 0.3),
           nrow = 2, byrow = T),
    matrix(c(0.1, 0.2, 0.3,
             0.1, 0.2, 0.3),
           nrow = 2, byrow = T),
    matrix(c(0.1, 0.3, 0.3,
             0.1, 0.3, 0.3),
           nrow = 2, byrow = T),
    matrix(c(0.1, 0.2, 0.3,
             0.3, 0.3, 0.3),
           nrow = 2, byrow = T)
  )

# set the seed
set.seed(20251107)


# Dose escalation stage -------------------------------------------------------

# Run the BOIN trials
boins_2 <-
  map(
    1:S,
    ~run_boin_trial(
      start_dose = 1, 
      max_cohorts = 10, 
      cohort_size = 3, 
      d0 = 0.3, 
      n_doses = k)
  )

# Identify the sufficiently safe dose
suff_safe_2 <- 
  map(
    boins_2,
    ~list_safe_doses(.x)
  )

# Save
save(
  boins_2,
  suff_safe_2,
  file = here::here("data", "escalate-sim-res-2.rda")
)


# Dose expansion stage --------------------------------------------------------

# Get the number of sufficiently safe doses
# Don't use this below - need to do everything for 3 doses, with NAs, and then figure out how to summarize later
n_safe_doses_2 <- 
  map(
    suff_safe_2,
    ~length(.x[!is.na(.x)]) 
  )

# generate dataframe with efficacy and toxicity results for each arm
# this is where p comes in 
# need to figure out how to integrate p here: need 1000 simulated datasets per p setting (total of 8 settings)
# Do we use the same S dose escalation results for each value of p in the next stages? That might be a more rigorous comparison? Also easier to program I think...
# This results in a list of length 8, for each value of p
# Then nested in each of those is a list of length S, containing a dose expansion dataframe based on the results of each of the S dose escalation trials
dose_expand_tab_2 <-
  map(
    p, 
    ~ map(
      suff_safe_2,
      function(a) tibble( 
        biomark_grp = rep(1:j, each = 3 * dose_expansion_n),
        dose = rep(rep(a, each = dose_expansion_n), 2) 
      ) |> 
        rowwise() |> 
        # Need to account for possible NAs in suff_safe (i.e. is dose 1 or 2 so there aren't 3 sufficiently safe doses)
        mutate(
          response = ifelse(is.na(dose), NA, rbinom(1, 1, .x[biomark_grp, dose])),
          toxicity = ifelse(is.na(dose), NA, rbinom(1, 1, d[dose]))
        )
    )
  )

# Now summarize the efficacy and toxicity results after every 5 patients
# This results in a list of length 8, for each value of p
# Then nested in each of those is a list of length S, containing a dose expansion dataframe based on the results of each of the S dose escalation trials
dose_expand_summ_2  <- 
  map(
    dose_expand_tab_2,
    ~ .x |> 
      map(function(df) df |> 
            # Should I filter the rows with NA for dose first?
            filter(!is.na(dose)) |> 
      group_by(biomark_grp, dose) |> 
      # Sum up the responses and toxicities within biomarker subgroup and dose
      mutate(
        response_cumsum = cumsum(response),
        toxicity_cumsum = cumsum(toxicity),
        n = row_number()
      ) |> 
      ungroup() |> 
      # Only keep rows that are at the interim analyses
      filter(n %in% n_interim) |> 
      # Add in the decision rules
      full_join(eff_rules_1 |> select(n, r) |> rename(r_eff = r), by = "n") |> 
      full_join(tox_rules_1 |> select(n, r) |> rename(r_tox = r), by = "n") |> 
      # Evaluate decisions to stop for efficacy or toxicity
      mutate(
        stop_eff = ifelse(response_cumsum <= r_eff, TRUE, FALSE),
        stop_tox = ifelse(toxicity_cumsum >= r_tox, TRUE, FALSE)
      ) |> 
      group_by(biomark_grp, dose) |> 
      # Identify first instance of stopping within biomarker/dose group
      mutate(
        first_eff_stop = case_when(
          stop_eff == TRUE & n == min(n) ~ TRUE, # if we stop at the first interim
          stop_eff == TRUE & lag(stop_eff) != TRUE ~ TRUE, 
          .default = FALSE
        ),
        first_tox_stop = case_when(
          stop_tox == TRUE & n == min(n) ~ TRUE, # if we stop at the first interim
          stop_tox == TRUE & lag(stop_tox) != TRUE ~ TRUE,
          .default = FALSE
        ),
        stop = case_when(
          first_eff_stop == TRUE & first_tox_stop == FALSE ~ "Futility",
          first_eff_stop == FALSE & first_tox_stop == TRUE ~ "Toxicity",
          first_eff_stop == TRUE & first_tox_stop == TRUE ~ "Futility & Toxicity"
        )
      ) |> 
      fill(stop, .direction = "down") |> 
      # Then limit columns to ones that occurred at or before stopping
      filter(
        # This keeps the last row for groups that didn't stop early
        (is.na(stop) & n == max(n)) | 
          # This keeps the row with early stopping decision for rows that did stop early
          (!is.na(stop) & is.na(lag(stop)))
      ) |> 
      ungroup() |> 
      select(
        -response, 
        -toxicity, 
        -stop_eff, 
        -stop_tox, 
        -first_eff_stop, 
        -first_tox_stop)
    )
  )

# Rank the doses according to safety and toxicity
# Add the posterior probabilities to the dataframe and calculate the euclidean distance
dose_rank_2 <- 
  map(
    dose_expand_summ_2,
    ~ .x |> 
      map(
        function(df) 
          if(any(is.na(df$stop[!is.na(df$dose)]))){
            df |> 
              # Need to filter out NAs of dose for situations where the MTD was dose 1 or 2 so there are some unneeded rows here with no info
              filter(!is.na(dose) & n == dose_expansion_n & is.na(stop)) |> 
              mutate(
                # Need to add in the null rates according to biomarker group
                # This is pretty manual and specific to this situation right now - would need to automate better for sims
                p0 = case_when(
                  biomark_grp == 1 ~ p01, 
                  biomark_grp == 2 ~ p02
                ),
                d0 = d0,
                post_eff = calc_posterior_v(y = response_cumsum, n = n, p0 = p0),
                post_tox = calc_posterior_v(y = toxicity_cumsum, n = n, p0 = d0, 
                                            direction = "less"),
                dist = sqrt((1 - post_eff)^2 + (1 - post_tox)^2)
              ) |> 
              # And only keep the dose in each biomarker group with the minimal distance
              group_by(
                biomark_grp
              ) |> 
              filter(
                dist == min(dist)
              ) |> 
              ungroup()
          } else {
            tibble(
              biomark_grp = NA,
              dose = NA,
              response_cumsum = NA,
              toxicity_cumsum = NA,
              n = NA, 
              r_eff = NA,
              r_tox = NA,
              stop = NA,
              p0 = NA,
              d0 = NA, 
              post_eff = NA, 
              post_tox = NA, 
              dist = NA
            )
          }
      )
  )

# Save
save(
  n_safe_doses_2,
  dose_expand_tab_2,
  dose_expand_summ_2,
  dose_rank_2,
  file = here::here("data", "expand-sim-res-2.rda")
)


# Efficacy stage --------------------------------------------------------------
# dataframe with efficacy results for control and treatment arm within biomarker subgroup for the stratified control arm design
# This results in a list of length 8, for each value of p
# Then nested in each of those is a list of length S, containing an efficacy dataframe based on the results of the dose expansion and ranking for each of the S simulated datasets 
efficacy_tab_2 <-
  map2(
    dose_rank_2,
    p,
    ~ .x |> 
      map(
        function(df) 
          if(any(!is.na(df$biomark_grp))){ # This is a check for whether no dose levels reached full enrollment in any arm
            tibble(
              biomark_grp = rep(df$biomark_grp, each = n_for_control + n_for_trt),
              dose = rep(df$dose, each = n_for_control + n_for_trt),
              rand_grp = rep(c(rep(0, n_for_control), rep(1, n_for_trt)), 
                             length(df$biomark_grp))
            ) |> 
              rowwise() |> 
              mutate(
                # The response rate is the null response rate for each biomarker-subgroup-specific control group, and the true efficacy rate for the optimal dose in each biomarker subgroup
                response = case_when(
                  rand_grp == 0 & biomark_grp == 1 ~ rbinom(1, 1, p01),
                  rand_grp == 0 & biomark_grp == 2 ~ rbinom(1, 1, p02), 
                  rand_grp == 1 ~ rbinom(1, 1, .y[biomark_grp, dose])
                ),
                toxicity = rbinom(1, 1, d[dose])
              ) |> 
              ungroup()
          } else {
            tibble(
              biomark_grp = NA,
              dose = NA, 
              rand_grp = NA,
              response = NA,
              toxicity = NA
            )
          }
      )
  )

# Now summarize the efficacy and toxicity results
efficacy_summ_2 <- 
  map(
    1:length(p),
    function(p_index){
      map(
        1:S,
        function(S_index) {
          if(any(!is.na(efficacy_tab_2[[p_index]][[S_index]]$biomark_grp))) { # Check for NA
            efficacy_tab_2[[p_index]][[S_index]] |> 
              # Need to add in the information for treated patients from the dose expansion stage
              bind_rows(
                dose_expand_tab_2[[p_index]][[S_index]] |> 
                  right_join(dose_rank_2[[p_index]][[S_index]] |> 
                               select(biomark_grp, dose),
                             by = join_by(biomark_grp, dose)) |> 
                  mutate(rand_grp = 1)) |> 
              group_by(biomark_grp, rand_grp) |> 
              # Sum up the responses within biomarker subgroup and randomization group
              mutate(
                response_cumsum = cumsum(response),
                toxicity_cumsum = cumsum(toxicity),
                n = row_number()
              ) |> 
              ungroup() |> 
              # Only keep rows that are at the interim analyses
              filter(n %in% n_interim_eff) |> 
              select(-response, -toxicity) |> 
              # Now get this into wide format with variables r0_eff, r1_eff, r0_tox, r1_tox
              pivot_wider(
                id_cols = c(biomark_grp, dose, n),
                names_from = rand_grp,
                names_glue = "{.value}_{rand_grp}",
                values_from = c(response_cumsum, toxicity_cumsum)
              ) |> 
              # Add in the decision rules 
              # Only keep response decision table rows that match to the observed number of responses in the control group
              left_join(
                eff_rules_2 |> 
                  select(n0, r0, r1) |> 
                  rename(n = n0, r0_eff = r0, r1_eff = r1),
                by = c("n" = "n", "response_cumsum_0" = "r0_eff")
              ) |>
              # And only keep toxicity decision table rows that match to the observed number of responses in the control group
              left_join(
                tox_rules_2 |>
                  select(n, r) |>
                  rename(r_tox = r),
                by = "n"
              ) |> 
              # Evaluate decisions to stop for efficacy or toxicity
              mutate(
                stop_eff = ifelse(response_cumsum_1 <= r1_eff, TRUE, FALSE),
                stop_tox = ifelse(toxicity_cumsum_1 >= r_tox, TRUE, FALSE)
              ) |> 
              group_by(biomark_grp) |> 
              # Identify first instance of stopping within biomarker/dose group
              mutate(
                first_eff_stop = case_when(
                  stop_eff == TRUE & n == min(n) ~ TRUE, # if we stop at the first interim
                  stop_eff == TRUE & lag(stop_eff) != TRUE ~ TRUE, 
                  .default = FALSE
                ),
                first_tox_stop = case_when(
                  stop_tox == TRUE & n == min(n) ~ TRUE, # if we stop at the first interim
                  stop_tox == TRUE & lag(stop_tox) != TRUE ~ TRUE,
                  .default = FALSE
                ),
                stop = case_when(
                  first_eff_stop == TRUE & first_tox_stop == FALSE ~ "Futility",
                  first_eff_stop == FALSE & first_tox_stop == TRUE ~ "Toxicity",
                  first_eff_stop == TRUE & first_tox_stop == TRUE ~ "Futility & Toxicity"
                )
              ) |> 
              fill(stop, .direction = "down") |> 
              # Then limit columns to ones that occurred at or before stopping
              filter(
                # This keeps the last row for groups that didn't stop early
                (is.na(stop) & n == max(n)) | 
                  # This keeps the row with early stopping decision for rows that did stop early
                  (!is.na(stop) & is.na(lag(stop)))
              ) |> 
              ungroup() |>  
              select(
                -stop_eff, 
                -stop_tox,
                -first_eff_stop, 
                -first_tox_stop)
          } else{
            tibble(
              biomark_grp = NA,
              dose = NA,
              n = NA, 
              response_cumsum_0 = NA,
              response_cumsum_1 = NA,
              toxicity_cumsum_0 = NA,
              toxicity_cumsum_1 = NA,
              r1_eff = NA,
              r_tox = NA, 
              stop = NA
            )
          }
        }
      )
    }
  )
  

# Save 
save(
  efficacy_tab_2,
  efficacy_summ_2,
  file = here::here("data", "efficacy-sim-res-2.rda")
)
  
 
  

