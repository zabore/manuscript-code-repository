# Functions needed for simulation study of seamless design


# Function to determine next dose level in the BOIN trial
get_next_dose <- function(decision, current_dose, eliminate, n_doses) {
  # If decision is escalate 
  # And not at highest dose
  # And the next highest dose hasn't been eliminated
  # Escalate
  if(decision == "escalate" & 
     current_dose < n_doses & 
     current_dose + 1 != eliminate) {
    current_dose <- current_dose + 1
  } else 
    # If decision is escalate and at highest dose
    # OR if decision is escalate and next dose has been eliminated
    # Stay at current dose
    if(decision == "escalate" & 
       current_dose == n_doses |
       (current_dose < n_doses & (current_dose + 1) == eliminate)) {
      current_dose <- current_dose
    } else 
      # If decision is deescalate 
      # And not at lowest dose
      # deescalate
      if(decision == "deescalate" & 
         current_dose != 1) {
        current_dose <- current_dose - 1
      } else
        # If decision is deescalate 
        # And at lowest dose
        # Stay at current dose
        if(decision == "deescalate" & 
           current_dose == 1) {
          current_dose <- current_dose
        } else 
          # If decision is eliminate 
          # And not at lowest dose
          # Deescalate
          if(decision == "eliminate" & 
             current_dose != 1) {
            eliminate <- current_dose
            current_dose <- current_dose - 1
          } else 
            # If deicision is eliminate 
            # And at lowest dose
            # End trial  
            if(decision == "eliminate" & 
               current_dose == 1) {
              eliminate <- current_dose
              current_dose <- 0
            } else
              # If decision is stay
              # Stay at current dose
              if(decision == "stay") {
                current_dose <- current_dose
              }
}


# Function to simulate a single BOIN trial and select the MTD
run_boin_trial <- function(start_dose, max_cohorts, cohort_size, d0, n_doses) {
  
  current_dose <- start_dose # initialize current dose at the starting dose
  cum_n <- 0 # initialize cumulative sample size 
  cum_n_d <- rep(0, n_doses) # initialize cumulative sample size per dose
  eliminate <- 0 # initialize lowest elimination dose
  cum_tox_d <- rep(0, n_doses) # initialize toxicities per dose level
  
  while(eliminate != 1 & 
        cum_n < max_cohorts * cohort_size &
        cum_n_d[current_dose] < max_cohort_size) {
    
    # increment cumulative sample size
    cum_n <- cum_n + cohort_size 
    
    # get vector of number of patients to add to cumulative count
    add_n_d <- rep(0, n_doses)
    add_n_d[current_dose] <- cohort_size
    
    # new cumulative patients per dose
    cum_n_d <- cum_n_d + add_n_d
    
    # generate number toxicities at current dose
    tox <- rbinom(1, cohort_size, d[current_dose]) 
    
    # get vector of number of toxicities to add to cumulative count
    add_tox_d <- rep(0, n_doses)
    add_tox_d[current_dose] <- tox
    
    # new cumulative toxicities per dose
    cum_tox_d <- cum_tox_d + add_tox_d
    
    # make a decision
    decision <- 
      case_when(
        cum_tox_d[current_dose] <= 
          boin_bounds_tab$escalate_if_number_of_dlt[
            boin_bounds_tab$number_of_patients_treated == 
              cum_n_d[current_dose]] ~ 
          "escalate",
        cum_tox_d[current_dose] >= 
          boin_bounds_tab$deescalate_if_number_of_dlt[
            boin_bounds_tab$number_of_patients_treated == 
              cum_n_d[current_dose]] ~ 
          "deescalate",
        cum_tox_d[current_dose] >= 
          boin_bounds_tab$eliminate_if_number_of_dlt[
            boin_bounds_tab$number_of_patients_treated == 
              cum_n_d[current_dose]] ~ 
          "eliminate",
        .default = "stay"
      ) 
    
    # Set next dose level
    current_dose <- get_next_dose(
      decision = decision, 
      current_dose = current_dose, 
      eliminate = eliminate, 
      n_doses = n_doses)
    
  }
  
  # Select the MTD
  sel_mtd <- 
    select.mtd(
      target = d0,
      npts = cum_n_d,
      ntox = cum_tox_d
    )
  
  return(list(cum_n_d = cum_n_d,
              cum_tox_d = cum_tox_d,
              sel_mtd = sel_mtd)
  )
  
}


# Function to identify the sufficiently safe doses
# This returns a vector of length 3
# A vector of (NA, NA, NA) means the MTD was below the lowest dose
list_safe_doses <- function(boin_sel_mtd_res) {
  suff_safe <- 
    case_when(
      boin_sel_mtd_res[["sel_mtd"]][["MTD"]] >= 3 & 
        boin_sel_mtd_res[["sel_mtd"]][["MTD"]] < 99 ~ 
        c(boin_sel_mtd_res[["sel_mtd"]][["MTD"]] - 2, 
          boin_sel_mtd_res[["sel_mtd"]][["MTD"]] - 1, 
          boin_sel_mtd_res[["sel_mtd"]][["MTD"]]),
      boin_sel_mtd_res[["sel_mtd"]][["MTD"]] == 2 ~ 
        c(NA, 
          boin_sel_mtd_res[["sel_mtd"]][["MTD"]] - 1,
          boin_sel_mtd_res[["sel_mtd"]][["MTD"]]),
      boin_sel_mtd_res[["sel_mtd"]][["MTD"]] == 1 ~ 
        c(NA,
          NA, 
          boin_sel_mtd_res[["sel_mtd"]][["MTD"]]),
      # returns 99 if MTD is below lowest level
      # return vectors of all NA in this case
      boin_sel_mtd_res[["sel_mtd"]][["MTD"]] == 99 ~ c(NA, NA, NA)
      # BOIN never says the MTD is higher than the highest level - always will return the highest dose I believe
      # We could add something here to check the DLT probability and consider saying the MTD was not reached if the highest dose was listed but the DLT probability is below some threshold
    )
  
  return(suff_safe)
}