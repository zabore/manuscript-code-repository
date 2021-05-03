# Function to create a single simulated dataset -------------------------------
# hr12 = HR compring high risk to low risk group for the transition from state 1 (alive and disease free) to state 2 (recurrence)
# hr13 = HR comparing high risk to low risk group for the transition from state 1 (alive and disease free) to state 3 (death)
# hr23 = HR comparing high risk to low risk group for the transition from state 2 (recurrence) to state 3 (death)
# mos_z0 = frequency of visits in the low risk group, in months
# mos_z1 = frequency of visits in the high risk group, in months
# n = total sample size, will be equally divided between the low risk and high risk group so must be an even number
# dv = is there an extra visit dependent on recurrence? If so
makesimdat <- 
  function(hr12 = 1.5, 
           hr13 = 1.5, 
           hr23 = 1, 
           mos_z0 = 6,
           mos_z1 = 6,
           n = 100,
           dv = FALSE) {
    
    if(dv == FALSE) {
      grp0 <- 
        simdat(
          n = n / 2,
          scale12 = 1 / .0008,
          scale13 = 1 / .0002,
          scale23 = 1 / .0016,
          visit.schedule = 30.4 * seq(mos_z0, 72, mos_z0),
          dependent.visit = NULL, 
          vital.lfu = c(0, 30.4 * 72)
        )
      
      grp1 <- 
        simdat(
          n = n / 2,
          scale12 = (1 / .0008) / hr12,
          scale13 = (1 / .0002) / hr13,
          scale23 = (1 / .0016) / hr23,
          visit.schedule = 30.4 * seq(mos_z1, 72, mos_z1),
          dependent.visit = NULL, 
          vital.lfu = c(0, 30.4 * 72)
        )
    } else if(dv == TRUE) {
      grp0 <- 
        simdat(
          n = n / 2,
          scale12 = 1 / .0008,
          scale13 = 1 / .0002,
          scale23 = 1 / .0016,
          visit.schedule = 30.4 * seq(mos_z0, 72, mos_z0),
          dependent.visit = c(30.4, 30.4/3), 
          vital.lfu = c(0, 30.4 * 72)
        )
      
      grp1 <- 
        simdat(
          n = n / 2,
          scale12 = (1 / .0008) / hr12,
          scale13 = (1 / .0002) / hr13,
          scale23 = (1 / .0016) / hr23,
          visit.schedule = 30.4 * seq(mos_z1, 72, mos_z1),
          dependent.visit = c(30.4, 30.4/3), 
          vital.lfu = c(0, 30.4 * 72)
        )
    }
    
    full_join(
      grp0 %>% mutate(z = 0), 
      grp1 %>% mutate(z = 1)
    )
  }


# Example usage ---------------------------------------------------------------
# Uncomment below to try
# df <- makesimdat(hr12 = 1.8,
#                  hr13 = 1.8,
#                  hr23 = 1,
#                  mos_z0 = 6,
#                  mos_z1 = 2,
#                  n = 500,
#                  dv = FALSE)
# head(df)