library(purrr)
library(tibble)


# Function to generate the simulated data for the stratified control arm design
sim_dat1 <- function(p, n) {
  if (length(p) == 2) {
    if (length(n) == 2 & is.matrix(n) == FALSE) {
      n <- matrix(n, nrow = 1)
    }
    y0 <- rbinom(n = 1, size = n[1, 1], prob = p[1])
    y1 <- rbinom(n = 1, size = n[1, 2], prob = p[2])
    if (length(n) > 2) {
      for (i in seq_len(nrow(n))[-1]) {
        y0 <- c(y0, y0[length(y0)] + rbinom(
          n = 1,
          size = n[i, 1] - n[i - 1, 1],
          prob = p[1]
        ))
        y1 <- c(y1, y1[length(y1)] + rbinom(
          n = 1,
          size = n[i, 2] - n[i - 1, 2],
          prob = p[2]
        ))
      }
    }
    return(tibble(
      n0 = n[, 1],
      n1 = n[, 2],
      y0 = y0,
      y1 = y1
    ))
  } else if (length(p) == 1) {
    y1 <- rbinom(n = 1, size = n[1], prob = p)
    if (length(n) > 1) {
      for (i in seq_along(n)[-1]) {
        y1 <- c(y1, y1[length(y1)] +
                  rbinom(n = 1, size = n[i] - n[i - 1], prob = p))
      }
    }
    return(tibble(n1 = n, y1 = y1))
  }
}


# simulation set up
set.seed(20220405)
n <- cbind(seq(10, 50, 10), seq(10, 50, 10))
nsim <- 1000


# generate the null datasets
s_sim_dat_trt0_null <- map(1:nsim, ~sim_dat1(p = c(0.1, 0.1), n = n))
s_sim_dat_trt1_null <- map(1:nsim, ~sim_dat1(p = c(0.1, 0.1), n = n))
s_sim_dat_trt2_null <- map(1:nsim, ~sim_dat1(p = c(0.1, 0.1), n = n))


# generate the alternative datasets
s_sim_dat_trt0_alt <- map(1:nsim, ~sim_dat1(p = c(0.1, 0.1), n = n))
s_sim_dat_trt1_alt <- map(1:nsim, ~sim_dat1(p = c(0.1, 0.2), n = n))
s_sim_dat_trt2_alt <- map(1:nsim, ~sim_dat1(p = c(0.1, 0.3), n = n))


# save all datasets
save(s_sim_dat_trt0_null, s_sim_dat_trt1_null, s_sim_dat_trt2_null, 
     s_sim_dat_trt0_alt, s_sim_dat_trt1_alt, s_sim_dat_trt2_alt,
     file = here::here("code", "stratified-control-group", "s-sim-dat.rda"))