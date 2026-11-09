# Source the function to generate a simulated dataset ----
# NOTE: Alter the file path to match where you have saved the file!
source(here::here("code", "msm-sim-study-generate-data-function.R"))

# Simulate a single dataset, broadly in line with the case study ----
dat <- simulate_sequential_msm_switch(
    n = 350,
    n_lines = 3,
    interval_length_fn = function(n, k) rlnorm(n, log(1.5), 0.5),
    cumulative_effect = FALSE,
    # Covariates
    vars_means = c(first_age_dx = 55,
                   spleen_size_through_imaging = 16,
                   any_cvd = 0.4,
                   ever_smoke = 0.5),
    vars_covs = matrix(c(264, -23, 3.3, 1,
                         -23, 65, .28, 1,
                         3.3, .28, .24, .04,
                         1, 1, .04, .25),
                       nrow = 4, ncol = 4, byrow = TRUE,
                       dimnames = list(c("first_age_dx",
                                         "spleen_size_through_imaging",
                                         "any_cvd",
                                         "ever_smoke"),
                                       c("first_age_dx",
                                         "spleen_size_through_imaging",
                                         "any_cvd",
                                         "ever_smoke"))),
    covariate_spec = data.frame(
      name = c("first_age_dx", "spleen_size_through_imaging", "any_cvd",
               "ever_smoke"),
      type = c("continuous", "continuous", "binary",  "binary"),
      beta = c(0.067, 0.053, 0.017, 0.47),
      clip_lower = c(18, 1, NA, NA),
      clip_upper = c(90, 40, NA, NA),
      standardize_for_trt = FALSE,
      stringsAsFactors = FALSE
    ),
    # Treatments
    treatment_names = c("Imatinib", "Dasatinib", "Nilotinib", "Others"),
    trt_coefs = matrix(c(-0.107, -0.015, 0.025, 0.158, -0.204,
                         -0.386, -0.001, -0.05, -0.544, 0.146,
                         -2.22, -0.03, -0.04, 1.88, -0.068),
                       nrow = 3, ncol = 5, byrow = TRUE),
    trt_history_coefs = matrix(0, 3, 3),                # replaced by trt_prev_coefs below
    trt_prev_coefs = list(                              # columns = previous treatment (Imatinib, Dasatinib, Nilotinib, Others)
      matrix(c(2.404, -0.276, 1.082, 0.794,             # line 2; rows = Dasatinib, Nilotinib, Others
               2.918,  1.692,  -Inf, 1.792,
               3.320,  3.549, 3.384, 4.013), 3, 4, byrow = TRUE),
      matrix(c(2.075, 0.571,  1.264, 0.000,             # line 3
               2.079, 0.901, -0.511, 0.410,
               4.706, 3.997,  4.405, 3.730), 3, 4, byrow = TRUE)),
    trt_prev_ref = list(c(0, 0, 0, 0), c(0, 0, 0, -Inf)), # line 3: Imatinib impossible after Others
    beta_treatments = c(0.037, -0.511, -0.127),
    confounding_strength = 1.0,
    # Time-varying confounder
    tv_conf = list(
      L0_mean = 0,
      L0_sd = 0.23,
      alpha_intercept = -0.82,
      alpha_lag_L = 0.5,
      alpha_lag_trt = c(-0.6, -0.9, -0.5),
      sigma_L = 1.38,
      gamma_L_on_trt = 0.3,
      beta_L = log(1.3),
      L_center = c(0, -1.1, -1.9)
    ),
    switch_spec = list(intercept = c(0.5, 0.3), beta_L = 0.6),
    # Outcome
    baseline_log_hazard = rep(-8.18, 3),
    # Censoring
    cens_rate = 0.15,
    cens_fn = function(n) pmin(rexp(n, 0.02), runif(n, 3, 15)),   # loss to follow-up + data lock 3-15 y after entry
    # Seed
    seed = 1185
    )

# Save the dataset
# NOTE: Alter the file path to match where you want to save the data file!
save(dat, file = here::here("data", "msm-sim-study-single-simulated-dataset.rda"))
