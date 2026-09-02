# 1. Sequential MSM data generation -------------------------------------------
#
# Arguments:
#   n                   Sample size (default 1000).
#   n_lines             Maximum number of treatment lines (default 3).
#   interval_length     Duration of each treatment interval in years
#                       (default 1).
#   interval_length_fn  A function to specify variable treatment duration
#                       across patients. If NULL uses interval_length.
#                       Example: interval_length_fn = function(n, k) {
#                                rlnorm(n, meanlog = log(1.5) - 0.2 * (k - 1),
#                                sdlog = 0.4)}
#   cumulative_effect  Logical, whether treatment effects are cumulative across
#                      lines
#
#   vars_means        Numeric vector of length P (MVN means for covariates).
#   vars_covs         P x P covariance matrix.
#   covariate_spec    Data frame with P rows. Required columns:
#                       name  - character, covariate name
#                       type  - "continuous" or "binary"
#                       beta  - numeric, log-HR for the outcome model
#                     Optional columns (defaults applied if missing):
#                       clip_lower         - lower bound for clipping (NA = none)
#                       clip_upper         - upper bound for clipping (NA = none)
#                       binary_threshold   - threshold for MVN -> binary (default 0.5)
#                       standardize_for_trt - logical, standardize for treatment
#                                            model? (default TRUE for continuous,
#                                            FALSE for binary)
#
#   treatment_names   Character vector of length J. First element is the
#                     reference treatment. Example: c("Imatinib", "Dasatinib",
#                     "Nilotinib").
#   beta_treatments   Numeric vector of length J-1. Log-HR for each
#                     non-reference treatment vs the reference. Order must
#                     match treatment_names[-1].
#   trt_coefs         (J-1) x (P+1) matrix. Columns = [intercept, cov_1, ...,
#                     cov_P]. Rows = each non-reference treatment vs reference.
#   trt_history_coefs (J-1) x (J-1) matrix. Effect of cumulative count of
#                     each prior non-reference treatment on current treatment
#                     log-odds.
#   confounding_strength  Scalar multiplier on covariate effects in the
#                     treatment model (not on intercept or history).
#
#   tv_conf           NULL or named list:
#                       L0_mean, L0_sd      - baseline L distribution
#                       alpha_intercept     - intercept in L_k model
#                       alpha_lag_L         - AR coefficient for L_{k-1}
#                       alpha_lag_trt       - numeric vector of length J-1,
#                                            effect of each previous non-ref
#                                            treatment on L_k
#                       sigma_L             - noise SD in L_k model
#                       gamma_L_on_trt      - scalar, effect of L_k on all
#                                            treatment log-odds
#                       beta_L              - effect of L_k on log-hazard
#
#   baseline_log_hazard  Numeric vector of length >= n_lines.
#   cens_rate            Rate of exponential censoring.
#   seed                 Integer seed.

simulate_sequential_msm <- function(
    n = 1000,
    n_lines = 3,
    interval_length = 1,
    interval_length_fn = NULL,
    cumulative_effect = FALSE,
    # Covariates
    vars_means,
    vars_covs,
    covariate_spec,
    # Treatments
    treatment_names,
    trt_coefs,
    trt_history_coefs,
    beta_treatments,
    confounding_strength = 1.0,
    # Time-varying confounder
    tv_conf = NULL,
    # Outcome
    baseline_log_hazard = rep(-3.5, n_lines),
    # Censoring
    cens_rate = 0.15,
    # Seed
    seed
) {

  set.seed(seed)

  # ---- Infer dimensions ----
  P <- length(vars_means) # number of covariates
  J <- length(treatment_names) # number of treatments
  ref_trt <- treatment_names[1] # reference treatment
  nonref_trts <- treatment_names[-1] # non-reference treatments

  # ---- Validate dimensions ----
  stopifnot(nrow(vars_covs) == P, ncol(vars_covs) == P)
  stopifnot(nrow(covariate_spec) == P)
  stopifnot(length(beta_treatments) == J - 1)
  stopifnot(nrow(trt_coefs) == J - 1, ncol(trt_coefs) == P + 1)
  stopifnot(nrow(trt_history_coefs) == J - 1, ncol(trt_history_coefs) == J - 1)
  stopifnot(length(baseline_log_hazard) >= n_lines)

  # ---- Fill in covariate_spec defaults ----
  # If clip_lower, clip_upper not specified, fill in NA;
  # If binary_threshold not specified, fill in 0.5;
  # If standardize_for_trt not specified, fill in standardize for all continuous
  if (!"clip_lower" %in% names(covariate_spec))
    covariate_spec$clip_lower <- NA_real_
  if (!"clip_upper" %in% names(covariate_spec))
    covariate_spec$clip_upper <- NA_real_
  if (!"binary_threshold" %in% names(covariate_spec))
    covariate_spec$binary_threshold <- ifelse(
      covariate_spec$type == "binary", 0.5, NA_real_
    )
  if (!"standardize_for_trt" %in% names(covariate_spec))
    covariate_spec$standardize_for_trt <- covariate_spec$type == "continuous"

  # Indicator for whether to use time-varying confounder and variable interval lengths
  use_tv <- !is.null(tv_conf)
  use_var_int <- !is.null(interval_length_fn)

  if (use_tv) {
    stopifnot(length(tv_conf$alpha_lag_trt) == J - 1)
  }

  # ---- Step 1: Baseline covariates from MVN ----
  # Draw the covariates from the multivariate normal distribution
  mvn_draw <- mvtnorm::rmvnorm(n, mean = vars_means, sigma = vars_covs)

  # Set up empty matrices for covariates and treatment model covariates
  covariates <- matrix(NA_real_, nrow = n, ncol = P)
  colnames(covariates) <- covariate_spec$name
  cov_trt <- matrix(NA_real_, nrow = n, ncol = P)

  for (p in seq_len(P)) {
    spec <- covariate_spec[p, ]
    raw <- mvn_draw[, p]

    # For continuous covariates, apply clipping and standardization if specified
    if (spec$type == "continuous") {
      if (!is.na(spec$clip_lower)) raw <- pmax(raw, spec$clip_lower)
      if (!is.na(spec$clip_upper)) raw <- pmin(raw, spec$clip_upper)
      covariates[, p] <- raw
      if (isTRUE(spec$standardize_for_trt)) {
        cov_trt[, p] <- (raw - vars_means[p]) / sqrt(vars_covs[p, p])
      } else {
        cov_trt[, p] <- raw
      }
    } else { # For binary covariates dichotomize based on threshold
      covariates[, p] <- as.integer(raw > spec$binary_threshold)
      cov_trt[, p] <- covariates[, p]
    }
  }

  # Treatment model design matrix: intercept + P covariates
  X_trt <- cbind(1, cov_trt)

  # ---- Step 2: Censoring time ----
  # Censoring times are from an exponential distribution with specified rate
  cens_time <- rexp(n, rate = cens_rate)

  # ---- Step 3: Initialize storage ----
  # Storage for treatment assignment
  treatment <- matrix(NA_character_, nrow = n, ncol = n_lines)
  colnames(treatment) <- paste0("treatment_", seq_len(n_lines))

  # Storage for interval lengths
  int_lengths <- matrix(NA_real_, nrow = n, ncol = n_lines)

  # Storage for time-varying confounder if applicable
  L <- NULL
  if (use_tv) {
    L <- matrix(NA_real_, nrow = n, ncol = n_lines)
  }

  # Storage for OS time and event indicator;
  # at_risk starts as TRUE for everyone and will change to FALSE
  # when they have an event or are censored
  os_yrs   <- rep(NA_real_, n)
  os_event <- rep(NA_integer_, n)
  at_risk  <- rep(TRUE, n)

  # ---- Step 4: Sequential simulation ----
  for (k in seq_len(n_lines)) {
    idx <- which(at_risk) # IDs for patients still at risk at the start of line k
    if (length(idx) == 0) break # Stop if no one is at risk
    nk <- length(idx) # the number of patients still at risk

    # Interval lengths
    if (use_var_int) {
      int_len_k <- interval_length_fn(nk, k)
      stopifnot(length(int_len_k) == nk, all(int_len_k > 0))
    } else {
      int_len_k <- rep(interval_length, nk)
    }
    int_lengths[idx, k] <- int_len_k

    # Interval start
    if (k == 1) {
      # all 0 at first line
      interval_start <- rep(0, nk)
    } else {
      # otherwise sum the intervals to get interval start for subsequent lines
      interval_start <- rowSums(int_lengths[idx, 1:(k - 1), drop = FALSE])
    }

    # -- (a) Time-varying confounder --
    if (use_tv) {
      if (k == 1) {
        # At line 1, draw from baseline normal distribution
        L[idx, k] <- rnorm(nk, mean = tv_conf$L0_mean, sd = tv_conf$L0_sd)
      } else {
        # Otherwise, draw from the AR(1) model with previous L and previous treatment
        prev_trt <- treatment[idx, k - 1]
        L_new <- tv_conf$alpha_intercept +
          tv_conf$alpha_lag_L * L[idx, k - 1]
        for (j in seq_along(nonref_trts)) {
          L_new <- L_new +
            tv_conf$alpha_lag_trt[j] * as.integer(prev_trt == nonref_trts[j])
        }
        L[idx, k] <- L_new + rnorm(nk, 0, tv_conf$sigma_L)
      }
    }

    # -- (b) Treatment assignment --
    # Cumulative counts of each non-reference treatment in prior lines
    n_prev <- matrix(0L, nrow = nk, ncol = J - 1)
    if (k > 1) {
      prev_mat <- treatment[idx, 1:(k - 1), drop = FALSE]
      for (j in seq_along(nonref_trts)) {
        n_prev[, j] <- rowSums(prev_mat == nonref_trts[j], na.rm = TRUE)
      }
    }

    X_k <- X_trt[idx, , drop = FALSE] # limit to those still at risk
    base_effects <- X_k[, -1, drop = FALSE] # remove intercept column

    # Linear predictors for each non-reference treatment
    eta <- matrix(0, nrow = nk, ncol = J - 1)
    for (j in seq_len(J - 1)) {
      eta[, j] <- X_k[, 1] * trt_coefs[j, 1] +
        as.numeric(base_effects %*% trt_coefs[j, -1]) * confounding_strength +
        as.numeric(n_prev %*% trt_history_coefs[j, ])
      if (use_tv) {
        eta[, j] <- eta[, j] + tv_conf$gamma_L_on_trt * L[idx, k]
      }
    }

    # Softmax probabilities
    exp_eta <- exp(eta) # exponentiate the linear predictor
    denom_sm <- 1 + rowSums(exp_eta) # denominator of the probability
    # Combine the reference category probability with the
    # non-reference probabilities, back to nk x k dimension
    probs <- cbind(1 / denom_sm, exp_eta / denom_sm)
    # This is the matrix of probabilities of getting each treatment
    # for each patient at this line
    colnames(probs) <- treatment_names

    # Multinomial draw
    # For each row, sample one treatment according to the probabilities in that row
    trt_draw <- apply(probs, 1, function(p) {
      sample(treatment_names, 1, prob = p)
    })
    treatment[idx, k] <- trt_draw

    # -- (c) Outcome: piecewise-exponential hazard --
    # Create dummy variable matrix for the non-reference treatments based on
    # preceding treatment assignment for this line
    trt_indicators <- matrix(0L, nrow = nk, ncol = J - 1)
    for (j in seq_along(nonref_trts)) {
      trt_indicators[, j] <- as.integer(trt_draw == nonref_trts[j])
    }

    # Will be same as above at first line, otherwise treatment effects accumulate
    if (cumulative_effect) {
      trt_effect <- n_prev + trt_indicators
    } else {
      trt_effect <- trt_indicators
    }

    # Linear predictor for the outcome, based on covariates and treatments
    linpred <- baseline_log_hazard[k] +
      as.numeric(covariates[idx, , drop = FALSE] %*% covariate_spec$beta) +
      as.numeric(trt_effect %*% beta_treatments)

    # Add in the effect of time-varying confounder if applicable
    if (use_tv) {
      linpred <- linpred + tv_conf$beta_L * L[idx, k]
    }

    hazard_k <- exp(linpred) # hazard of having the event
    T_k <- rexp(nk, rate = hazard_k) # event time from exponential distribution
    # How much censoring time is left from start of current line
    remaining_cens <- cens_time[idx] - interval_start

    # -- (d) Interval resolution --
    # Indicator for if the event happened before the end of the treatment
    # interval and before the censoring time
    event_here <- (T_k < int_len_k) & (T_k <= remaining_cens) &
      (remaining_cens > 0)
    # Indicator for if the patient was censored during this interval
    cens_here  <- (remaining_cens < int_len_k) &
      (remaining_cens < T_k) &
      (remaining_cens > 0) & (!event_here)

    # IDs of patients who had the event in this interval
    ev_idx <- idx[event_here]
    # event time of patients who had the event in this interval
    os_yrs[ev_idx]   <- interval_start[event_here] + T_k[event_here]
    # event indicator of patients who had the event in this interval
    os_event[ev_idx]  <- 1L
    # change at_risk to FALSE for patients who had the event in this interval
    at_risk[ev_idx]   <- FALSE

    # IDs of patients who were censored in this interval
    cn_idx <- idx[cens_here]
    # censoring time of patients who were censored in this interval
    os_yrs[cn_idx]   <- cens_time[cn_idx]
    # censoring indicator of patients who were censored in this interval
    os_event[cn_idx]  <- 0L
    # change at_risk to FALSE for patients who were censored in this interval
    at_risk[cn_idx]   <- FALSE
  }

  # Patients surviving all intervals
  still_alive <- which(at_risk)
  # Patients surviving all intervals are censored at the first of their
  # censoring time or the end of their last treatment
  max_fu <- rowSums(int_lengths[still_alive, , drop = FALSE], na.rm = TRUE)
  os_yrs[still_alive]  <- pmin(cens_time[still_alive], max_fu)
  os_event[still_alive] <- 0L

  # ---- Step 5: Build wide-format data ----
  data_wide <- data.frame(patient_id = seq_len(n), stringsAsFactors = FALSE)
  for (p in seq_len(P)) {
    data_wide[[covariate_spec$name[p]]] <- covariates[, p]
  }
  for (k in seq_len(n_lines)) {
    data_wide[[paste0("treatment_", k)]] <- factor(
      treatment[, k], levels = treatment_names
    )
    data_wide[[paste0("interval_length_", k)]] <- int_lengths[, k]
  }
  if (use_tv) {
    for (k in seq_len(n_lines)) {
      data_wide[[paste0("L_", k)]] <- L[, k]
    }
  }
  data_wide$os_yrs  <- os_yrs
  data_wide$os_event <- os_event
  data_wide$n_lines_received <- rowSums(!is.na(treatment))

  # ---- Step 6: Build long-format (counting-process) data ----
  long_list <- vector("list", n_lines)

  for (k in seq_len(n_lines)) {
    has_trt <- !is.na(data_wide[[paste0("treatment_", k)]])
    sub <- data_wide[has_trt, ]
    ns <- nrow(sub)

    if (k == 1) {
      int_start <- rep(0, ns)
    } else {
      int_cols <- paste0("interval_length_", seq_len(k - 1))
      int_start <- rowSums(sub[, int_cols, drop = FALSE])
    }
    int_end <- int_start + sub[[paste0("interval_length_", k)]]

    # Cumulative treatment counts through this line
    cum_trt <- matrix(0L, nrow = ns, ncol = J - 1)
    cum_names <- paste0("cum_", nonref_trts)
    colnames(cum_trt) <- cum_names
    for (j in seq_along(nonref_trts)) {
      for (jj in seq_len(k)) {
        cum_trt[, j] <- cum_trt[, j] +
          as.integer(sub[[paste0("treatment_", jj)]] == nonref_trts[j])
      }
    }

    row_k <- data.frame(
      patient_id = sub$patient_id,
      line       = k,
      tstart     = int_start,
      tstop      = pmin(int_end, sub$os_yrs),
      event      = as.integer(sub$os_event == 1L &
                                sub$os_yrs >  int_start &
                                sub$os_yrs <= int_end),
      treatment  = factor(sub[[paste0("treatment_", k)]],
                          levels = treatment_names),
      stringsAsFactors = FALSE
    )

    for (j in seq_along(cum_names)) {
      row_k[[cum_names[j]]] <- cum_trt[, j]
    }
    for (p in seq_len(P)) {
      row_k[[covariate_spec$name[p]]] <- sub[[covariate_spec$name[p]]]
    }
    if (use_tv) {
      row_k$L <- sub[[paste0("L_", k)]]
    }

    long_list[[k]] <- row_k
  }

  data_long <- do.call(rbind, long_list)
  data_long <- data_long[order(data_long$patient_id, data_long$line), ]
  rownames(data_long) <- NULL

  # ---- Step 7: Return ----
  true_params <- list(
    beta_treatments     = beta_treatments,
    beta_covariates     = setNames(covariate_spec$beta, covariate_spec$name),
    covariate_names     = covariate_spec$name,
    treatment_names     = treatment_names,
    nonref_treatments   = nonref_trts,
    cumulative_effect   = cumulative_effect,
    baseline_log_hazard = baseline_log_hazard[seq_len(n_lines)],
    confounding_strength = confounding_strength,
    interval_length     = interval_length,
    n_lines             = n_lines,
    cens_rate           = cens_rate
  )
  if (use_tv) true_params$tv_conf <- tv_conf

  list(
    data_wide   = data_wide,
    data_long   = data_long,
    true_params = true_params
  )
}



################################################################################
# Example dataset generation using default simulation parameters ----
################################################################################

# ############################
# Setup the function inputs
# ############################

# Mean vector for baseline covariates
vars_means1 <- c(mean_age = 55.0492508, mean_female = 0.4243323,
                mean_comorbid = 0.5341246, mean_smoke = 0.4658754)

# Covariance matrix for baseline covariates
vars_covs1 <- matrix(
  c(269.06477015, -0.05578908, 4.15481897,  1.00848242,
    -0.05578908,  0.24500141, -0.01303518, -0.02565529,
    4.15481897, -0.01303518, 0.24957609, 0.03018581,
    1.00848242, -0.02565529, 0.03018581, 0.24957609),
  nrow = 4, ncol = 4, byrow = TRUE,
  dimnames = list(c("first_age_dx", "female", "any_comorbid", "ever_smoke"),
                  c("first_age_dx", "female", "any_comorbid", "ever_smoke")))

# Covariate specification: P = 4
covariate_spec1 <- data.frame(
  name     = c("first_age_dx", "female", "any_comorbid", "ever_smoke"),
  type     = c("continuous",   "binary", "binary",       "binary"),
  beta     = c(log(1.03),      log(0.81), log(1.98),     log(1.21)),
  clip_lower = c(18, NA, NA, NA),
  clip_upper = c(90, NA, NA, NA),
  stringsAsFactors = FALSE
)

# Treatment names: J = 3 (first is reference)
treatment_names1 <- c("Imatinib", "Dasatinib", "Nilotinib")

# Treatment model coefficients: (J-1) x (P+1)
trt_coefs1 <- matrix(c(
  -0.3,  0.15, -0.10,  0.30,  0.10,   # Dasatinib vs Imatinib
  -0.8,  0.10,  0.20,  0.15, -0.20    # Nilotinib vs Imatinib
), nrow = 2, byrow = TRUE)

# Treatment history coefficients: (J-1) x (J-1)
trt_history_coefs1 <- matrix(c(
  0.8, -0.3,   # effect of prior D/N counts on current D log-odds
  -0.3,  0.7    # effect of prior D/N counts on current N log-odds
), nrow = 2, byrow = TRUE)

# True causal treatment effects: J-1 = 2
beta_treatments1 <- c(log(1.71), log(0.81))

# Number of treatment lines
n_lines <- 3



# #####################
# Execute the function
# #####################
dat1 <- simulate_sequential_msm(
  n = 1000,
  n_lines = n_lines,
  interval_length = 1,
  interval_length_fn = NULL,
  cumulative_effect = FALSE,
  # Covariates
  vars_means = vars_means1,
  vars_covs = vars_covs1,
  covariate_spec = covariate_spec1,
  # Treatments
  treatment_names = treatment_names1,
  trt_coefs = trt_coefs1,
  trt_history_coefs = trt_history_coefs1,
  beta_treatments = beta_treatments1,
  confounding_strength = 1.0,
  # Time-varying confounder
  tv_conf = NULL,
  # Outcome
  baseline_log_hazard = rep(-3.5, n_lines),
  # Censoring
  cens_rate = 0.15,
  # Seed
  seed = 1098
)
