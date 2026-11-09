# 1. Sequential MSM data generation -------------------------------------------
# Generates one simulated dataset of multi-line treatment, a time-varying
# confounder and a survival outcome. This is a self-contained variant of
# simulate_sequential_msm() from the simulation study (simulation-functions-v8.R),
# written for the case-study example data. It differs from that function in three
# ways; everything else (covariates, treatment model, outcome model, output
# format) is unchanged.
#
#   1. SWITCHING. In the original, every patient who is alive and uncensored at the
#      end of line k automatically starts line k+1, so the share reaching later lines
#      is fixed by the hazard, censoring and interval length. Here (when switch_spec
#      is supplied) such a patient switches with probability
#           expit( switch_spec$intercept[k] + switch_spec$beta_L * (L_{k+1} - L_center[k+1]) )
#      and otherwise stays on the current treatment until death or censoring.
#      L_{k+1} (e.g. BCR-ABL1 measured during line k) is generated before the
#      switching decision, so it exists even for patients who do not switch.
#
#   2. PREVIOUS-TREATMENT EFFECTS (trt_prev_coefs, trt_prev_ref). Treatment at line
#      k >= 2 can depend on the treatment received at line k-1, INCLUDING the
#      reference treatment. The original only has counts of prior non-reference
#      treatments, so a patient on the reference treatment could not be made to
#      switch away from it. Impossible sequences are obtained with -Inf.
#
#   3. CENTERING OF L in the treatment and switching models (tv_conf$L_center), so
#      that the treatment mix at each line does not depend on the level of L.
#
# The trt_sampler hook of the original (used by compute_msm_truth) is not included.
#
# Arguments that are new or changed (all others as in simulate_sequential_msm()):
#   trt_prev_coefs  NULL, a (J-1) x J matrix, or a list of n_lines-1 such matrices
#                   (one per line 2..n_lines). Entry [j, p] is added to the log-odds of
#                   non-reference treatment j when the previous treatment is
#                   treatment_names[p]. Use -Inf for an impossible transition.
#   trt_prev_ref    NULL, a length-J vector, or a list of n_lines-1 such vectors.
#                   Entry [p] is the log-weight of the REFERENCE treatment when the
#                   previous treatment is treatment_names[p] (default 0). Use -Inf to make
#                   the reference treatment impossible after treatment p.
#   tv_conf         as before, plus optional L_center (length n_lines): the reference
#                   level of L subtracted before it enters the treatment and switching
#                   models. The hazard still uses L itself.
#   switch_spec     NULL (original behaviour) or list(intercept = numeric(n_lines - 1),
#                   beta_L = scalar, default 0).
#
# The last line is terminal: patients on it stay until death or censoring. The
# interval (interval_length / interval_length_fn) is the time on a line before a
# switch can occur.
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


simulate_sequential_msm_switch <- function(
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
    trt_prev_coefs = NULL,
    trt_prev_ref = NULL,
    beta_treatments,
    confounding_strength = 1.0,
    # Time-varying confounder
    tv_conf = NULL,
    # Switching between lines
    switch_spec = NULL,
    # Outcome
    baseline_log_hazard = rep(-3.5, n_lines),
    # Censoring
    cens_rate = 0.15,
    cens_fn = NULL,
    seed
) {

  set.seed(seed)

  P <- length(vars_means)
  J <- length(treatment_names)
  nonref_trts <- treatment_names[-1]

  stopifnot(nrow(vars_covs) == P, ncol(vars_covs) == P)
  stopifnot(nrow(covariate_spec) == P)
  stopifnot(length(beta_treatments) == J - 1)
  stopifnot(nrow(trt_coefs) == J - 1, ncol(trt_coefs) == P + 1)
  stopifnot(nrow(trt_history_coefs) == J - 1, ncol(trt_history_coefs) == J - 1)
  stopifnot(length(baseline_log_hazard) >= n_lines)

  # ---- covariate_spec defaults ----
  if (!"clip_lower" %in% names(covariate_spec)) covariate_spec$clip_lower <- NA_real_
  if (!"clip_upper" %in% names(covariate_spec)) covariate_spec$clip_upper <- NA_real_
  if (!"binary_threshold" %in% names(covariate_spec))
    covariate_spec$binary_threshold <- ifelse(covariate_spec$type == "binary", 0.5, NA_real_)
  if (!"standardize_for_trt" %in% names(covariate_spec))
    covariate_spec$standardize_for_trt <- covariate_spec$type == "continuous"

  use_tv     <- !is.null(tv_conf)
  use_var_int <- !is.null(interval_length_fn)
  use_switch <- !is.null(switch_spec)
  if (use_tv) stopifnot(length(tv_conf$alpha_lag_trt) == J - 1)
  L_center <- if (use_tv && !is.null(tv_conf$L_center)) tv_conf$L_center else rep(0, n_lines)
  stopifnot(length(L_center) == n_lines)
  if (use_switch) {
    stopifnot(length(switch_spec$intercept) == n_lines - 1)
    if (is.null(switch_spec$beta_L)) switch_spec$beta_L <- 0
  }

  # previous-treatment effects, one entry per line 2..n_lines
  as_line_list <- function(x, default) {
    if (n_lines < 2) return(list())
    if (is.null(x)) return(rep(list(default), n_lines - 1))
    if (is.list(x)) { stopifnot(length(x) == n_lines - 1); return(x) }
    rep(list(x), n_lines - 1)
  }
  prev_coefs <- as_line_list(trt_prev_coefs, matrix(0, J - 1, J))
  prev_ref   <- as_line_list(trt_prev_ref, rep(0, J))
  for (m in seq_along(prev_coefs)) {
    stopifnot(nrow(prev_coefs[[m]]) == J - 1, ncol(prev_coefs[[m]]) == J, length(prev_ref[[m]]) == J)
  }

  # ---- Step 1: baseline covariates ----
  mvn_draw   <- mvtnorm::rmvnorm(n, mean = vars_means, sigma = vars_covs)
  covariates <- matrix(NA_real_, nrow = n, ncol = P)
  colnames(covariates) <- covariate_spec$name
  cov_trt <- matrix(NA_real_, nrow = n, ncol = P)

  for (p in seq_len(P)) {
    spec <- covariate_spec[p, ]
    raw  <- mvn_draw[, p]
    if (spec$type == "continuous") {
      if (!is.na(spec$clip_lower)) raw <- pmax(raw, spec$clip_lower)
      if (!is.na(spec$clip_upper)) raw <- pmin(raw, spec$clip_upper)
      covariates[, p] <- raw
      cov_trt[, p] <- if (isTRUE(spec$standardize_for_trt))
        (raw - vars_means[p]) / sqrt(vars_covs[p, p]) else raw
    } else {
      covariates[, p] <- as.integer(raw > spec$binary_threshold)
      cov_trt[, p]    <- covariates[, p]
    }
  }
  X_trt <- cbind(1, cov_trt)

  # ---- Step 2: censoring times ----
  cens_time <- if (is.null(cens_fn)) rexp(n, rate = cens_rate) else cens_fn(n)
  stopifnot(length(cens_time) == n, all(cens_time > 0))

  # ---- Step 3: storage ----
  treatment   <- matrix(NA_character_, nrow = n, ncol = n_lines)
  int_lengths <- matrix(NA_real_, nrow = n, ncol = n_lines)
  L <- if (use_tv) matrix(NA_real_, nrow = n, ncol = n_lines) else NULL
  os_yrs   <- rep(NA_real_, n)
  os_event <- rep(NA_integer_, n)
  at_risk  <- rep(TRUE, n)

  # ---- Step 4: sequential simulation ----
  for (k in seq_len(n_lines)) {
    idx <- which(at_risk)
    if (length(idx) == 0) break
    nk  <- length(idx)

    # time on this line before a switch can occur; the last line (switch mode) has no end
    terminal <- use_switch && k == n_lines
    if (terminal) {
      int_len_k <- rep(Inf, nk)
    } else if (use_var_int) {
      int_len_k <- interval_length_fn(nk, k)
      stopifnot(length(int_len_k) == nk, all(int_len_k > 0))
    } else {
      int_len_k <- rep(interval_length, nk)
    }
    int_lengths[idx, k] <- int_len_k

    interval_start <- if (k == 1) rep(0, nk) else
      rowSums(int_lengths[idx, 1:(k - 1), drop = FALSE])

    # -- (a) time-varying confounder at the start of line k --
    # (in switch mode, L_k for k >= 2 was generated at the end of line k-1)
    if (use_tv && k == 1) {
      L[idx, 1] <- rnorm(nk, mean = tv_conf$L0_mean, sd = tv_conf$L0_sd)
    } else if (use_tv && !use_switch) {
      prev_trt <- treatment[idx, k - 1]
      L_new <- tv_conf$alpha_intercept + tv_conf$alpha_lag_L * L[idx, k - 1]
      for (j in seq_along(nonref_trts))
        L_new <- L_new + tv_conf$alpha_lag_trt[j] * as.integer(prev_trt == nonref_trts[j])
      L[idx, k] <- L_new + rnorm(nk, 0, tv_conf$sigma_L)
    }

    # -- (b) treatment assignment --
    n_prev <- matrix(0L, nrow = nk, ncol = J - 1)
    if (k > 1) {
      prev_mat <- treatment[idx, 1:(k - 1), drop = FALSE]
      for (j in seq_along(nonref_trts))
        n_prev[, j] <- rowSums(prev_mat == nonref_trts[j], na.rm = TRUE)
    }

    X_k <- X_trt[idx, , drop = FALSE]
    base_effects <- X_k[, -1, drop = FALSE]
    eta <- matrix(0, nrow = nk, ncol = J - 1)
    for (j in seq_len(J - 1)) {
      eta[, j] <- X_k[, 1] * trt_coefs[j, 1] +
        as.numeric(base_effects %*% trt_coefs[j, -1]) * confounding_strength +
        as.numeric(n_prev %*% trt_history_coefs[j, ])
      if (use_tv) eta[, j] <- eta[, j] + tv_conf$gamma_L_on_trt * (L[idx, k] - L_center[k])
    }

    ref_w <- rep(1, nk)                               # weight of the reference treatment
    if (k > 1) {
      p_idx <- match(treatment[idx, k - 1], treatment_names)
      for (j in seq_len(J - 1)) eta[, j] <- eta[, j] + prev_coefs[[k - 1]][j, p_idx]
      ref_w <- exp(prev_ref[[k - 1]][p_idx])
    }
    w_all <- cbind(ref_w, exp(eta))
    denom_sm <- rowSums(w_all)
    stopifnot(all(is.finite(denom_sm)), all(denom_sm > 0))
    probs <- w_all / denom_sm
    colnames(probs) <- treatment_names
    trt_draw <- apply(probs, 1, function(p) sample(treatment_names, 1, prob = p))
    treatment[idx, k] <- trt_draw

    # -- (b2) switch mode: the time-varying confounder measured DURING line k --
    # Generated for every patient who enters line k (not only those who survive the
    # interval or switch), so that, as in a cohort with complete data, a value exists
    # for everyone who reaches the line. It depends on L_k and on the treatment now
    # being received, and drives the switching decision below.
    if (use_switch && use_tv && !terminal) {
      L_next <- tv_conf$alpha_intercept + tv_conf$alpha_lag_L * L[idx, k]
      for (j in seq_along(nonref_trts))
        L_next <- L_next + tv_conf$alpha_lag_trt[j] * as.integer(trt_draw == nonref_trts[j])
      L[idx, k + 1] <- L_next + rnorm(nk, 0, tv_conf$sigma_L)
    }

    # -- (c) outcome: piecewise-exponential hazard --
    trt_indicators <- matrix(0L, nrow = nk, ncol = J - 1)
    for (j in seq_along(nonref_trts))
      trt_indicators[, j] <- as.integer(trt_draw == nonref_trts[j])
    trt_effect <- if (cumulative_effect) n_prev + trt_indicators else trt_indicators

    linpred <- baseline_log_hazard[k] +
      as.numeric(covariates[idx, , drop = FALSE] %*% covariate_spec$beta) +
      as.numeric(trt_effect %*% beta_treatments)
    if (use_tv) linpred <- linpred + tv_conf$beta_L * L[idx, k]

    T_k <- rexp(nk, rate = exp(linpred))
    remaining_cens <- cens_time[idx] - interval_start

    # -- (d) resolution within the interval --
    event_here <- (T_k < int_len_k) & (T_k <= remaining_cens) & (remaining_cens > 0)
    cens_here  <- (remaining_cens < int_len_k) & (remaining_cens < T_k) &
      (remaining_cens > 0) & (!event_here)

    ev_idx <- idx[event_here]
    os_yrs[ev_idx]   <- interval_start[event_here] + T_k[event_here]
    os_event[ev_idx] <- 1L
    at_risk[ev_idx]  <- FALSE

    cn_idx <- idx[cens_here]
    os_yrs[cn_idx]   <- cens_time[cn_idx]
    os_event[cn_idx] <- 0L
    at_risk[cn_idx]  <- FALSE

    # -- (e) switch mode: survivors of the interval switch or stay --
    if (use_switch) {
      surv <- !event_here & !cens_here
      if (!terminal && any(surv)) {
        s_pos <- which(surv)                       # positions within idx
        s_idx <- idx[s_pos]

        # L_{k+1} was generated in step (b2) above
        lp_sw <- switch_spec$intercept[k]
        if (use_tv) lp_sw <- lp_sw + switch_spec$beta_L * (L[s_idx, k + 1] - L_center[k + 1])
        goes  <- runif(length(s_idx)) < plogis(lp_sw)

        stay_pos <- s_pos[!goes]
        if (length(stay_pos) > 0) {
          st_idx   <- idx[stay_pos]
          st_event <- T_k[stay_pos] <= remaining_cens[stay_pos]
          os_yrs[st_idx]   <- ifelse(st_event, interval_start[stay_pos] + T_k[stay_pos], cens_time[st_idx])
          os_event[st_idx] <- as.integer(st_event)
          at_risk[st_idx]  <- FALSE
          int_lengths[st_idx, k] <- os_yrs[st_idx] - interval_start[stay_pos]
        }
      }
      if (terminal) {
        # everyone left on the last line is resolved above; record the observed duration
        int_lengths[idx, k] <- os_yrs[idx] - interval_start
      }
    }
  }

  # patients still at risk (original behaviour only): censored at the first of their
  # censoring time and the end of their last interval
  still_alive <- which(at_risk)
  if (length(still_alive) > 0) {
    max_fu <- rowSums(int_lengths[still_alive, , drop = FALSE], na.rm = TRUE)
    os_yrs[still_alive]   <- pmin(cens_time[still_alive], max_fu)
    os_event[still_alive] <- 0L
  }

  # ---- Step 5: wide-format data ----
  data_wide <- data.frame(patient_id = seq_len(n), stringsAsFactors = FALSE)
  for (p in seq_len(P)) data_wide[[covariate_spec$name[p]]] <- covariates[, p]
  for (k in seq_len(n_lines)) {
    data_wide[[paste0("treatment_", k)]] <- factor(treatment[, k], levels = treatment_names)
    data_wide[[paste0("interval_length_", k)]] <- int_lengths[, k]
  }
  if (use_tv) for (k in seq_len(n_lines)) data_wide[[paste0("L_", k)]] <- L[, k]
  data_wide$os_yrs   <- os_yrs
  data_wide$os_event <- os_event
  data_wide$n_lines_received <- rowSums(!is.na(treatment))

  # ---- Step 6: long-format (counting-process) data ----
  long_list <- vector("list", n_lines)
  for (k in seq_len(n_lines)) {
    has_trt <- !is.na(data_wide[[paste0("treatment_", k)]])
    sub <- data_wide[has_trt, ]
    ns  <- nrow(sub)

    int_start <- if (k == 1) rep(0, ns) else
      rowSums(sub[, paste0("interval_length_", seq_len(k - 1)), drop = FALSE])
    int_end <- int_start + sub[[paste0("interval_length_", k)]]

    cum_trt <- matrix(0L, nrow = ns, ncol = J - 1)
    cum_names <- paste0("cum_", nonref_trts)
    colnames(cum_trt) <- cum_names
    for (j in seq_along(nonref_trts))
      for (jj in seq_len(k))
        cum_trt[, j] <- cum_trt[, j] + as.integer(sub[[paste0("treatment_", jj)]] == nonref_trts[j])

    row_k <- data.frame(
      patient_id = sub$patient_id,
      line       = k,
      tstart     = int_start,
      tstop      = pmin(int_end, sub$os_yrs),
      event      = as.integer(sub$os_event == 1L & sub$os_yrs > int_start &
                                sub$os_yrs <= int_end + 1e-9),
      treatment  = factor(sub[[paste0("treatment_", k)]], levels = treatment_names),
      stringsAsFactors = FALSE
    )
    for (j in seq_along(cum_names)) row_k[[cum_names[j]]] <- cum_trt[, j]
    for (p in seq_len(P)) row_k[[covariate_spec$name[p]]] <- sub[[covariate_spec$name[p]]]
    if (use_tv) row_k$L <- sub[[paste0("L_", k)]]
    long_list[[k]] <- row_k
  }
  data_long <- do.call(rbind, long_list)
  data_long <- data_long[order(data_long$patient_id, data_long$line), ]
  rownames(data_long) <- NULL

  # ---- Step 7: return ----
  true_params <- list(
    beta_treatments = beta_treatments,
    beta_covariates = setNames(covariate_spec$beta, covariate_spec$name),
    covariate_names = covariate_spec$name,
    treatment_names = treatment_names,
    nonref_treatments = nonref_trts,
    cumulative_effect = cumulative_effect,
    baseline_log_hazard = baseline_log_hazard[seq_len(n_lines)],
    confounding_strength = confounding_strength,
    interval_length = interval_length,
    n_lines = n_lines,
    cens_rate = cens_rate
  )
  if (use_tv) true_params$tv_conf <- tv_conf

  list(data_wide = data_wide, data_long = data_long, true_params = true_params,
       version = "2026-11-10: next-line measurement for everyone entering a line; cens_fn option")
}
