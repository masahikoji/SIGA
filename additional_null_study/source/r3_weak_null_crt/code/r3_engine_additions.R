# Additional functions for the R3 simulation (weak-null heterogeneity and the
# variance-corrected reference test).  Base R only.  This file is sourced
# after siga_pair_path_engine.R and does not start a simulation.
#
# Notation (manuscript):
#   T_n          = 0.5 * sum_i z_i r_i                       (observed statistic)
#   V_S          = SIGA-S variance estimate                  (sampling)
#   V_R(d)       = V_S + n d'(Psi_hat - Pi_hat) d / 16       (conditional reference)
#   lambda(d)    = sqrt(V_R(d) / V_S)
#   corrected RT: compare lambda * T_n with the regenerated statistics T*_j.
#   With d = 0 (hence lambda = 1) the corrected test is the reference test.

options(stringsAsFactors = FALSE, warn = 1)

# -----------------------------------------------------------------------------
# Stratum-specific effect deviations
# -----------------------------------------------------------------------------

# Data-based estimator (Supplementary Appendix C.2):
#   d_hat_s(b) = mean(Y | A=1, S=s) - mean(Y | A=0, S=s) - b,   0 if an arm is empty.
# Also returns a plug-in variance estimate of d_hat_s, used by the debiased
# quadratic form below.
estimate_stratum_effect_deviations <- function(y, A, stratum_id, J, boundary) {
  y <- as.numeric(y)
  A <- as.integer(A)
  stratum_id <- as.integer(stratum_id)
  n1 <- tabulate(stratum_id[A == 1L], nbins = J)
  n0 <- tabulate(stratum_id[A == 0L], nbins = J)
  s1 <- numeric(J); s0 <- numeric(J); q1 <- numeric(J); q0 <- numeric(J)
  if (any(A == 1L)) {
    tmp <- rowsum(cbind(y[A == 1L], y[A == 1L]^2), group = stratum_id[A == 1L], reorder = FALSE)
    idx <- as.integer(rownames(tmp)); s1[idx] <- tmp[, 1L]; q1[idx] <- tmp[, 2L]
  }
  if (any(A == 0L)) {
    tmp <- rowsum(cbind(y[A == 0L], y[A == 0L]^2), group = stratum_id[A == 0L], reorder = FALSE)
    idx <- as.integer(rownames(tmp)); s0[idx] <- tmp[, 1L]; q0[idx] <- tmp[, 2L]
  }
  both <- n1 > 0L & n0 > 0L
  d <- numeric(J)
  d[both] <- s1[both] / n1[both] - s0[both] / n0[both] - boundary
  # unbiased within-arm variances (0 when an arm has fewer than two observations)
  v1 <- numeric(J); v0 <- numeric(J)
  ok1 <- n1 > 1L; ok0 <- n0 > 1L
  v1[ok1] <- pmax((q1[ok1] - s1[ok1]^2 / n1[ok1]) / (n1[ok1] - 1L), 0)
  v0[ok0] <- pmax((q0[ok0] - s0[ok0]^2 / n0[ok0]) / (n0[ok0] - 1L), 0)
  var_d <- numeric(J)
  var_d[both] <- v1[both] / n1[both] + v0[both] / n0[both]
  list(d = d, var_d = var_d, both_arms = both, n1 = n1, n0 = n0)
}

# Correction matrix M = Psi_hat - Pi_hat, quadratic form and lambda.
# debias = TRUE subtracts sum_s M_ss Var_hat(d_hat_s), the expected contribution of
# estimation noise in d_hat to the quadratic form.
corrected_variance <- function(sampling_variance_result, d, var_d = NULL, calibration,
                               debias = FALSE, floor_fraction = NULL) {
  d <- as.numeric(d)
  if (length(d) != calibration$J) stop("d must have one entry per joint profile.", call. = FALSE)
  counts <- sampling_variance_result$decomposition$counts
  pi_hat <- counts / sum(counts)
  M <- calibration$psi - diag(pi_hat, nrow = length(pi_hat))
  quad <- drop(crossprod(d, M %*% d))
  if (isTRUE(debias)) {
    if (is.null(var_d)) stop("var_d is required when debias = TRUE.", call. = FALSE)
    quad <- quad - sum(diag(M) * as.numeric(var_d))
  }
  correction <- calibration$n * quad / 16
  v_s <- sampling_variance_result$variance
  v_raw <- v_s + correction
  eps_n <- floor_fraction %||% (1 / calibration$n)    # manuscript: epsilon_n = 1/n
  v <- max(v_raw, eps_n * v_s, .Machine$double.eps)
  list(
    variance = v,
    variance_raw = v_raw,
    correction = correction,
    quadratic_form = quad,
    lambda = sqrt(v / v_s),
    floored = v_raw < eps_n * v_s
  )
}

# -----------------------------------------------------------------------------
# Reference randomisation test with rescaled critical values
# -----------------------------------------------------------------------------

# scores: n x L matrix of retained residual vectors (one column per boundary x analysis)
# lambda: L x C matrix of scale factors; lambda[l, 1] should be 1 for the uncorrected
#         reference test.  Column c compares lambda[l, c] * T_n with the regenerated
#         statistics of column l.
# Returns L x C matrices of inclusive plus-one Monte Carlo p-values (two-sided,
# upper, lower) and the empirical mean and variance of the regenerated statistics.
rt_corrected_pvalues_R <- function(X, scores, z_obs, lambda, B = 1999L, pbc = 0.80,
                                   weights_ = NULL, seed = 1, tolerance = 1e-12) {
  X <- check_binary_factor_matrix(X)
  scores <- as.matrix(scores)
  storage.mode(scores) <- "double"
  lambda <- as.matrix(lambda)
  z_obs <- as.integer(z_obs)
  n <- nrow(X); K <- ncol(X); L <- ncol(scores); C <- ncol(lambda)
  B <- as.integer(B)
  if (nrow(scores) != n) stop("scores must have n rows.", call. = FALSE)
  if (nrow(lambda) != L) stop("lambda must have one row per score column.", call. = FALSE)
  if (any(!is.finite(lambda)) || any(lambda <= 0)) stop("lambda must be finite and positive.", call. = FALSE)
  if (length(z_obs) != n) stop("z_obs must have length n.", call. = FALSE)
  if (B < 1L) stop("B must be positive.", call. = FALSE)
  weights <- as.numeric(weights_ %||% rep(1, K + 1L))
  validate_minimization_inputs(K, pbc, weights)
  set.seed(normalise_seed(seed))

  t_obs <- drop(0.5 * crossprod(z_obs, scores))
  t_rand <- matrix(0, nrow = B, ncol = L, dimnames = list(NULL, colnames(scores)))
  overall <- integer(B)
  marginal <- matrix(0L, nrow = B, ncol = 2L * K)

  for (i in seq_len(n)) {
    score_plus <- weights[1L] * abs(overall + 1L)
    score_minus <- weights[1L] * abs(overall - 1L)
    cols <- integer(K)
    for (j in seq_len(K)) {
      cols[j] <- 2L * (j - 1L) + X[i, j] + 1L
      current <- marginal[, cols[j]]
      score_plus <- score_plus + weights[j + 1L] * abs(current + 1L)
      score_minus <- score_minus + weights[j + 1L] * abs(current - 1L)
    }
    u <- runif(B)
    z <- integer(B)
    better_plus <- score_plus < score_minus
    better_minus <- score_plus > score_minus
    ties <- !(better_plus | better_minus)
    z[better_plus] <- ifelse(u[better_plus] < pbc, 1L, -1L)
    z[better_minus] <- ifelse(u[better_minus] < 1 - pbc, 1L, -1L)
    z[ties] <- ifelse(u[ties] < 0.5, 1L, -1L)
    overall <- overall + z
    for (j in seq_len(K)) marginal[, cols[j]] <- marginal[, cols[j]] + z
    t_rand <- t_rand + 0.5 * tcrossprod(z, scores[i, ])
  }

  p_two <- matrix(NA_real_, L, C); p_upper <- p_two; p_lower <- p_two
  for (l in seq_len(L)) {
    tr <- t_rand[, l]
    tie_signed <- abs(tr - t_obs[l]) <= tolerance              # T* = T_n
    tie_abs <- tie_signed | abs(tr + t_obs[l]) <= tolerance    # |T*| = |T_n|
    for (c in seq_len(C)) {
      thr <- lambda[l, c] * t_obs[l]
      # Regenerated statistics tied with the observed statistic always count as at least as
      # extreme (inclusive convention of the reference test); this clause is immaterial when
      # lambda = 1 and keeps the correction from excluding ties of lattice-valued statistics.
      p_two[l, c] <- (1 + sum((abs(tr) + tolerance >= abs(thr)) | tie_abs)) / (B + 1)
      p_upper[l, c] <- (1 + sum((tr + tolerance >= thr) | tie_signed)) / (B + 1)
      p_lower[l, c] <- (1 + sum((tr - tolerance <= thr) | tie_signed)) / (B + 1)
    }
  }
  dimnames(p_two) <- dimnames(p_upper) <- dimnames(p_lower) <-
    list(colnames(scores), colnames(lambda))
  list(
    statistic = t_obs,
    two_sided = p_two,
    greater = p_upper,
    less = p_lower,
    randomization_mean = colMeans(t_rand),
    randomization_variance = if (B > 1L) apply(t_rand, 2L, var) else rep(NA_real_, L),
    B = B
  )
}

# -----------------------------------------------------------------------------
# Effect-deviation directions
# -----------------------------------------------------------------------------

# Directions are defined on the J joint strata and then centred to pi'd = 0.
make_direction <- function(type, patterns, profile_prob, calibration = NULL,
                           beta_master = c(0.50, 0.40, 0.30, 0.20, 0.10)) {
  patterns <- check_binary_factor_matrix(patterns)
  K <- ncol(patterns)
  p_marg <- profile_marginals(patterns, profile_prob)
  v <- switch(
    type,
    zero = rep(0, nrow(patterns)),
    first_factor = patterns[, 1L] - p_marg[1L],
    main_effects = drop(sweep(patterns, 2L, p_marg, "-") %*% beta_master[seq_len(K)]),
    interaction = apply(2L * patterns - 1L, 1L, prod),
    min_ratio = {
      if (is.null(calibration)) stop("min_ratio needs a calibration object.", call. = FALSE)
      calibration$directions$min_direction
    },
    max_ratio = {
      if (is.null(calibration)) stop("max_ratio needs a calibration object.", call. = FALSE)
      calibration$directions$max_direction
    },
    stop("Unknown direction type: ", type, call. = FALSE)
  )
  weighted_center(v, profile_prob)
}

scale_centred_direction <- function(v, profile_prob, max_abs_d) {
  v <- weighted_center(v, profile_prob)
  m <- max(abs(v))
  if (m <= 0 || max_abs_d <= 0) return(rep(0, length(v)))
  v * (max_abs_d / m)
}

# -----------------------------------------------------------------------------
# Outcome models with a separate true marginal effect and stratum deviations
# -----------------------------------------------------------------------------
#  Every model carries:
#    $type             engine type used by generate_outcome_from_model()
#    $true_effect      marginal effect Delta = E{Y(1)-Y(0)}
#    $stratum_effect   Delta_s = E{Y(1)-Y(0) | S=s}
#    $d0               Delta_s - Delta (pi-weighted mean zero)
#  so that the model-based deviation at a null value b is  d(b) = stratum_effect - b.

make_continuous_model_r3 <- function(patterns, profile_prob, true_effect, d0,
                                     outcome_sd = 1.0, individual_effect_sd = 0.25) {
  d0 <- weighted_center(d0, profile_prob)
  mu0 <- pattern_baseline_mean(patterns, profile_prob)
  list(
    type = "continuous_realistic",
    boundary = true_effect,                      # field name used by the engine generator
    true_effect = true_effect,
    d0 = d0,
    stratum_effect = true_effect + d0,
    delta = true_effect + d0,                    # engine field
    mu0 = mu0,
    outcome_sd = outcome_sd,
    individual_effect_sd = if (max(abs(d0)) > 0) individual_effect_sd else 0,
    achieved_effect = sum(profile_prob * (true_effect + d0))
  )
}

make_binary_direct_model_r3 <- function(patterns, profile_prob, true_effect, d0,
                                        target_control_risk = 0.60,
                                        prob_floor = 0.02) {
  baseline <- calibrate_logistic_baseline(
    patterns = patterns, profile_prob = profile_prob,
    target_control_risk = target_control_risk
  )
  d0 <- weighted_center(d0, profile_prob)
  p0 <- baseline$p0
  shrink <- 1
  repeat {
    p1 <- p0 + true_effect + shrink * d0
    if (all(p1 > prob_floor & p1 < 1 - prob_floor)) break
    shrink <- shrink * 0.9
    if (shrink < 1e-3) stop("Binary heterogeneity cannot be made valid: true effect too large.", call. = FALSE)
  }
  d0 <- shrink * d0
  if (any(p0 <= prob_floor | p0 >= 1 - prob_floor)) stop("Invalid control risks.", call. = FALSE)
  list(
    type = "binary_direct_probability",
    boundary = true_effect,
    true_effect = true_effect,
    d0 = d0,
    stratum_effect = p1 - p0,
    p0 = p0,
    p1 = p1,
    heterogeneity_shrink = shrink,
    achieved_control_risk = sum(profile_prob * p0),
    achieved_treatment_risk = sum(profile_prob * p1),
    achieved_effect = sum(profile_prob * (p1 - p0))
  )
}

make_binary_common_log_odds_model_r3 <- function(patterns, profile_prob, true_effect,
                                                 target_control_risk = 0.60) {
  m <- calibrate_common_log_odds_model(
    patterns = patterns, profile_prob = profile_prob,
    target_control_risk = target_control_risk,
    target_risk_difference = true_effect
  )
  m$true_effect <- true_effect
  m$stratum_effect <- m$p1 - m$p0
  m$d0 <- m$stratum_effect - true_effect
  m$boundary <- true_effect
  m
}

model_based_deviation <- function(model, boundary) {
  as.numeric(model$stratum_effect) - boundary
}

# -----------------------------------------------------------------------------
# Optional: lattice mixture of continuity-corrected normal probabilities
# (Supplementary Appendix D) for the unadjusted binary superiority test at b = 0
# -----------------------------------------------------------------------------
lattice_mixture_pvalue <- function(statistic, variance, n, y_plus, treated_prob,
                                   alternative = c("two.sided", "greater", "less"),
                                   tolerance = 1e-12) {
  alternative <- match.arg(alternative)
  sdv <- sqrt(max(variance, .Machine$double.eps))
  total <- 0
  wsum <- 0
  for (k in which(treated_prob > 0) - 1L) {
    lo <- max(0L, k - (n - y_plus)); hi <- min(k, y_plus)
    if (lo > hi) next
    x <- lo:hi
    t <- x - k * y_plus / n
    g <- pnorm((t + 0.5) / sdv) - pnorm((t - 0.5) / sdv)
    G <- sum(g)
    if (!(G > 0)) next
    sel <- switch(
      alternative,
      two.sided = abs(t) + tolerance >= abs(statistic),
      greater = t + tolerance >= statistic,
      less = t - tolerance <= statistic
    )
    total <- total + treated_prob[k + 1L] * sum(g[sel]) / G
    wsum <- wsum + treated_prob[k + 1L]
  }
  if (wsum <= 0) return(NA_real_)
  min(max(total / wsum, 0), 1)
}

# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------
binomial_se <- function(p, n) sqrt(pmax(p * (1 - p), 0) / n)

paired_difference_se <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y)
  d <- x[keep] - y[keep]
  if (length(d) < 2L) return(NA_real_)
  sd(d) / sqrt(length(d))
}

# -----------------------------------------------------------------------------
# Optional compiled kernel for the regeneration loop (rt_kernel.c, base R .Call)
# -----------------------------------------------------------------------------
# load_rt_kernel() compiles rt_kernel.c with `R CMD SHLIB` on first use (a C compiler
# is required: Xcode command-line tools on macOS, build-essential on Linux) and returns
# TRUE when the compiled kernel is available.  rt_corrected_pvalues_R() is used otherwise.
.r3_kernel_state <- new.env(parent = emptyenv())
load_rt_kernel <- function(source_dir = NULL, quiet = FALSE) {
  if (isTRUE(.r3_kernel_state$loaded)) return(TRUE)
  if (isTRUE(.r3_kernel_state$failed)) return(FALSE)
  source_dir <- source_dir %||% (if (exists("BASE_DIR", inherits = TRUE)) get("BASE_DIR", inherits = TRUE) else getwd())
  src <- file.path(source_dir, "rt_kernel.c")
  if (!file.exists(src)) { .r3_kernel_state$failed <- TRUE; return(FALSE) }
  build_dir <- Sys.getenv("R3_BUILD_DIR", unset = "")
  if (!nzchar(build_dir)) build_dir <- file.path(source_dir, "build", paste0(R.version$platform, "_R", getRversion()))
  dir.create(build_dir, showWarnings = FALSE, recursive = TRUE)
  lib <- file.path(build_dir, paste0("rt_kernel", .Platform$dynlib.ext))
  if (!file.exists(lib) || file.mtime(lib) < file.mtime(src)) {
    ok <- tryCatch({
      file.copy(src, file.path(build_dir, "rt_kernel.c"), overwrite = TRUE)
      old <- setwd(build_dir); on.exit(setwd(old), add = TRUE)
      status <- system2(file.path(R.home("bin"), "R"), c("CMD", "SHLIB", "rt_kernel.c"),
                        stdout = if (quiet) FALSE else "", stderr = if (quiet) FALSE else "")
      status == 0L && file.exists(lib)
    }, error = function(e) FALSE)
    if (!isTRUE(ok)) { .r3_kernel_state$failed <- TRUE; if (!quiet) message("rt_kernel.c could not be compiled; using the pure-R kernel."); return(FALSE) }
  }
  ok <- tryCatch({ dyn.load(lib); TRUE }, error = function(e) FALSE)
  if (!ok) { .r3_kernel_state$failed <- TRUE; return(FALSE) }
  .r3_kernel_state$loaded <- TRUE
  .r3_kernel_state$lib <- lib
  TRUE
}

rt_corrected_pvalues_C <- function(X, scores, z_obs, lambda, B = 1999L, pbc = 0.80,
                                   weights_ = NULL, seed = 1, tolerance = 1e-12) {
  X <- check_binary_factor_matrix(X)
  scores <- as.matrix(scores); storage.mode(scores) <- "double"
  lambda <- as.matrix(lambda); z_obs <- as.integer(z_obs)
  n <- nrow(X); K <- ncol(X); L <- ncol(scores); C <- ncol(lambda); B <- as.integer(B)
  if (nrow(scores) != n) stop("scores must have n rows.", call. = FALSE)
  if (nrow(lambda) != L) stop("lambda must have one row per score column.", call. = FALSE)
  if (any(!is.finite(lambda)) || any(lambda <= 0)) stop("lambda must be finite and positive.", call. = FALSE)
  weights <- as.numeric(weights_ %||% rep(1, K + 1L)); validate_minimization_inputs(K, pbc, weights)
  set.seed(normalise_seed(seed))
  t_obs <- drop(0.5 * crossprod(z_obs, scores))
  t_rand <- .Call("C_rt_regenerate", X, scores, B, as.numeric(pbc), weights, PACKAGE = "rt_kernel")
  colnames(t_rand) <- colnames(scores)
  p_two <- matrix(NA_real_, L, C); p_upper <- p_two; p_lower <- p_two
  for (l in seq_len(L)) {
    tr <- t_rand[, l]
    tie_signed <- abs(tr - t_obs[l]) <= tolerance              # T* = T_n
    tie_abs <- tie_signed | abs(tr + t_obs[l]) <= tolerance    # |T*| = |T_n|
    for (c in seq_len(C)) {
      thr <- lambda[l, c] * t_obs[l]
      # Regenerated statistics tied with the observed statistic always count as at least as
      # extreme (inclusive convention of the reference test); this clause is immaterial when
      # lambda = 1 and keeps the correction from excluding ties of lattice-valued statistics.
      p_two[l, c] <- (1 + sum((abs(tr) + tolerance >= abs(thr)) | tie_abs)) / (B + 1)
      p_upper[l, c] <- (1 + sum((tr + tolerance >= thr) | tie_signed)) / (B + 1)
      p_lower[l, c] <- (1 + sum((tr - tolerance <= thr) | tie_signed)) / (B + 1)
    }
  }
  dimnames(p_two) <- dimnames(p_upper) <- dimnames(p_lower) <- list(colnames(scores), colnames(lambda))
  list(statistic = t_obs, two_sided = p_two, greater = p_upper, less = p_lower,
       randomization_mean = colMeans(t_rand),
       randomization_variance = if (B > 1L) apply(t_rand, 2L, var) else rep(NA_real_, L), B = B)
}

# Dispatcher: uses the compiled kernel when R3_KERNEL is "auto" (default) or "C" and it loads.
rt_corrected_pvalues <- function(..., kernel = Sys.getenv("R3_KERNEL", unset = "auto")) {
  kernel <- tolower(kernel)
  if (kernel %in% c("auto", "c") && load_rt_kernel(quiet = TRUE)) return(rt_corrected_pvalues_C(...))
  if (kernel == "c") stop("R3_KERNEL=C was requested but rt_kernel.c could not be compiled or loaded.", call. = FALSE)
  rt_corrected_pvalues_R(...)
}
