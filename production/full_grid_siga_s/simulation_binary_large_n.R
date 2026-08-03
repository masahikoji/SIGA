#!/usr/bin/env Rscript

# =============================================================================
# PURE-R binary superiority, non-inferiority, and equivalence under
# biased-coin Pocock--Simon minimization
#
# Proposed method:
#   Stratum-imbalance Gaussian approximation (SIGA).
# Reference method:
#   Fixed-score randomization test using regenerated allocation paths from the
#   same minimization algorithm.
#
# Nonzero risk-difference boundaries use W_i(b) = Y_i - b A_i.  These are
# weak-null procedures and are not claimed to be finite-sample exact Fisher
# tests.  The program evaluates their type I error directly at the lower and
# upper null boundaries.
#
# Final benchmark configuration (matched to the continuous-outcome study):
#   PWRT_N_OUTER=100000
#   PWRT_N_RERAND=4999
#   PWRT_N_CALIBRATION=100000
#
# The reference randomization test is implemented by the same fixed-score
# algorithm as in the continuous-outcome benchmark: condition on the ordered
# factor path, regenerate allocations from the same minimization rule, hold the
# boundary score fixed, and use an inclusive plus-one Monte Carlo p-value.
#
# Pure-R tuning:
#   PWRT_PURE_R_CALIBRATION_BATCH=2000 (increase if memory permits)
#
# Small smoke test:
#   PWRT_N_OUTER=20 PWRT_N_RERAND=99 PWRT_N_CALIBRATION=1000 \
#     Rscript simulation_binary_SIGA_weak_null_pureR_100k_4999_all_in_one.R
#
# Sharded run, for example shard 3 of 36:
#   PWRT_N_SHARDS=36 PWRT_SHARD_ID=3 Rscript simulation_binary_SIGA_weak_null_pureR_100k_4999_all_in_one.R
#
# Aggregate after all shards finish:
#   PWRT_MODE=aggregate PWRT_N_SHARDS=36 Rscript simulation_binary_SIGA_weak_null_pureR_100k_4999_all_in_one.R
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

script_directory <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(
      sub("^--file=", "", file_arg[1L]),
      winslash = "/", mustWork = FALSE
    )))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

BASE_DIR <- script_directory()
# =============================================================================
# Pure-R engine and analysis functions (embedded; no external source file)
# =============================================================================

# Pure-R helper functions for binary superiority, non-inferiority, and
# equivalence under biased-coin Pocock--Simon minimization.
#
# No C++ or Rcpp is used.  The randomization-test engine is vectorized across
# regenerated allocation paths, while retaining the sequential minimization
# rule across participants.
#
# Nonzero risk-difference boundaries use W_i(b) = Y_i - b A_i.  The resulting
# reference procedure is a weak-null randomization procedure rather than a
# finite-sample exact Fisher test.

`%||%` <- function(x, y) if (is.null(x)) y else x

check_binary_factor_matrix <- function(X) {
  X <- as.matrix(X)
  storage.mode(X) <- "integer"
  if (anyNA(X) || any(!X %in% c(0L, 1L))) {
    stop("All minimization factors must be coded 0/1.", call. = FALSE)
  }
  X
}

joint_stratum_id <- function(X) {
  X <- check_binary_factor_matrix(X)
  K <- ncol(X)
  as.integer(1L + drop(X %*% (2L ^ (0:(K - 1L)))))
}

nearest_psd <- function(M, tolerance = 1e-12) {
  M <- (M + t(M)) / 2
  ee <- eigen(M, symmetric = TRUE)
  cutoff <- tolerance * max(1, max(abs(ee$values)))
  values <- pmax(ee$values, 0)
  values[values < cutoff] <- 0
  out <- ee$vectors %*% (values * t(ee$vectors))
  (out + t(out)) / 2
}

normalise_seed <- function(seed) {
  seed <- as.double(seed)
  if (!is.finite(seed)) stop("seed must be finite.", call. = FALSE)
  as.integer(abs(seed) %% 2147483646 + 1)
}

validate_minimization_inputs <- function(K, pbc, weights) {
  if (!(pbc >= 0.5 && pbc <= 1)) {
    stop("pbc must lie in [0.5, 1].", call. = FALSE)
  }
  if (length(weights) != K + 1L) {
    stop("weights must have length K + 1: overall plus one per factor.", call. = FALSE)
  }
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("weights must be finite and nonnegative.", call. = FALSE)
  }
  invisible(TRUE)
}

# Generate one allocation path for a fixed ordered factor matrix.
ps_assign_R <- function(X,
                        pbc = 0.80,
                        weights_ = NULL,
                        seed = 1) {
  X <- check_binary_factor_matrix(X)
  n <- nrow(X)
  K <- ncol(X)
  weights <- as.numeric(weights_ %||% rep(1, K + 1L))
  validate_minimization_inputs(K, pbc, weights)
  set.seed(normalise_seed(seed))

  z <- integer(n)
  overall <- 0L
  marginal <- matrix(0L, nrow = K, ncol = 2L)

  for (i in seq_len(n)) {
    score_plus <- weights[1L] * abs(overall + 1L)
    score_minus <- weights[1L] * abs(overall - 1L)
    for (j in seq_len(K)) {
      level <- X[i, j] + 1L
      current <- marginal[j, level]
      score_plus <- score_plus + weights[j + 1L] * abs(current + 1L)
      score_minus <- score_minus + weights[j + 1L] * abs(current - 1L)
    }

    u <- runif(1L)
    zi <- if (score_plus < score_minus) {
      if (u < pbc) 1L else -1L
    } else if (score_plus > score_minus) {
      if (u < 1 - pbc) 1L else -1L
    } else {
      if (u < 0.5) 1L else -1L
    }

    z[i] <- zi
    overall <- overall + zi
    for (j in seq_len(K)) {
      marginal[j, X[i, j] + 1L] <- marginal[j, X[i, j] + 1L] + zi
    }
  }
  z
}

# Allocation-only calibration for independent binary factors.  Paths are
# processed in vectorized batches to avoid an R loop over calibration paths.
precompute_design_R <- function(B0,
                                n,
                                probs,
                                pbc = 0.80,
                                weights_ = NULL,
                                seed = 1,
                                batch_size = NULL,
                                progress = TRUE) {
  B0 <- as.integer(B0)
  n <- as.integer(n)
  probs <- as.numeric(probs)
  K <- length(probs)
  if (B0 < 2L) stop("B0 must be at least 2.", call. = FALSE)
  if (n < 2L) stop("n must be at least 2.", call. = FALSE)
  if (K < 1L || K > 20L) stop("K must be between 1 and 20.", call. = FALSE)
  if (any(!is.finite(probs)) || any(probs <= 0 | probs >= 1)) {
    stop("Each binary factor probability must lie strictly between 0 and 1.", call. = FALSE)
  }
  weights <- as.numeric(weights_ %||% rep(1, K + 1L))
  validate_minimization_inputs(K, pbc, weights)
  J <- 2L ^ K

  if (is.null(batch_size)) {
    env_batch <- suppressWarnings(as.integer(Sys.getenv(
      "PWRT_PURE_R_CALIBRATION_BATCH", unset = "2000"
    )))
    batch_size <- if (is.na(env_batch) || env_batch < 1L) 2000L else env_batch
  }
  batch_size <- max(1L, min(as.integer(batch_size), B0))
  set.seed(normalise_seed(seed))

  sum_u <- numeric(J)
  sum_uu <- matrix(0, nrow = J, ncol = J)
  treated_count <- numeric(n + 1L)
  completed <- 0L
  start <- proc.time()[3L]

  while (completed < B0) {
    m <- min(batch_size, B0 - completed)
    rows <- seq_len(m)
    overall <- integer(m)
    marginal <- matrix(0L, nrow = m, ncol = 2L * K)
    counts <- matrix(0L, nrow = m, ncol = J)
    imbalance <- matrix(0L, nrow = m, ncol = J)

    for (i in seq_len(n)) {
      score_plus <- weights[1L] * abs(overall + 1L)
      score_minus <- weights[1L] * abs(overall - 1L)
      joint0 <- integer(m)
      x_list <- vector("list", K)

      for (j in seq_len(K)) {
        xj <- rbinom(m, size = 1L, prob = probs[j])
        x_list[[j]] <- xj
        joint0 <- joint0 + xj * (2L ^ (j - 1L))
        col_j <- 2L * (j - 1L) + xj + 1L
        idx_j <- rows + (col_j - 1L) * m
        current <- marginal[idx_j]
        score_plus <- score_plus + weights[j + 1L] * abs(current + 1L)
        score_minus <- score_minus + weights[j + 1L] * abs(current - 1L)
      }

      u <- runif(m)
      zi <- integer(m)
      better_plus <- score_plus < score_minus
      better_minus <- score_plus > score_minus
      ties <- !(better_plus | better_minus)
      zi[better_plus] <- ifelse(u[better_plus] < pbc, 1L, -1L)
      zi[better_minus] <- ifelse(u[better_minus] < 1 - pbc, 1L, -1L)
      zi[ties] <- ifelse(u[ties] < 0.5, 1L, -1L)

      overall <- overall + zi
      for (j in seq_len(K)) {
        col_j <- 2L * (j - 1L) + x_list[[j]] + 1L
        idx_j <- rows + (col_j - 1L) * m
        marginal[idx_j] <- marginal[idx_j] + zi
      }
      idx_s <- rows + joint0 * m
      counts[idx_s] <- counts[idx_s] + 1L
      imbalance[idx_s] <- imbalance[idx_s] + zi
    }

    U <- matrix(0, nrow = m, ncol = J)
    nonzero <- counts > 0L
    U[nonzero] <- imbalance[nonzero] / sqrt(counts[nonzero])
    sum_u <- sum_u + colSums(U)
    sum_uu <- sum_uu + crossprod(U)
    treated <- as.integer((n + overall) / 2L)
    treated_count <- treated_count + tabulate(treated + 1L, nbins = n + 1L)

    completed <- completed + m
    if (isTRUE(progress)) {
      message(
        "  pure-R calibration: ", completed, "/", B0,
        " paths; elapsed ", sprintf("%.1f min", (proc.time()[3L] - start) / 60)
      )
    }
  }

  mean_u <- sum_u / B0
  gamma <- (sum_uu - B0 * tcrossprod(mean_u)) / (B0 - 1L)
  list(
    gamma = gamma,
    mean_u = mean_u,
    treated_count_prob = treated_count / B0,
    B0 = B0,
    n = n,
    K = K,
    J = J,
    probs = probs,
    pbc = pbc,
    weights = weights,
    seed = seed,
    elapsed_seconds = proc.time()[3L] - start,
    engine = "pure_R_vectorized"
  )
}

# Reference p-values for one or more fixed score vectors.  The B regenerated
# allocation paths are evolved simultaneously; the treatment rule remains
# sequential in participant order.
rt_pvalues_R <- function(X,
                         scores,
                         z_obs,
                         B = 1999L,
                         pbc = 0.80,
                         weights_ = NULL,
                         seed = 1,
                         tolerance = 1e-12) {
  X <- check_binary_factor_matrix(X)
  scores <- as.matrix(scores)
  storage.mode(scores) <- "double"
  z_obs <- as.integer(z_obs)
  n <- nrow(X)
  K <- ncol(X)
  L <- ncol(scores)
  B <- as.integer(B)
  if (nrow(scores) != n) stop("scores must have n rows.", call. = FALSE)
  if (length(z_obs) != n) stop("z_obs must have length n.", call. = FALSE)
  if (B < 1L) stop("B must be positive.", call. = FALSE)
  weights <- as.numeric(weights_ %||% rep(1, K + 1L))
  validate_minimization_inputs(K, pbc, weights)
  set.seed(normalise_seed(seed))

  t_obs <- drop(0.5 * crossprod(z_obs, scores))
  t_rand <- matrix(0, nrow = B, ncol = L)
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
    for (ell in seq_len(L)) {
      t_rand[, ell] <- t_rand[, ell] + 0.5 * z * scores[i, ell]
    }
  }

  abs_obs <- matrix(abs(t_obs), nrow = B, ncol = L, byrow = TRUE)
  obs_mat <- matrix(t_obs, nrow = B, ncol = L, byrow = TRUE)
  p_two <- (1 + colSums(abs(t_rand) + tolerance >= abs_obs)) / (B + 1)
  p_upper <- (1 + colSums(t_rand + tolerance >= obs_mat)) / (B + 1)
  p_lower <- (1 + colSums(t_rand - tolerance <= obs_mat)) / (B + 1)

  list(
    statistic = t_obs,
    two_sided = p_two,
    greater = p_upper,
    less = p_lower,
    B = B,
    engine = "pure_R_vectorized"
  )
}

calibrate_siga_binary <- function(n,
                                  factor_prob,
                                  pbc = 0.80,
                                  weights = rep(1, length(factor_prob) + 1L),
                                  B0 = 100000L,
                                  seed = 20260724,
                                  cpp_file = NULL,
                                  compile = FALSE,
                                  batch_size = NULL) {
  raw <- precompute_design_R(
    B0 = as.integer(B0),
    n = as.integer(n),
    probs = as.numeric(factor_prob),
    pbc = as.numeric(pbc),
    weights_ = as.numeric(weights),
    seed = as.double(seed),
    batch_size = batch_size,
    progress = TRUE
  )
  raw$gamma <- nearest_psd(as.matrix(raw$gamma))
  class(raw) <- c("siga_binary_calibration", "list")
  raw
}

all_binary_patterns <- function(K) {
  out <- as.matrix(expand.grid(rep(list(c(0L, 1L)), K)))
  storage.mode(out) <- "integer"
  colnames(out) <- paste0("X", seq_len(K))
  out
}

pattern_probabilities <- function(patterns, factor_prob) {
  factor_prob <- as.numeric(factor_prob)
  probs <- apply(patterns, 1L, function(x) {
    prod(ifelse(x == 1L, factor_prob, 1 - factor_prob))
  })
  probs / sum(probs)
}

calibrate_binary_outcome_model <- function(factor_prob,
                                           target_control_risk,
                                           target_risk_difference,
                                           beta,
                                           interaction = 0.20,
                                           root_interval = c(-30, 30),
                                           tolerance = 1e-12) {
  K <- length(factor_prob)
  beta <- as.numeric(beta)[seq_len(K)]
  patterns <- all_binary_patterns(K)
  pi_x <- pattern_probabilities(patterns, factor_prob)
  nonintercept <- drop(patterns %*% beta)
  if (K >= 2L) {
    nonintercept <- nonintercept + interaction * patterns[, 1L] * patterns[, 2L]
  }

  mean_control <- function(alpha0) {
    sum(pi_x * plogis(alpha0 + nonintercept))
  }
  intercept <- uniroot(
    function(a) mean_control(a) - target_control_risk,
    interval = root_interval,
    tol = tolerance
  )$root

  p0_x <- plogis(intercept + nonintercept)
  marginal_rd <- function(theta) {
    sum(pi_x * (plogis(intercept + nonintercept + theta) - p0_x))
  }
  theta <- uniroot(
    function(th) marginal_rd(th) - target_risk_difference,
    interval = root_interval,
    tol = tolerance
  )$root

  p1_x <- plogis(intercept + nonintercept + theta)
  structure(
    list(
      intercept = intercept,
      theta = theta,
      beta = beta,
      interaction = interaction,
      target_control_risk = target_control_risk,
      target_risk_difference = target_risk_difference,
      achieved_control_risk = sum(pi_x * p0_x),
      achieved_treatment_risk = sum(pi_x * p1_x),
      achieved_risk_difference = sum(pi_x * (p1_x - p0_x)),
      patterns = patterns,
      pattern_prob = pi_x
    ),
    class = c("binary_outcome_model", "list")
  )
}

generate_factors <- function(n, factor_prob) {
  K <- length(factor_prob)
  X <- matrix(0L, nrow = n, ncol = K)
  for (j in seq_len(K)) {
    X[, j] <- rbinom(n, size = 1L, prob = factor_prob[j])
  }
  colnames(X) <- paste0("X", seq_len(K))
  X
}

binary_event_probabilities <- function(X, model) {
  X <- check_binary_factor_matrix(X)
  eta0 <- model$intercept + drop(X %*% model$beta)
  if (ncol(X) >= 2L) {
    eta0 <- eta0 + model$interaction * X[, 1L] * X[, 2L]
  }
  list(p0 = plogis(eta0), p1 = plogis(eta0 + model$theta))
}

generate_binary_outcome <- function(X, A, model) {
  pp <- binary_event_probabilities(X, model)
  prob <- ifelse(A == 1L, pp$p1, pp$p0)
  rbinom(length(A), size = 1L, prob = prob)
}

make_boundary_scores <- function(y, A, X, boundaries) {
  y <- as.numeric(y)
  A <- as.numeric(A)
  X <- as.matrix(X)
  boundaries <- as.numeric(boundaries)
  n <- length(y)
  if (length(A) != n || nrow(X) != n) {
    stop("y, A, and X must have the same number of participants.", call. = FALSE)
  }

  C <- cbind(`(Intercept)` = 1, X)
  qr_C <- qr(C)
  unadjusted <- matrix(NA_real_, nrow = n, ncol = length(boundaries))
  adjusted <- matrix(NA_real_, nrow = n, ncol = length(boundaries))

  for (j in seq_along(boundaries)) {
    W <- y - boundaries[j] * A
    unadjusted[, j] <- W - mean(W)
    adjusted[, j] <- qr.resid(qr_C, W)
    adjusted[, j] <- adjusted[, j] - mean(adjusted[, j])
  }
  colnames(unadjusted) <- paste0("b", seq_along(boundaries))
  colnames(adjusted) <- paste0("b", seq_along(boundaries))
  list(unadjusted = unadjusted, adjusted = adjusted)
}

score_decomposition <- function(score, X, J) {
  score <- as.numeric(score)
  id <- joint_stratum_id(X)
  counts <- tabulate(id, nbins = J)
  sums <- numeric(J)
  tmp <- rowsum(score, group = id, reorder = FALSE)
  sums[as.integer(rownames(tmp))] <- drop(tmp)
  means <- numeric(J)
  nonempty <- counts > 0L
  means[nonempty] <- sums[nonempty] / counts[nonempty]
  residual <- score - means[id]
  list(id = id, counts = counts, means = means, residual = residual)
}

siga_variance <- function(score, X, calibration, truncate_kappa = TRUE) {
  if (!inherits(calibration, "siga_binary_calibration")) {
    stop("calibration must be produced by calibrate_siga_binary().", call. = FALSE)
  }
  X <- check_binary_factor_matrix(X)
  if (nrow(X) != calibration$n || ncol(X) != calibration$K) {
    stop("Trial dimensions do not match the allocation calibration.", call. = FALSE)
  }

  dec <- score_decomposition(score, X, calibration$J)
  sqrt_counts <- sqrt(dec$counts)
  Omega <- tcrossprod(sqrt_counts) * calibration$gamma
  nonempty <- dec$counts > 0L
  Jplus <- sum(nonempty)
  trace_low <- sum(diag(Omega)[nonempty] / dec$counts[nonempty])
  denominator <- calibration$n - Jplus
  kappa <- if (denominator > 0) {
    (calibration$n - trace_low) / denominator
  } else {
    0
  }
  if (truncate_kappa && kappa < 0) kappa <- 0

  low <- drop(crossprod(dec$means, Omega %*% dec$means))
  within <- sum(dec$residual^2)
  variance <- 0.25 * (low + kappa * within)
  if (!is.finite(variance) || variance <= 0) variance <- .Machine$double.eps

  list(
    variance = variance,
    kappa = kappa,
    low_component = 0.25 * low,
    within_component = 0.25 * kappa * within,
    decomposition = dec
  )
}

gaussian_pvalue <- function(statistic, variance,
                            alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  z <- statistic / sqrt(variance)
  switch(
    alternative,
    two.sided = 2 * pnorm(-abs(z)),
    greater = pnorm(z, lower.tail = FALSE),
    less = pnorm(z)
  )
}

lattice_normal_pvalue <- function(statistic,
                                  total_success,
                                  variance,
                                  treated_count_prob,
                                  alternative = c("two.sided", "greater", "less"),
                                  tolerance = 1e-12) {
  alternative <- match.arg(alternative)
  n <- length(treated_count_prob) - 1L
  S <- as.integer(total_success)
  if (S < 0L || S > n) stop("total_success must be between 0 and n.")
  if (S == 0L || S == n) return(1)

  wk <- as.numeric(treated_count_prob)
  wk <- wk / sum(wk)
  sdv <- sqrt(max(variance, .Machine$double.eps))
  ans <- 0

  for (k in 0:n) {
    if (wk[k + 1L] <= 0) next
    lo <- max(0L, k - (n - S))
    hi <- min(k, S)
    rr <- lo:hi
    support <- rr - k * S / n
    mass <- pnorm((support + 0.5) / sdv) -
      pnorm((support - 0.5) / sdv)
    normalizer <- sum(mass)
    if (!is.finite(normalizer) || normalizer <= 0) next
    keep <- switch(
      alternative,
      two.sided = abs(support) + tolerance >= abs(statistic),
      greater = support + tolerance >= statistic,
      less = support - tolerance <= statistic
    )
    ans <- ans + wk[k + 1L] * sum(mass[keep]) / normalizer
  }
  min(max(ans, 0), 1)
}

siga_component_pvalue <- function(score,
                                  z,
                                  X,
                                  calibration,
                                  alternative,
                                  use_binary_lattice = FALSE,
                                  y = NULL) {
  statistic <- 0.5 * sum(as.numeric(z) * as.numeric(score))
  vv <- siga_variance(score, X, calibration)
  if (use_binary_lattice) {
    if (is.null(y)) stop("Supply y when use_binary_lattice=TRUE.")
    p <- lattice_normal_pvalue(
      statistic = statistic,
      total_success = sum(y),
      variance = vv$variance,
      treated_count_prob = calibration$treated_count_prob,
      alternative = alternative
    )
  } else {
    p <- gaussian_pvalue(statistic, vv$variance, alternative)
  }
  list(
    statistic = statistic,
    variance = vv$variance,
    z = statistic / sqrt(vv$variance),
    p = p,
    kappa = vv$kappa
  )
}

reference_component_pvalue <- function(rt, column, alternative) {
  switch(
    alternative,
    two.sided = rt$two_sided[column],
    greater = rt$greater[column],
    less = rt$less[column]
  )
}

objective_specification <- function(objective,
                                    ni_margin = 0.10,
                                    equivalence_limits = c(-0.10, 0.10)) {
  objective <- match.arg(
    objective,
    c("superiority", "noninferiority", "equivalence")
  )
  switch(
    objective,
    superiority = list(boundaries = 0, alternatives = "two.sided"),
    noninferiority = list(
      boundaries = -abs(ni_margin),
      alternatives = "greater"
    ),
    equivalence = list(
      boundaries = sort(as.numeric(equivalence_limits)),
      alternatives = c("greater", "less")
    )
  )
}

combine_objective_pvalues <- function(component_p, objective) {
  if (objective == "equivalence") max(component_p) else component_p[1L]
}

wilson_interval <- function(x, n, conf.level = 0.95) {
  if (!is.finite(n) || n <= 0) {
    return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  z <- qnorm(1 - (1 - conf.level) / 2)
  phat <- x / n
  denominator <- 1 + z^2 / n
  center <- (phat + z^2 / (2 * n)) / denominator
  half <- z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2)) /
    denominator
  c(
    estimate = phat,
    lower = max(0, center - half),
    upper = min(1, center + half)
  )
}


# =============================================================================
# Environment helpers and paths
# =============================================================================

env_integer <- function(name, default) {
  value <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) return(as.integer(default))
  out <- suppressWarnings(as.integer(value))
  if (is.na(out)) stop(name, " must be an integer.", call. = FALSE)
  out
}

env_numeric <- function(name, default) {
  value <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) return(as.numeric(default))
  out <- suppressWarnings(as.numeric(value))
  if (!is.finite(out)) stop(name, " must be numeric.", call. = FALSE)
  out
}

env_integer_vector <- function(name) {
  value <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) return(integer())
  out <- suppressWarnings(as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1L]])))
  if (anyNA(out)) stop(name, " must contain comma-separated integers.", call. = FALSE)
  out
}

ensure_directory <- function(path, attempts = 5L, wait_seconds = 2) {
  path <- path.expand(path)
  for (attempt in seq_len(attempts)) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (dir.exists(path)) {
      return(normalizePath(path, winslash = "/", mustWork = TRUE))
    }
    if (attempt < attempts) Sys.sleep(wait_seconds)
  }
  stop("Could not create or access directory: ", path, call. = FALSE)
}

PROJECT_DIR <- path.expand(Sys.getenv(
  "PWRT_PROJECT_DIR",
  unset = "~/Dropbox/Tex/03_Post_PhD/029_exact_minimization/02_program"
))
OUTPUT_DIR <- path.expand(Sys.getenv(
  "PWRT_OUTPUT_DIR",
  unset = file.path(PROJECT_DIR, "binary_siga_weak_null_pureR_100k_4999_output")
))
PROJECT_DIR <- ensure_directory(PROJECT_DIR)
OUTPUT_DIR <- ensure_directory(OUTPUT_DIR)
CHUNK_DIR <- ensure_directory(file.path(OUTPUT_DIR, "scenario_shards"))
CALIBRATION_DIR <- ensure_directory(file.path(
  PROJECT_DIR, "binary_siga_allocation_calibration_cache_pureR"
))

# =============================================================================
# Simulation settings
# =============================================================================

RUN_MODE <- tolower(trimws(Sys.getenv("PWRT_MODE", unset = "run")))
if (!RUN_MODE %in% c("run", "aggregate")) {
  stop("PWRT_MODE must be 'run' or 'aggregate'.", call. = FALSE)
}

N_OUTER <- env_integer("PWRT_N_OUTER", 100000L)
N_RERANDOMIZATIONS <- env_integer("PWRT_N_RERAND", 4999L)
N_CALIBRATION <- env_integer("PWRT_N_CALIBRATION", 100000L)
N_SHARDS <- env_integer("PWRT_N_SHARDS", 1L)
SHARD_ID <- env_integer("PWRT_SHARD_ID", 1L)
OUTER_BATCH_SIZE <- env_integer("PWRT_OUTER_BATCH", 10L)
P_BIASED_COIN <- env_numeric("PWRT_P_BIASED_COIN", 0.80)
BASE_SEED <- env_integer("PWRT_SEED", 20260726L)
SCENARIO_FILTER <- env_integer_vector("PWRT_SCENARIO_IDS")
BINARY_ANALYSIS_MODE <- tolower(trimws(Sys.getenv(
  "PWRT_BINARY_ANALYSIS", unset = "joint"
)))
if (!BINARY_ANALYSIS_MODE %in% c("joint", "unadjusted", "adjusted")) {
  stop(
    "PWRT_BINARY_ANALYSIS must be 'joint', 'unadjusted', or 'adjusted'.",
    call. = FALSE
  )
}
ALLOW_PARTIAL <- identical(Sys.getenv("PWRT_ALLOW_PARTIAL", unset = "0"), "1")
STRICT_FINAL <- identical(Sys.getenv("PWRT_STRICT_FINAL", unset = "0"), "1")
BENCHMARK_VERSION <- "binary_siga_100k_4999_v1"

physical_cores <- parallel::detectCores(logical = FALSE)
if (is.na(physical_cores)) physical_cores <- parallel::detectCores(logical = TRUE)
if (is.na(physical_cores)) physical_cores <- 1L
DEFAULT_CORES <- if (.Platform$OS.type == "windows") {
  1L
} else {
  max(1L, min(8L, physical_cores - 1L))
}
N_CORES <- env_integer("PWRT_CORES", DEFAULT_CORES)
if (.Platform$OS.type == "windows") N_CORES <- 1L
USE_PARALLEL <- .Platform$OS.type == "unix" && N_CORES > 1L

ALPHA_SUPERIORITY <- 0.05
ALPHA_NONINFERIORITY <- 0.025
ALPHA_EQUIVALENCE <- 0.05
NONINFERIORITY_MARGIN <- 0.10
EQUIVALENCE_LIMITS <- c(-0.10, 0.10)
TARGET_CONTROL_RISK <- 0.60
BINARY_BETA_MASTER <- c(0.35, -0.25, 0.20, -0.15, 0.10)
BINARY_INTERACTION <- 0.20
FACTOR_PREVALENCE_MASTER <- c(0.50, 0.40, 0.30, 0.20, 0.10)

validate_settings <- function() {
  stopifnot(
    N_OUTER >= 1L,
    N_RERANDOMIZATIONS >= 1L,
    N_CALIBRATION >= 2L,
    N_SHARDS >= 1L,
    SHARD_ID >= 1L,
    SHARD_ID <= N_SHARDS,
    OUTER_BATCH_SIZE >= 1L,
    P_BIASED_COIN > 0.5,
    P_BIASED_COIN < 1,
    N_CORES >= 1L
  )
  expected <- c(
    outer_trials = 100000L,
    rerandomizations = 4999L,
    calibration_paths = 100000L
  )
  observed <- c(
    outer_trials = N_OUTER,
    rerandomizations = N_RERANDOMIZATIONS,
    calibration_paths = N_CALIBRATION
  )
  mismatch <- observed != expected
  if (any(mismatch)) {
    detail <- paste0(
      names(observed)[mismatch], "=", observed[mismatch],
      " (expected ", expected[mismatch], ")",
      collapse = "; "
    )
    if (STRICT_FINAL) {
      stop("Final binary benchmark configuration mismatch: ", detail, call. = FALSE)
    }
    warning("Non-final binary benchmark configuration: ", detail)
  }
  invisible(TRUE)
}
validate_settings()

# =============================================================================
# Reproducible seeds and checkpoint utilities
# =============================================================================

seed_value <- function(scenario_id, replicate_id, stream_id) {
  modulus <- 2147483646
  value <- as.double(BASE_SEED) +
    1000003 * as.double(stream_id) +
    104729 * as.double(scenario_id) +
    1009 * as.double(replicate_id)
  as.integer(value %% modulus + 1)
}

calibration_seed <- function(total_n, factor_count) {
  modulus <- 2147483646
  value <- as.double(BASE_SEED) + 7000003 +
    1009 * as.double(total_n) + 104729 * as.double(factor_count)
  as.integer(value %% modulus + 1)
}

append_csv <- function(data, path, attempts = 5L, wait_seconds = 2) {
  if (!nrow(data)) return(invisible(NULL))
  parent <- ensure_directory(dirname(path))
  last_error <- NULL

  for (attempt in seq_len(attempts)) {
    tmp <- NULL
    ok <- tryCatch({
      exists <- file.exists(path)
      tmp <- tempfile(
        pattern = paste0(basename(path), ".batch_"),
        tmpdir = parent,
        fileext = ".tmp"
      )
      write.table(
        data,
        file = tmp,
        sep = ",",
        row.names = FALSE,
        col.names = !exists,
        append = FALSE,
        quote = TRUE,
        qmethod = "double"
      )
      if (exists) {
        if (!file.append(path, tmp)) stop("file.append returned FALSE")
        unlink(tmp)
      } else if (!file.rename(tmp, path)) {
        if (!file.copy(tmp, path, overwrite = FALSE)) {
          stop("Could not move the first checkpoint batch into place")
        }
        unlink(tmp)
      }
      TRUE
    }, error = function(e) {
      last_error <<- conditionMessage(e)
      if (!is.null(tmp) && file.exists(tmp)) unlink(tmp)
      FALSE
    })
    if (isTRUE(ok)) return(invisible(NULL))
    if (attempt < attempts) Sys.sleep(wait_seconds)
  }
  stop(
    "Failed to append checkpoint: ", path,
    if (!is.null(last_error)) paste0("; ", last_error) else "",
    call. = FALSE
  )
}

format_elapsed <- function(seconds) {
  if (!is.finite(seconds)) return("NA")
  if (seconds < 60) return(sprintf("%.1f seconds", seconds))
  if (seconds < 3600) return(sprintf("%.1f minutes", seconds / 60))
  sprintf("%.2f hours", seconds / 3600)
}

# =============================================================================
# Scenario grid
# =============================================================================

make_scenarios <- function() {
  designs <- data.frame(
    factor_count = c(2L, 5L),
    target_per_group = c(500L, 1000L),
    superiority_power_rd = c(0.10, 0.07),
    noninferiority_power_rd = c(0.00, -0.03),
    equivalence_power_rd = c(0.00, 0.04),
    stringsAsFactors = FALSE
  )

  out <- list()
  pos <- 0L
  for (d in seq_len(nrow(designs))) {
    roles <- data.frame(
      objective = c(
        "superiority", "superiority",
        "noninferiority", "noninferiority",
        "equivalence", "equivalence", "equivalence"
      ),
      scenario_role = c(
        "type1", "power",
        "type1", "power",
        "type1_lower", "type1_upper", "power"
      ),
      true_risk_difference = c(
        0, designs$superiority_power_rd[d],
        -NONINFERIORITY_MARGIN, designs$noninferiority_power_rd[d],
        EQUIVALENCE_LIMITS[1L], EQUIVALENCE_LIMITS[2L],
        designs$equivalence_power_rd[d]
      ),
      alpha = c(
        ALPHA_SUPERIORITY, ALPHA_SUPERIORITY,
        ALPHA_NONINFERIORITY, ALPHA_NONINFERIORITY,
        ALPHA_EQUIVALENCE, ALPHA_EQUIVALENCE, ALPHA_EQUIVALENCE
      ),
      stringsAsFactors = FALSE
    )
    for (r in seq_len(nrow(roles))) {
      pos <- pos + 1L
      out[[pos]] <- data.frame(
        factor_count = designs$factor_count[d],
        target_per_group = designs$target_per_group[d],
        total_n = 2L * designs$target_per_group[d],
        objective = roles$objective[r],
        scenario_role = roles$scenario_role[r],
        true_risk_difference = roles$true_risk_difference[r],
        alpha = roles$alpha[r],
        stringsAsFactors = FALSE
      )
    }
  }
  tab <- do.call(rbind, out)
  tab$scenario_id <- seq_len(nrow(tab))
  tab$scenario_label <- paste0(
    "K", tab$factor_count,
    "_npg", tab$target_per_group,
    "_", tab$objective,
    "_", tab$scenario_role
  )
  if (length(SCENARIO_FILTER)) {
    tab <- tab[tab$scenario_id %in% SCENARIO_FILTER, , drop = FALSE]
  }
  rownames(tab) <- NULL
  tab
}

SCENARIOS <- make_scenarios()
if (!nrow(SCENARIOS)) stop("No scenarios remain after filtering.", call. = FALSE)

factor_probabilities <- function(factor_count) {
  FACTOR_PREVALENCE_MASTER[seq_len(factor_count)]
}

scenario_specification <- function(scenario) {
  objective_specification(
    objective = scenario$objective,
    ni_margin = NONINFERIORITY_MARGIN,
    equivalence_limits = EQUIVALENCE_LIMITS
  )
}

# Calibrate each distinct binary outcome model once, before parallel workers run.
OUTCOME_MODELS <- setNames(
  lapply(seq_len(nrow(SCENARIOS)), function(i) {
    scenario <- SCENARIOS[i, ]
    calibrate_binary_outcome_model(
      factor_prob = factor_probabilities(scenario$factor_count),
      target_control_risk = TARGET_CONTROL_RISK,
      target_risk_difference = scenario$true_risk_difference,
      beta = BINARY_BETA_MASTER,
      interaction = BINARY_INTERACTION
    )
  }),
  SCENARIOS$scenario_id
)

# =============================================================================
# Allocation-only calibration cache
# =============================================================================

calibration_cache_file <- function(scenario) {
  probs <- factor_probabilities(scenario$factor_count)
  prob_tag <- paste(formatC(probs, format = "f", digits = 3), collapse = "-")
  file.path(
    CALIBRATION_DIR,
    sprintf(
      "binary_siga_pureR_K%d_n%d_p%s_pbc%.3f_B%d.rds",
      scenario$factor_count, scenario$total_n, prob_tag,
      P_BIASED_COIN, N_CALIBRATION
    )
  )
}

get_design_calibration <- function(scenario) {
  cache <- calibration_cache_file(scenario)
  if (file.exists(cache)) return(readRDS(cache))

  lock <- paste0(cache, ".lock")
  have_lock <- dir.create(lock, showWarnings = FALSE)
  if (!have_lock) {
    for (attempt in seq_len(720L)) {
      if (file.exists(cache)) return(readRDS(cache))
      Sys.sleep(5)
      have_lock <- dir.create(lock, showWarnings = FALSE)
      if (have_lock) break
    }
  }
  if (!have_lock) stop("Could not acquire calibration lock: ", lock, call. = FALSE)
  on.exit(unlink(lock, recursive = TRUE, force = TRUE), add = TRUE)

  if (file.exists(cache)) return(readRDS(cache))
  message(
    "Creating allocation-only calibration: K=", scenario$factor_count,
    ", n=", scenario$total_n, ", B0=", N_CALIBRATION
  )
  calibration <- calibrate_siga_binary(
    n = scenario$total_n,
    factor_prob = factor_probabilities(scenario$factor_count),
    pbc = P_BIASED_COIN,
    weights = rep(1, scenario$factor_count + 1L),
    B0 = N_CALIBRATION,
    seed = calibration_seed(scenario$total_n, scenario$factor_count)
  )

  tmp <- tempfile(
    pattern = "binary_calibration_", tmpdir = CALIBRATION_DIR, fileext = ".rds"
  )
  saveRDS(calibration, tmp)
  if (!file.rename(tmp, cache)) {
    if (!file.copy(tmp, cache, overwrite = TRUE)) {
      unlink(tmp)
      stop("Could not save calibration: ", cache, call. = FALSE)
    }
    unlink(tmp)
  }
  calibration
}

# =============================================================================
# One trial and its two analyses
# =============================================================================

simulate_trial <- function(scenario, replicate_id, outcome_model) {
  set.seed(seed_value(scenario$scenario_id, replicate_id, 1L))
  X <- generate_factors(
    n = scenario$total_n,
    factor_prob = factor_probabilities(scenario$factor_count)
  )
  z <- ps_assign_R(
    X = X,
    pbc = P_BIASED_COIN,
    weights_ = rep(1, scenario$factor_count + 1L),
    seed = as.double(seed_value(scenario$scenario_id, replicate_id, 2L))
  )
  A <- as.integer((z + 1L) / 2L)
  set.seed(seed_value(scenario$scenario_id, replicate_id, 3L))
  y <- generate_binary_outcome(X, A, outcome_model)
  list(X = X, z = as.integer(z), A = A, y = y)
}

analyze_proposed <- function(trial, scenario, calibration, scores, spec) {
  nb <- length(spec$boundaries)
  component_u <- rep(NA_real_, nb)
  component_a <- rep(NA_real_, nb)
  z_u <- rep(NA_real_, nb)
  z_a <- rep(NA_real_, nb)

  for (j in seq_along(spec$boundaries)) {
    if (BINARY_ANALYSIS_MODE %in% c("joint", "unadjusted")) {
      use_lattice <- scenario$objective == "superiority" &&
        abs(spec$boundaries[j]) < 1e-14
      result_u <- siga_component_pvalue(
        score = scores$unadjusted[, j],
        z = trial$z,
        X = trial$X,
        calibration = calibration,
        alternative = spec$alternatives[j],
        use_binary_lattice = use_lattice,
        y = trial$y
      )
      component_u[j] <- result_u$p
      z_u[j] <- result_u$z
    }

    if (BINARY_ANALYSIS_MODE %in% c("joint", "adjusted")) {
      result_a <- siga_component_pvalue(
        score = scores$adjusted[, j],
        z = trial$z,
        X = trial$X,
        calibration = calibration,
        alternative = spec$alternatives[j],
        use_binary_lattice = FALSE
      )
      component_a[j] <- result_a$p
      z_a[j] <- result_a$z
    }
  }

  list(
    p_unadjusted = if (BINARY_ANALYSIS_MODE %in% c("joint", "unadjusted")) {
      combine_objective_pvalues(component_u, scenario$objective)
    } else {
      NA_real_
    },
    p_adjusted = if (BINARY_ANALYSIS_MODE %in% c("joint", "adjusted")) {
      combine_objective_pvalues(component_a, scenario$objective)
    } else {
      NA_real_
    },
    component_unadjusted = component_u,
    component_adjusted = component_a,
    z_unadjusted = z_u,
    z_adjusted = z_a
  )
}

analyze_reference <- function(trial, scenario, scores, spec, replicate_id) {
  nb <- length(spec$boundaries)
  component_u <- rep(NA_real_, nb)
  component_a <- rep(NA_real_, nb)

  if (BINARY_ANALYSIS_MODE == "joint") {
    score_matrix <- cbind(scores$unadjusted, scores$adjusted)
  } else if (BINARY_ANALYSIS_MODE == "unadjusted") {
    score_matrix <- scores$unadjusted
  } else {
    score_matrix <- scores$adjusted
  }

  rt <- rt_pvalues_R(
    X = trial$X,
    scores = score_matrix,
    z_obs = trial$z,
    B = N_RERANDOMIZATIONS,
    pbc = P_BIASED_COIN,
    weights_ = rep(1, scenario$factor_count + 1L),
    seed = as.double(seed_value(scenario$scenario_id, replicate_id, 4L))
  )

  for (j in seq_len(nb)) {
    if (BINARY_ANALYSIS_MODE %in% c("joint", "unadjusted")) {
      component_u[j] <- reference_component_pvalue(
        rt, j, spec$alternatives[j]
      )
    }
    if (BINARY_ANALYSIS_MODE == "joint") {
      component_a[j] <- reference_component_pvalue(
        rt, nb + j, spec$alternatives[j]
      )
    } else if (BINARY_ANALYSIS_MODE == "adjusted") {
      component_a[j] <- reference_component_pvalue(
        rt, j, spec$alternatives[j]
      )
    }
  }

  list(
    p_unadjusted = if (BINARY_ANALYSIS_MODE %in% c("joint", "unadjusted")) {
      combine_objective_pvalues(component_u, scenario$objective)
    } else {
      NA_real_
    },
    p_adjusted = if (BINARY_ANALYSIS_MODE %in% c("joint", "adjusted")) {
      combine_objective_pvalues(component_a, scenario$objective)
    } else {
      NA_real_
    },
    component_unadjusted = component_u,
    component_adjusted = component_a
  )
}

component_value <- function(x, position) {
  if (length(x) >= position) x[position] else NA_real_
}

run_replicate <- function(scenario, replicate_id, calibration, outcome_model) {
  generation_start <- proc.time()[3L]
  trial <- simulate_trial(scenario, replicate_id, outcome_model)
  generation_seconds <- proc.time()[3L] - generation_start

  spec <- scenario_specification(scenario)
  score_start <- proc.time()[3L]
  scores <- make_boundary_scores(
    y = trial$y,
    A = trial$A,
    X = trial$X,
    boundaries = spec$boundaries
  )
  score_seconds <- proc.time()[3L] - score_start

  proposed_start <- proc.time()[3L]
  proposed <- analyze_proposed(trial, scenario, calibration, scores, spec)
  proposed_seconds <- proc.time()[3L] - proposed_start

  reference_start <- proc.time()[3L]
  reference <- analyze_reference(trial, scenario, scores, spec, replicate_id)
  reference_seconds <- proc.time()[3L] - reference_start

  data.frame(
    scenario_id = scenario$scenario_id,
    scenario_label = scenario$scenario_label,
    replicate = replicate_id,
    factor_count = scenario$factor_count,
    target_per_group = scenario$target_per_group,
    total_n = scenario$total_n,
    objective = scenario$objective,
    scenario_role = scenario$scenario_role,
    true_risk_difference = scenario$true_risk_difference,
    alpha = scenario$alpha,
    benchmark_version = BENCHMARK_VERSION,
    binary_analysis_mode = BINARY_ANALYSIS_MODE,
    configured_outer_trials = N_OUTER,
    rerandomizations_per_trial = N_RERANDOMIZATIONS,
    allocation_calibration_paths = N_CALIBRATION,
    configured_shards = N_SHARDS,
    achieved_control_risk = outcome_model$achieved_control_risk,
    achieved_treatment_risk = outcome_model$achieved_treatment_risk,
    achieved_risk_difference = outcome_model$achieved_risk_difference,
    p_proposed_unadjusted = proposed$p_unadjusted,
    p_reference_unadjusted = reference$p_unadjusted,
    p_proposed_adjusted = proposed$p_adjusted,
    p_reference_adjusted = reference$p_adjusted,
    p_proposed_unadjusted_lower = component_value(proposed$component_unadjusted, 1L),
    p_proposed_unadjusted_upper = component_value(proposed$component_unadjusted, 2L),
    p_reference_unadjusted_lower = component_value(reference$component_unadjusted, 1L),
    p_reference_unadjusted_upper = component_value(reference$component_unadjusted, 2L),
    p_proposed_adjusted_lower = component_value(proposed$component_adjusted, 1L),
    p_proposed_adjusted_upper = component_value(proposed$component_adjusted, 2L),
    p_reference_adjusted_lower = component_value(reference$component_adjusted, 1L),
    p_reference_adjusted_upper = component_value(reference$component_adjusted, 2L),
    reject_proposed_unadjusted = if (is.na(proposed$p_unadjusted)) NA else proposed$p_unadjusted <= scenario$alpha,
    reject_reference_unadjusted = if (is.na(reference$p_unadjusted)) NA else reference$p_unadjusted <= scenario$alpha,
    reject_proposed_adjusted = if (is.na(proposed$p_adjusted)) NA else proposed$p_adjusted <= scenario$alpha,
    reject_reference_adjusted = if (is.na(reference$p_adjusted)) NA else reference$p_adjusted <= scenario$alpha,
    data_generation_seconds = generation_seconds,
    score_construction_seconds = score_seconds,
    proposed_analysis_seconds = proposed_seconds,
    reference_analysis_seconds = reference_seconds,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# Sharded execution and checkpointing
# =============================================================================

scenario_paths <- function(scenario) {
  directory <- ensure_directory(file.path(
    CHUNK_DIR,
    sprintf("scenario_%02d_%s", scenario$scenario_id, scenario$scenario_label)
  ))
  list(
    directory = directory,
    detail = file.path(
      directory,
      sprintf("shard_%04d_of_%04d.csv", SHARD_ID, N_SHARDS)
    ),
    timing = file.path(
      directory,
      sprintf("shard_%04d_of_%04d_timing.csv", SHARD_ID, N_SHARDS)
    )
  )
}

read_completed_replicates <- function(path) {
  if (!file.exists(path)) return(integer())
  dat <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(dat) || !"replicate" %in% names(dat)) return(integer())
  unique(as.integer(dat$replicate))
}

run_scenario_shard <- function(scenario) {
  paths <- scenario_paths(scenario)
  assigned <- seq.int(from = SHARD_ID, to = N_OUTER, by = N_SHARDS)
  completed <- read_completed_replicates(paths$detail)
  pending <- setdiff(assigned, completed)

  message(
    "Scenario ", scenario$scenario_id, " [", scenario$scenario_label, "]",
    "; shard ", SHARD_ID, "/", N_SHARDS,
    "; assigned=", length(assigned),
    "; completed=", length(completed),
    "; pending=", length(pending)
  )
  if (!length(pending)) return(invisible(NULL))

  calibration <- get_design_calibration(scenario)
  outcome_model <- OUTCOME_MODELS[[as.character(scenario$scenario_id)]]
  batches <- split(pending, ceiling(seq_along(pending) / OUTER_BATCH_SIZE))
  session_start <- proc.time()[3L]

  for (batch_index in seq_along(batches)) {
    ids <- batches[[batch_index]]
    batch_start <- proc.time()[3L]
    worker <- function(id) {
      run_replicate(scenario, id, calibration, outcome_model)
    }
    pieces <- if (USE_PARALLEL) {
      parallel::mclapply(
        ids, worker,
        mc.cores = N_CORES,
        mc.preschedule = FALSE,
        mc.set.seed = FALSE
      )
    } else {
      lapply(ids, worker)
    }
    batch <- do.call(rbind, pieces)
    batch <- batch[order(batch$replicate), , drop = FALSE]
    append_csv(batch, paths$detail)

    message(
      "  batch ", batch_index, "/", length(batches),
      "; replicates ", min(ids), "--", max(ids),
      "; elapsed ", format_elapsed(proc.time()[3L] - batch_start)
    )
  }

  timing <- data.frame(
    scenario_id = scenario$scenario_id,
    scenario_label = scenario$scenario_label,
    shard_id = SHARD_ID,
    n_shards = N_SHARDS,
    assigned_replicates = length(assigned),
    completed_replicates = length(unique(c(completed, pending))),
    wall_clock_seconds = proc.time()[3L] - session_start,
    stringsAsFactors = FALSE
  )
  write.csv(timing, paths$timing, row.names = FALSE)
  invisible(timing)
}

# =============================================================================
# Aggregation: no difference, agreement, or time-ratio columns
# =============================================================================

read_scenario_details <- function(scenario) {
  directory <- file.path(
    CHUNK_DIR,
    sprintf("scenario_%02d_%s", scenario$scenario_id, scenario$scenario_label)
  )
  shard_pattern <- sprintf("^shard_[0-9]{4}_of_%04d\\.csv$", N_SHARDS)
  files <- list.files(
    directory,
    pattern = shard_pattern,
    full.names = TRUE
  )
  if (!length(files)) stop("No shard files found for ", scenario$scenario_label)
  expected_names <- sprintf("shard_%04d_of_%04d.csv", seq_len(N_SHARDS), N_SHARDS)
  observed_names <- basename(files)
  missing_names <- setdiff(expected_names, observed_names)
  if (!ALLOW_PARTIAL && length(missing_names)) {
    stop(
      scenario$scenario_label, " is missing shard file(s): ",
      paste(missing_names, collapse = ", "), call. = FALSE
    )
  }
  pieces <- lapply(files, function(path) read.csv(path, stringsAsFactors = FALSE))
  dat <- do.call(rbind, pieces)
  dat <- dat[order(dat$replicate), , drop = FALSE]

  duplicate_id <- duplicated(dat$replicate)
  if (any(duplicate_id)) {
    duplicated_reps <- unique(dat$replicate[duplicate_id])
    for (id in duplicated_reps) {
      rows <- dat[dat$replicate == id, , drop = FALSE]
      check_cols <- setdiff(names(rows), c(
        "data_generation_seconds", "score_construction_seconds",
        "proposed_analysis_seconds", "reference_analysis_seconds"
      ))
      if (nrow(unique(rows[, check_cols, drop = FALSE])) > 1L) {
        stop("Conflicting duplicate replicate ", id, " in ", scenario$scenario_label)
      }
    }
    dat <- dat[!duplicated(dat$replicate), , drop = FALSE]
  }

  if (!ALLOW_PARTIAL && nrow(dat) != N_OUTER) {
    stop(
      scenario$scenario_label, " has ", nrow(dat),
      " unique replicates; expected ", N_OUTER, "."
    )
  }
  if (!ALLOW_PARTIAL && !setequal(as.integer(dat$replicate), seq_len(N_OUTER))) {
    stop(
      scenario$scenario_label,
      " does not contain exactly replicate IDs 1 through ", N_OUTER, ".",
      call. = FALSE
    )
  }
  config_columns <- c(
    "configured_outer_trials", "rerandomizations_per_trial",
    "allocation_calibration_paths", "configured_shards"
  )
  if (all(config_columns %in% names(dat))) {
    expected_config <- c(N_OUTER, N_RERANDOMIZATIONS, N_CALIBRATION, N_SHARDS)
    for (j in seq_along(config_columns)) {
      nm <- config_columns[j]
      values <- unique(as.numeric(dat[[nm]]))
      if (length(values) != 1L || !is.finite(values) || values != expected_config[j]) {
        stop(
          "Configuration column ", nm, " disagrees in ",
          scenario$scenario_label, call. = FALSE
        )
      }
    }
  }
  dat
}

scenario_summary <- function(dat, scenario) {
  rejection_summary <- function(column) {
    x <- as.logical(dat[[column]])
    x <- x[!is.na(x)]
    if (!length(x)) {
      return(c(probability = NA_real_, lower = NA_real_, upper = NA_real_))
    }
    ci <- wilson_interval(sum(x), length(x))
    c(
      probability = unname(ci["estimate"]),
      lower = unname(ci["lower"]),
      upper = unname(ci["upper"])
    )
  }

  pu <- rejection_summary("reject_proposed_unadjusted")
  ru <- rejection_summary("reject_reference_unadjusted")
  pa <- rejection_summary("reject_proposed_adjusted")
  ra <- rejection_summary("reject_reference_adjusted")

  data.frame(
    scenario_id = scenario$scenario_id,
    scenario_label = scenario$scenario_label,
    factor_count = scenario$factor_count,
    target_per_group = scenario$target_per_group,
    total_n = scenario$total_n,
    objective = scenario$objective,
    scenario_role = scenario$scenario_role,
    true_risk_difference = scenario$true_risk_difference,
    alpha = scenario$alpha,
    benchmark_version = BENCHMARK_VERSION,
    n_outer = nrow(dat),
    rerandomizations_per_trial = N_RERANDOMIZATIONS,
    allocation_calibration_paths = N_CALIBRATION,
    proposed_unadjusted = unname(pu["probability"]),
    proposed_unadjusted_lower_95 = unname(pu["lower"]),
    proposed_unadjusted_upper_95 = unname(pu["upper"]),
    reference_unadjusted = unname(ru["probability"]),
    reference_unadjusted_lower_95 = unname(ru["lower"]),
    reference_unadjusted_upper_95 = unname(ru["upper"]),
    proposed_adjusted = unname(pa["probability"]),
    proposed_adjusted_lower_95 = unname(pa["lower"]),
    proposed_adjusted_upper_95 = unname(pa["upper"]),
    reference_adjusted = unname(ra["probability"]),
    reference_adjusted_lower_95 = unname(ra["lower"]),
    reference_adjusted_upper_95 = unname(ra["upper"]),
    score_construction_minutes = sum(dat$score_construction_seconds) / 60,
    proposed_analysis_minutes = sum(dat$proposed_analysis_seconds) / 60,
    reference_analysis_minutes = sum(dat$reference_analysis_seconds) / 60,
    stringsAsFactors = FALSE
  )
}

objective_label <- function(x) {
  switch(
    x,
    superiority = "Superiority",
    noninferiority = "Non-inferiority",
    equivalence = "Equivalence"
  )
}

role_label <- function(x) {
  switch(
    x,
    type1 = "Type I error",
    type1_lower = "Type I error at lower limit",
    type1_upper = "Type I error at upper limit",
    power = "Power"
  )
}

write_latex_results <- function(summary, path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\small",
    "\\begin{tabular}{rrlllrrrr}",
    "\\toprule",
    "Factors & $n$/group & Objective & Scenario & True RD & Proposed U (\\%) & RT U (\\%) & Proposed A (\\%) & RT A (\\%)\\\\",
    "\\midrule"
  ), con)
  for (i in seq_len(nrow(summary))) {
    row <- summary[i, ]
    line <- sprintf(
      "%d & %d & %s & %s & %.3f & %.2f & %.2f & %.2f & %.2f\\\\",
      row$factor_count,
      row$target_per_group,
      objective_label(row$objective),
      role_label(row$scenario_role),
      row$true_risk_difference,
      100 * row$proposed_unadjusted,
      100 * row$reference_unadjusted,
      100 * row$proposed_adjusted,
      100 * row$reference_adjusted
    )
    writeLines(line, con)
  }
  writeLines(c(
    "\\bottomrule",
    "\\end{tabular}",
    paste0("\\caption{Binary operating characteristics of SIGA and the reference randomization test based on $M=", format(N_OUTER, big.mark = ",", scientific = FALSE), "$ outer trials, $B=", format(N_RERANDOMIZATIONS, big.mark = ",", scientific = FALSE), "$ regenerated allocation paths per reference test, and $B_0=", format(N_CALIBRATION, big.mark = ",", scientific = FALSE), "$ allocation-only calibration paths. U and A denote unadjusted and baseline-adjusted scores, respectively.}"),
    "\\label{tab:binary-siga-results}",
    "\\end{table}"
  ), con)
}

write_latex_timing <- function(summary, path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\small",
    "\\begin{tabular}{rrlllrr}",
    "\\toprule",
    "Factors & $n$/group & Objective & Scenario & True RD & Proposed time (min) & RT time (min)\\\\",
    "\\midrule"
  ), con)
  for (i in seq_len(nrow(summary))) {
    row <- summary[i, ]
    line <- sprintf(
      "%d & %d & %s & %s & %.3f & %.2f & %.2f\\\\",
      row$factor_count,
      row$target_per_group,
      objective_label(row$objective),
      role_label(row$scenario_role),
      row$true_risk_difference,
      row$proposed_analysis_minutes,
      row$reference_analysis_minutes
    )
    writeLines(line, con)
  }
  writeLines(c(
    "\\bottomrule",
    "\\end{tabular}",
    paste0("\\caption{Cumulative analysis time for the binary benchmark over $M=", format(N_OUTER, big.mark = ",", scientific = FALSE), "$ outer trials. Each reference test used $B=", format(N_RERANDOMIZATIONS, big.mark = ",", scientific = FALSE), "$ regenerated allocation paths. Times exclude data generation and score construction and are reported in minutes.}"),
    "\\label{tab:binary-siga-timing}",
    "\\end{table}"
  ), con)
}

write_simulation_configuration <- function(path_csv, path_tex) {
  configuration <- data.frame(
    benchmark_version = BENCHMARK_VERSION,
    outer_trials_per_scenario = N_OUTER,
    regenerated_allocation_paths_per_reference_test = N_RERANDOMIZATIONS,
    allocation_only_calibration_paths = N_CALIBRATION,
    biased_coin_probability = P_BIASED_COIN,
    number_of_shards = N_SHARDS,
    scenario_ids = paste(SCENARIOS$scenario_id, collapse = ","),
    reference_statistic = "fixed boundary-specific treatment score",
    reference_pvalue = "inclusive plus-one Monte Carlo p-value",
    score_recomputed_for_rerandomization = FALSE,
    binary_analysis_mode = BINARY_ANALYSIS_MODE,
    stringsAsFactors = FALSE
  )
  write.csv(configuration, path_csv, row.names = FALSE)

  con <- file(path_tex, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\small",
    "\\begin{tabular}{lr}",
    "\\toprule",
    "Configuration item & Value\\\\",
    "\\midrule",
    sprintf("Outer trials per scenario, $M$ & %s\\\\", format(N_OUTER, big.mark = ",", scientific = FALSE)),
    sprintf("Regenerated allocation paths per reference test, $B$ & %s\\\\", format(N_RERANDOMIZATIONS, big.mark = ",", scientific = FALSE)),
    sprintf("Allocation-only calibration paths, $B_0$ & %s\\\\", format(N_CALIBRATION, big.mark = ",", scientific = FALSE)),
    sprintf("Biased-coin probability & %.2f\\\\", P_BIASED_COIN),
    sprintf("Computational shards & %d\\\\", N_SHARDS),
    "\\bottomrule",
    "\\end{tabular}",
    "\\caption{Final binary simulation configuration. The reference randomization test uses the same fixed-score rerandomization algorithm as the continuous-outcome benchmark.}",
    "\\label{tab:binary-simulation-configuration}",
    "\\end{table}"
  ), con)
  invisible(configuration)
}

aggregate_results <- function() {
  summaries <- vector("list", nrow(SCENARIOS))
  for (i in seq_len(nrow(SCENARIOS))) {
    scenario <- SCENARIOS[i, ]
    message("Aggregating ", scenario$scenario_label)
    dat <- read_scenario_details(scenario)
    summaries[[i]] <- scenario_summary(dat, scenario)
  }
  summary <- do.call(rbind, summaries)
  summary <- summary[order(summary$scenario_id), , drop = FALSE]

  csv_path <- file.path(OUTPUT_DIR, "binary_siga_weak_null_summary.csv")
  tex_path <- file.path(OUTPUT_DIR, "binary_siga_weak_null_results.tex")
  timing_tex_path <- file.path(OUTPUT_DIR, "binary_siga_weak_null_timing.tex")
  config_csv_path <- file.path(OUTPUT_DIR, "binary_simulation_configuration.csv")
  config_tex_path <- file.path(OUTPUT_DIR, "binary_simulation_configuration.tex")
  write.csv(summary, csv_path, row.names = FALSE)
  write_latex_results(summary, tex_path)
  write_latex_timing(summary, timing_tex_path)
  write_simulation_configuration(config_csv_path, config_tex_path)

  message("Wrote: ", csv_path)
  message("Wrote: ", tex_path)
  message("Wrote: ", timing_tex_path)
  invisible(summary)
}

# =============================================================================
# Main
# =============================================================================

if (RUN_MODE == "run") {
  for (i in seq_len(nrow(SCENARIOS))) {
    run_scenario_shard(SCENARIOS[i, ])
  }
} else {
  aggregate_results()
}
