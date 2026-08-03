# Core functions for the additional SIGA pair-path simulations.
# Base R only. This file defines functions and does not start a simulation.

options(stringsAsFactors = FALSE, warn = 1)

`%||%` <- function(x, y) if (is.null(x)) y else x

normalise_seed <- function(seed) {
  seed <- as.double(seed)
  if (!is.finite(seed)) stop("seed must be finite.", call. = FALSE)
  as.integer(abs(seed) %% 2147483646 + 1)
}

check_binary_factor_matrix <- function(X) {
  X <- as.matrix(X)
  storage.mode(X) <- "integer"
  if (anyNA(X) || any(!X %in% c(0L, 1L))) {
    stop("All minimization factors must be coded 0/1.", call. = FALSE)
  }
  X
}

all_binary_patterns <- function(K) {
  K <- as.integer(K)
  if (K < 1L || K > 20L) stop("K must be between 1 and 20.", call. = FALSE)
  out <- as.matrix(expand.grid(
    rep(list(c(0L, 1L)), K),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ))
  storage.mode(out) <- "integer"
  colnames(out) <- paste0("X", seq_len(K))
  out
}

joint_stratum_id <- function(X) {
  X <- check_binary_factor_matrix(X)
  K <- ncol(X)
  as.integer(1L + drop(X %*% (2L ^ (0:(K - 1L)))))
}

nearest_psd <- function(M, tolerance = 1e-12) {
  M <- as.matrix(M)
  M <- (M + t(M)) / 2
  ee <- eigen(M, symmetric = TRUE)
  cutoff <- tolerance * max(1, max(abs(ee$values)))
  values <- pmax(ee$values, 0)
  values[values < cutoff] <- 0
  out <- ee$vectors %*% (values * t(ee$vectors))
  (out + t(out)) / 2
}

profile_prob_independent <- function(factor_prob, patterns = NULL) {
  factor_prob <- as.numeric(factor_prob)
  if (any(!is.finite(factor_prob)) || any(factor_prob <= 0 | factor_prob >= 1)) {
    stop("Each factor probability must be strictly between 0 and 1.", call. = FALSE)
  }
  patterns <- patterns %||% all_binary_patterns(length(factor_prob))
  prob <- apply(patterns, 1L, function(x) {
    prod(ifelse(x == 1L, factor_prob, 1 - factor_prob))
  })
  prob / sum(prob)
}

profile_prob_latent_correlated <- function(factor_prob,
                                           latent_strength = 1.0,
                                           mixing_prob = 0.5,
                                           patterns = NULL) {
  factor_prob <- as.numeric(factor_prob)
  if (any(!is.finite(factor_prob)) || any(factor_prob <= 0 | factor_prob >= 1)) {
    stop("Each factor probability must be strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.finite(latent_strength) || latent_strength < 0) {
    stop("latent_strength must be finite and nonnegative.", call. = FALSE)
  }
  if (!is.finite(mixing_prob) || mixing_prob <= 0 || mixing_prob >= 1) {
    stop("mixing_prob must be strictly between 0 and 1.", call. = FALSE)
  }

  K <- length(factor_prob)
  patterns <- patterns %||% all_binary_patterns(K)
  intercept <- numeric(K)
  q0 <- numeric(K)
  q1 <- numeric(K)

  for (j in seq_len(K)) {
    target <- factor_prob[j]
    objective <- function(a) {
      (1 - mixing_prob) * plogis(a - latent_strength) +
        mixing_prob * plogis(a + latent_strength) - target
    }
    intercept[j] <- uniroot(objective, interval = c(-30, 30), tol = 1e-13)$root
    q0[j] <- plogis(intercept[j] - latent_strength)
    q1[j] <- plogis(intercept[j] + latent_strength)
  }

  component_prob <- function(q) {
    apply(patterns, 1L, function(x) prod(ifelse(x == 1L, q, 1 - q)))
  }
  prob <- (1 - mixing_prob) * component_prob(q0) + mixing_prob * component_prob(q1)
  prob <- prob / sum(prob)

  attr(prob, "latent_intercept") <- intercept
  attr(prob, "conditional_prob_0") <- q0
  attr(prob, "conditional_prob_1") <- q1
  prob
}

profile_marginals <- function(patterns, profile_prob) {
  patterns <- check_binary_factor_matrix(patterns)
  profile_prob <- as.numeric(profile_prob)
  drop(crossprod(profile_prob, patterns))
}

profile_correlations <- function(patterns, profile_prob) {
  patterns <- check_binary_factor_matrix(patterns)
  profile_prob <- as.numeric(profile_prob)
  K <- ncol(patterns)
  mu <- profile_marginals(patterns, profile_prob)
  second <- crossprod(patterns * sqrt(profile_prob), patterns * sqrt(profile_prob))
  covariance <- second - tcrossprod(mu)
  sdv <- sqrt(pmax(diag(covariance), 0))
  denominator <- outer(sdv, sdv)
  correlation <- covariance / denominator
  correlation[!is.finite(correlation)] <- 0
  diag(correlation) <- 1
  correlation
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

rt_pvalues_R <- function(X,
                         scores,
                         z_obs,
                         B = 1999L,
                         pbc = 0.80,
                         weights_ = NULL,
                         seed = 1,
                         tolerance = 1e-12,
                         return_randomization_moments = TRUE) {
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

  out <- list(
    statistic = t_obs,
    two_sided = p_two,
    greater = p_upper,
    less = p_lower,
    B = B,
    engine = "pure_R_vectorized"
  )
  if (isTRUE(return_randomization_moments)) {
    out$randomization_mean <- colMeans(t_rand)
    out$randomization_variance <- if (B > 1L) {
      apply(t_rand, 2L, var)
    } else {
      rep(NA_real_, L)
    }
  }
  out
}

covariance_from_sums <- function(sum_x, sum_xx, n) {
  if (n < 2L) stop("At least two observations are required.", call. = FALSE)
  (sum_xx - tcrossprod(sum_x) / n) / (n - 1L)
}

cross_covariance_from_sums <- function(sum_x, sum_y, sum_xy, n) {
  if (n < 2L) stop("At least two observations are required.", call. = FALSE)
  (sum_xy - tcrossprod(sum_x, sum_y) / n) / (n - 1L)
}

calibrate_pair_path_design_R <- function(B0,
                                         n,
                                         patterns,
                                         profile_prob,
                                         pbc = 0.80,
                                         weights_ = NULL,
                                         seed = 1,
                                         batch_size = 1000L,
                                         progress = TRUE) {
  B0 <- as.integer(B0)
  n <- as.integer(n)
  patterns <- check_binary_factor_matrix(patterns)
  profile_prob <- as.numeric(profile_prob)
  K <- ncol(patterns)
  J <- nrow(patterns)
  weights <- as.numeric(weights_ %||% rep(1, K + 1L))

  if (B0 < 2L) stop("B0 must be at least 2.", call. = FALSE)
  if (n < 2L) stop("n must be at least 2.", call. = FALSE)
  if (J != 2L ^ K) stop("patterns must contain all 2^K binary profiles.", call. = FALSE)
  if (length(profile_prob) != J || any(!is.finite(profile_prob)) ||
      any(profile_prob <= 0)) {
    stop("profile_prob must contain positive finite probabilities for all profiles.", call. = FALSE)
  }
  profile_prob <- profile_prob / sum(profile_prob)
  validate_minimization_inputs(K, pbc, weights)
  batch_size <- max(1L, min(as.integer(batch_size), B0))
  set.seed(normalise_seed(seed))

  sum_u <- numeric(J)
  sum_uu <- matrix(0, J, J)

  sum_d <- lapply(seq_len(3L), function(i) numeric(J))
  sum_dd <- lapply(seq_len(3L), function(i) matrix(0, J, J))
  sum_d01 <- matrix(0, J, J)
  sum_d02 <- matrix(0, J, J)
  sum_d12 <- matrix(0, J, J)

  sum_g01 <- numeric(J)
  sum_g02 <- numeric(J)
  sum_g01g01 <- matrix(0, J, J)
  sum_g02g02 <- matrix(0, J, J)
  sum_g01g02 <- matrix(0, J, J)

  treated_count <- numeric(n + 1L)
  completed <- 0L
  start <- proc.time()[3L]

  while (completed < B0) {
    m <- min(batch_size, B0 - completed)
    rows <- seq_len(m)
    overall <- lapply(seq_len(3L), function(i) integer(m))
    marginal <- lapply(seq_len(3L), function(i) {
      matrix(0L, nrow = m, ncol = 2L * K)
    })
    counts <- matrix(0L, nrow = m, ncol = J)
    imbalance <- lapply(seq_len(3L), function(i) {
      matrix(0L, nrow = m, ncol = J)
    })
    pair01 <- matrix(0L, nrow = m, ncol = J)
    pair02 <- matrix(0L, nrow = m, ncol = J)

    for (i in seq_len(n)) {
      sid <- sample.int(J, size = m, replace = TRUE, prob = profile_prob)
      X <- patterns[sid, , drop = FALSE]
      cols <- matrix(0L, nrow = m, ncol = K)
      for (j in seq_len(K)) {
        cols[, j] <- 2L * (j - 1L) + X[, j] + 1L
      }

      z_copy <- vector("list", 3L)
      for (copy in seq_len(3L)) {
        score_plus <- weights[1L] * abs(overall[[copy]] + 1L)
        score_minus <- weights[1L] * abs(overall[[copy]] - 1L)
        for (j in seq_len(K)) {
          idx_j <- rows + (cols[, j] - 1L) * m
          current <- marginal[[copy]][idx_j]
          score_plus <- score_plus + weights[j + 1L] * abs(current + 1L)
          score_minus <- score_minus + weights[j + 1L] * abs(current - 1L)
        }

        u <- runif(m)
        z <- integer(m)
        better_plus <- score_plus < score_minus
        better_minus <- score_plus > score_minus
        ties <- !(better_plus | better_minus)
        z[better_plus] <- ifelse(u[better_plus] < pbc, 1L, -1L)
        z[better_minus] <- ifelse(u[better_minus] < 1 - pbc, 1L, -1L)
        z[ties] <- ifelse(u[ties] < 0.5, 1L, -1L)

        z_copy[[copy]] <- z
        overall[[copy]] <- overall[[copy]] + z
        for (j in seq_len(K)) {
          idx_j <- rows + (cols[, j] - 1L) * m
          marginal[[copy]][idx_j] <- marginal[[copy]][idx_j] + z
        }
      }

      idx_s <- rows + (sid - 1L) * m
      counts[idx_s] <- counts[idx_s] + 1L
      for (copy in seq_len(3L)) {
        imbalance[[copy]][idx_s] <- imbalance[[copy]][idx_s] + z_copy[[copy]]
      }
      pair01[idx_s] <- pair01[idx_s] + z_copy[[1L]] * z_copy[[2L]]
      pair02[idx_s] <- pair02[idx_s] + z_copy[[1L]] * z_copy[[3L]]
    }

    U <- matrix(0, nrow = m, ncol = J)
    nonzero <- counts > 0L
    U[nonzero] <- imbalance[[1L]][nonzero] / sqrt(counts[nonzero])
    D <- lapply(imbalance, function(x) x / sqrt(n))
    G01 <- pair01 / sqrt(n)
    G02 <- pair02 / sqrt(n)

    sum_u <- sum_u + colSums(U)
    sum_uu <- sum_uu + crossprod(U)
    for (copy in seq_len(3L)) {
      sum_d[[copy]] <- sum_d[[copy]] + colSums(D[[copy]])
      sum_dd[[copy]] <- sum_dd[[copy]] + crossprod(D[[copy]])
    }
    sum_d01 <- sum_d01 + crossprod(D[[1L]], D[[2L]])
    sum_d02 <- sum_d02 + crossprod(D[[1L]], D[[3L]])
    sum_d12 <- sum_d12 + crossprod(D[[2L]], D[[3L]])

    sum_g01 <- sum_g01 + colSums(G01)
    sum_g02 <- sum_g02 + colSums(G02)
    sum_g01g01 <- sum_g01g01 + crossprod(G01)
    sum_g02g02 <- sum_g02g02 + crossprod(G02)
    sum_g01g02 <- sum_g01g02 + crossprod(G01, G02)

    treated <- as.integer((n + overall[[1L]]) / 2L)
    treated_count <- treated_count + tabulate(treated + 1L, nbins = n + 1L)

    completed <- completed + m
    if (isTRUE(progress)) {
      message(
        "  three-copy calibration: ", completed, "/", B0,
        " paths; elapsed ", sprintf("%.1f min", (proc.time()[3L] - start) / 60)
      )
    }
  }

  gamma <- covariance_from_sums(sum_u, sum_uu, B0)
  sigma_d <- lapply(seq_len(3L), function(copy) {
    covariance_from_sums(sum_d[[copy]], sum_dd[[copy]], B0)
  })
  psi01 <- covariance_from_sums(sum_g01, sum_g01g01, B0)
  psi02 <- covariance_from_sums(sum_g02, sum_g02g02, B0)
  psi <- nearest_psd((psi01 + psi02) / 2)

  cross_d01 <- cross_covariance_from_sums(sum_d[[1L]], sum_d[[2L]], sum_d01, B0)
  cross_d02 <- cross_covariance_from_sums(sum_d[[1L]], sum_d[[3L]], sum_d02, B0)
  cross_d12 <- cross_covariance_from_sums(sum_d[[2L]], sum_d[[3L]], sum_d12, B0)
  cross_g0102 <- cross_covariance_from_sums(sum_g01, sum_g02, sum_g01g02, B0)

  structure(
    list(
      gamma = nearest_psd(gamma),
      mean_u = sum_u / B0,
      sigma_d = lapply(sigma_d, nearest_psd),
      psi = psi,
      psi01 = nearest_psd(psi01),
      psi02 = nearest_psd(psi02),
      mean_g01 = sum_g01 / B0,
      mean_g02 = sum_g02 / B0,
      cross_d01 = cross_d01,
      cross_d02 = cross_d02,
      cross_d12 = cross_d12,
      cross_g0102 = cross_g0102,
      treated_count_prob = treated_count / B0,
      profile_prob = profile_prob,
      patterns = patterns,
      B0 = B0,
      n = n,
      K = K,
      J = J,
      pbc = pbc,
      weights = weights,
      seed = seed,
      elapsed_seconds = proc.time()[3L] - start,
      engine = "pure_R_three_copy_vectorized"
    ),
    class = c("siga_pair_calibration", "list")
  )
}

null_space_basis_for_profile_mean <- function(profile_prob) {
  profile_prob <- as.numeric(profile_prob)
  J <- length(profile_prob)
  if (J < 2L || any(!is.finite(profile_prob)) || any(profile_prob <= 0)) {
    stop("profile_prob must contain at least two positive values.", call. = FALSE)
  }
  qfit <- qr(matrix(profile_prob, ncol = 1L), LAPACK = FALSE)
  qfull <- qr.Q(qfit, complete = TRUE)
  qfull[, -1L, drop = FALSE]
}

generalized_pair_directions <- function(psi, profile_prob) {
  psi <- nearest_psd(psi)
  profile_prob <- as.numeric(profile_prob)
  J <- length(profile_prob)
  if (!all(dim(psi) == c(J, J))) {
    stop("psi and profile_prob dimensions do not agree.", call. = FALSE)
  }

  Q <- null_space_basis_for_profile_mean(profile_prob)
  Pi <- diag(profile_prob, nrow = J, ncol = J)
  A <- crossprod(Q, psi %*% Q)
  C <- crossprod(Q, Pi %*% Q)
  R <- chol((C + t(C)) / 2)
  Rinv <- backsolve(R, diag(ncol(R)))
  M <- crossprod(Rinv, A %*% Rinv)
  ee <- eigen((M + t(M)) / 2, symmetric = TRUE)
  order_index <- order(ee$values)

  make_direction <- function(index) {
    y <- ee$vectors[, index]
    x <- drop(Rinv %*% y)
    d <- drop(Q %*% x)
    d <- d / max(abs(d))
    d <- d - sum(profile_prob * d) / sum(profile_prob)
    d <- d / max(abs(d))
    numerator <- drop(crossprod(d, psi %*% d))
    denominator <- sum(profile_prob * d^2)
    list(
      direction = d,
      ratio = numerator / denominator,
      numerator = numerator,
      denominator = denominator
    )
  }

  min_result <- make_direction(order_index[1L])
  max_result <- make_direction(order_index[length(order_index)])
  list(
    min_ratio = min_result$ratio,
    min_direction = min_result$direction,
    max_ratio = max_result$ratio,
    max_direction = max_result$direction,
    all_ratios = sort(ee$values),
    basis = Q
  )
}

scale_direction <- function(direction, target_max_abs = 1) {
  direction <- as.numeric(direction)
  if (!length(direction) || any(!is.finite(direction))) {
    stop("direction must be a finite nonempty vector.", call. = FALSE)
  }
  m <- max(abs(direction))
  if (m <= 0) return(rep(0, length(direction)))
  direction * (target_max_abs / m)
}

scale_direction_for_binary_delta <- function(direction,
                                             boundary,
                                             target_max_abs_d = 0.80,
                                             max_abs_delta = 0.90,
                                             safety = 0.98) {
  direction <- scale_direction(direction, 1)
  upper <- max_abs_delta
  lower <- -max_abs_delta
  candidate <- Inf
  positive <- direction > 0
  negative <- direction < 0
  if (any(positive)) {
    candidate <- min(candidate, min((upper - boundary) / direction[positive]))
  }
  if (any(negative)) {
    candidate <- min(candidate, min((lower - boundary) / direction[negative]))
  }
  if (!is.finite(candidate) || candidate <= 0) {
    stop("No positive heterogeneity scale satisfies the binary bounds.", call. = FALSE)
  }
  scale <- min(target_max_abs_d, safety * candidate)
  d <- scale * direction
  delta <- boundary + d
  if (any(abs(delta) >= 1)) {
    stop("Scaled binary risk differences are outside (-1,1).", call. = FALSE)
  }
  d
}

weighted_center <- function(x, weight) {
  x <- as.numeric(x)
  weight <- as.numeric(weight)
  x - sum(weight * x) / sum(weight)
}

pattern_baseline_mean <- function(patterns,
                                  profile_prob,
                                  beta_master = c(0.50, 0.40, 0.30, 0.20, 0.10),
                                  interaction = 0.25) {
  patterns <- check_binary_factor_matrix(patterns)
  K <- ncol(patterns)
  beta <- as.numeric(beta_master)[seq_len(K)]
  mu <- drop(patterns %*% beta)
  if (K >= 2L) mu <- mu + interaction * patterns[, 1L] * patterns[, 2L]
  weighted_center(mu, profile_prob)
}

calibrate_logistic_baseline <- function(patterns,
                                        profile_prob,
                                        target_control_risk = 0.60,
                                        beta_master = c(0.35, -0.25, 0.20, -0.15, 0.10),
                                        interaction = 0.20) {
  patterns <- check_binary_factor_matrix(patterns)
  profile_prob <- as.numeric(profile_prob)
  K <- ncol(patterns)
  beta <- as.numeric(beta_master)[seq_len(K)]
  eta_no_intercept <- drop(patterns %*% beta)
  if (K >= 2L) {
    eta_no_intercept <- eta_no_intercept + interaction * patterns[, 1L] * patterns[, 2L]
  }
  objective <- function(alpha0) {
    sum(profile_prob * plogis(alpha0 + eta_no_intercept)) - target_control_risk
  }
  alpha0 <- uniroot(objective, interval = c(-30, 30), tol = 1e-13)$root
  list(
    intercept = alpha0,
    eta_no_intercept = eta_no_intercept,
    p0 = plogis(alpha0 + eta_no_intercept),
    achieved_control_risk = sum(profile_prob * plogis(alpha0 + eta_no_intercept))
  )
}

calibrate_common_log_odds_model <- function(patterns,
                                            profile_prob,
                                            target_control_risk,
                                            target_risk_difference,
                                            beta_master = c(0.35, -0.25, 0.20, -0.15, 0.10),
                                            interaction = 0.20) {
  baseline <- calibrate_logistic_baseline(
    patterns = patterns,
    profile_prob = profile_prob,
    target_control_risk = target_control_risk,
    beta_master = beta_master,
    interaction = interaction
  )
  objective <- function(theta) {
    p1 <- plogis(baseline$intercept + baseline$eta_no_intercept + theta)
    sum(profile_prob * (p1 - baseline$p0)) - target_risk_difference
  }
  theta <- uniroot(objective, interval = c(-30, 30), tol = 1e-13)$root
  p1 <- plogis(baseline$intercept + baseline$eta_no_intercept + theta)
  list(
    type = "binary_common_log_odds",
    p0 = baseline$p0,
    p1 = p1,
    theta = theta,
    d = p1 - baseline$p0 - target_risk_difference,
    achieved_control_risk = sum(profile_prob * baseline$p0),
    achieved_treatment_risk = sum(profile_prob * p1),
    achieved_effect = sum(profile_prob * (p1 - baseline$p0))
  )
}

make_continuous_homogeneous_model <- function(patterns,
                                              profile_prob,
                                              boundary,
                                              outcome_sd = 1.0) {
  mu0 <- pattern_baseline_mean(patterns, profile_prob)
  list(
    type = "continuous_realistic",
    boundary = boundary,
    d = rep(0, nrow(patterns)),
    delta = rep(boundary, nrow(patterns)),
    mu0 = mu0,
    outcome_sd = outcome_sd,
    individual_effect_sd = 0,
    achieved_effect = boundary
  )
}

make_continuous_pure_model <- function(patterns,
                                       profile_prob,
                                       boundary,
                                       direction,
                                       target_max_abs_d = 1.0,
                                       shared_noise_sd = 0.05) {
  d <- scale_direction(direction, target_max_abs_d)
  d <- weighted_center(d, profile_prob)
  list(
    type = "continuous_symmetric",
    boundary = boundary,
    d = d,
    delta = boundary + d,
    shared_noise_sd = shared_noise_sd,
    achieved_effect = sum(profile_prob * (boundary + d))
  )
}

make_continuous_realistic_model <- function(patterns,
                                            profile_prob,
                                            boundary,
                                            direction,
                                            target_max_abs_d = 0.50,
                                            outcome_sd = 1.0,
                                            individual_effect_sd = 0.25) {
  d <- scale_direction(direction, target_max_abs_d)
  d <- weighted_center(d, profile_prob)
  mu0 <- pattern_baseline_mean(patterns, profile_prob)
  list(
    type = "continuous_realistic",
    boundary = boundary,
    d = d,
    delta = boundary + d,
    mu0 = mu0,
    outcome_sd = outcome_sd,
    individual_effect_sd = individual_effect_sd,
    achieved_effect = sum(profile_prob * (boundary + d))
  )
}

make_binary_homogeneous_model <- function(patterns,
                                          profile_prob,
                                          boundary,
                                          target_control_risk = 0.60) {
  baseline <- calibrate_logistic_baseline(
    patterns = patterns,
    profile_prob = profile_prob,
    target_control_risk = target_control_risk
  )
  p1 <- baseline$p0 + boundary
  if (any(p1 <= 0 | p1 >= 1)) {
    stop("The constant risk-difference model produces invalid treatment risks.", call. = FALSE)
  }
  list(
    type = "binary_direct_probability",
    boundary = boundary,
    p0 = baseline$p0,
    p1 = p1,
    d = rep(0, nrow(patterns)),
    achieved_control_risk = sum(profile_prob * baseline$p0),
    achieved_treatment_risk = sum(profile_prob * p1),
    achieved_effect = sum(profile_prob * (p1 - baseline$p0))
  )
}

make_binary_symmetric_model <- function(patterns,
                                        profile_prob,
                                        boundary,
                                        direction,
                                        target_max_abs_d = 0.80,
                                        max_abs_delta = 0.90) {
  d <- scale_direction_for_binary_delta(
    direction = direction,
    boundary = boundary,
    target_max_abs_d = target_max_abs_d,
    max_abs_delta = max_abs_delta
  )
  d <- weighted_center(d, profile_prob)
  delta <- boundary + d
  p0 <- 0.5 - 0.5 * delta
  p1 <- 0.5 + 0.5 * delta
  if (any(p0 <= 0 | p0 >= 1 | p1 <= 0 | p1 >= 1)) {
    stop("The symmetric binary stress model produces invalid probabilities.", call. = FALSE)
  }
  list(
    type = "binary_direct_probability",
    boundary = boundary,
    p0 = p0,
    p1 = p1,
    d = p1 - p0 - boundary,
    achieved_control_risk = sum(profile_prob * p0),
    achieved_treatment_risk = sum(profile_prob * p1),
    achieved_effect = sum(profile_prob * (p1 - p0))
  )
}

generate_profile_sequence <- function(n, patterns, profile_prob) {
  sid <- sample.int(nrow(patterns), size = n, replace = TRUE, prob = profile_prob)
  list(id = sid, X = patterns[sid, , drop = FALSE])
}

generate_outcome_from_model <- function(profile_id, A, z, model) {
  profile_id <- as.integer(profile_id)
  A <- as.integer(A)
  z <- as.integer(z)
  n <- length(A)
  if (length(profile_id) != n || length(z) != n) {
    stop("profile_id, A and z lengths do not agree.", call. = FALSE)
  }

  if (model$type == "continuous_symmetric") {
    epsilon <- rnorm(n, sd = model$shared_noise_sd)
    return(epsilon + 0.5 * z * model$delta[profile_id])
  }
  if (model$type == "continuous_realistic") {
    y0 <- model$mu0[profile_id] + rnorm(n, sd = model$outcome_sd)
    eta <- if (model$individual_effect_sd > 0) {
      rnorm(n, sd = model$individual_effect_sd)
    } else {
      rep(0, n)
    }
    return(y0 + A * (model$delta[profile_id] + eta))
  }
  if (model$type == "binary_direct_probability" ||
      model$type == "binary_common_log_odds") {
    p <- ifelse(A == 1L, model$p1[profile_id], model$p0[profile_id])
    return(rbinom(n, size = 1L, prob = p))
  }
  stop("Unknown outcome model type: ", model$type, call. = FALSE)
}

make_boundary_scores <- function(y, A, X, boundary) {
  y <- as.numeric(y)
  A <- as.numeric(A)
  X <- check_binary_factor_matrix(X)
  if (length(y) != length(A) || length(y) != nrow(X)) {
    stop("y, A and X must contain the same number of participants.", call. = FALSE)
  }
  W <- y - boundary * A
  unadjusted <- W - mean(W)
  C <- cbind(`(Intercept)` = 1, X)
  adjusted <- qr.resid(qr(C), W)
  adjusted <- adjusted - mean(adjusted)
  list(W = W, unadjusted = unadjusted, adjusted = adjusted)
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

siga_sampling_variance <- function(score,
                                   X,
                                   calibration,
                                   truncate_kappa = TRUE) {
  if (!inherits(calibration, "siga_pair_calibration")) {
    stop("calibration must be produced by calibrate_pair_path_design_R().", call. = FALSE)
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
  kappa_raw <- if (denominator > 0L) {
    (calibration$n - trace_low) / denominator
  } else {
    0
  }
  kappa <- if (truncate_kappa) max(kappa_raw, 0) else kappa_raw

  low <- drop(crossprod(dec$means, Omega %*% dec$means))
  within <- sum(dec$residual^2)
  variance <- 0.25 * (low + kappa * within)
  if (!is.finite(variance) || variance <= 0) variance <- .Machine$double.eps

  list(
    variance = variance,
    kappa = kappa,
    kappa_raw = kappa_raw,
    low_component = 0.25 * low,
    within_component = 0.25 * kappa * within,
    decomposition = dec
  )
}

siga_randomization_variance <- function(sampling_variance_result,
                                        d,
                                        calibration) {
  d <- as.numeric(d)
  if (length(d) != calibration$J) {
    stop("d must have one entry per joint profile.", call. = FALSE)
  }
  counts <- sampling_variance_result$decomposition$counts
  pi_hat <- counts / sum(counts)
  correction_matrix <- calibration$psi - diag(pi_hat, nrow = length(pi_hat))
  scalar_gap <- drop(crossprod(d, correction_matrix %*% d))
  correction <- calibration$n * scalar_gap / 16
  variance_raw <- sampling_variance_result$variance + correction
  variance <- if (is.finite(variance_raw) && variance_raw > 0) {
    variance_raw
  } else {
    .Machine$double.eps
  }
  list(
    variance = variance,
    variance_raw = variance_raw,
    correction = correction,
    scalar_gap = scalar_gap,
    pi_hat = pi_hat,
    truncated = !is.finite(variance_raw) || variance_raw <= 0
  )
}

gaussian_pvalue <- function(statistic,
                            variance,
                            alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  z <- statistic / sqrt(max(variance, .Machine$double.eps))
  p <- switch(
    alternative,
    two.sided = 2 * pnorm(-abs(z)),
    greater = pnorm(z, lower.tail = FALSE),
    less = pnorm(z)
  )
  list(z = z, p = min(max(p, 0), 1))
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

paired_difference_interval <- function(x, y, conf.level = 0.95) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (!length(x)) {
    return(c(estimate = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  delta <- x - y
  estimate <- mean(delta)
  se <- if (length(delta) > 1L) sd(delta) / sqrt(length(delta)) else NA_real_
  z <- qnorm(1 - (1 - conf.level) / 2)
  c(
    estimate = estimate,
    se = se,
    lower = estimate - z * se,
    upper = estimate + z * se
  )
}

matrix_max_abs <- function(M) {
  if (!length(M)) return(NA_real_)
  max(abs(M))
}

matrix_frobenius <- function(M) {
  sqrt(sum(M^2))
}
