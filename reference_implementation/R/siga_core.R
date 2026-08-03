# Core functions for the SIGA numerical studies.
#
# This file implements the exact equal-weight absolute-imbalance (range)
# Pocock-Simon allocation rule used in the manuscript, reusable one-path and
# three-path allocation calibrations, SIGA-S, SIGA-R, the fixed-score
# randomization test, and the lattice-normal mixture used only for unadjusted
# binary superiority with SIGA-S.
#
# The implementation uses base R only.

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

normalise_seed <- function(seed) {
  seed <- as.double(seed)
  if (length(seed) != 1L || !is.finite(seed)) {
    stop("seed must be one finite number.", call. = FALSE)
  }
  as.integer(abs(seed) %% 2147483646 + 1)
}

seed_value <- function(base_seed, scenario_id, replicate_id, stream_id) {
  modulus <- 2147483646
  value <- as.double(base_seed) +
    1000003 * as.double(stream_id) +
    104729 * as.double(scenario_id) +
    1009 * as.double(replicate_id)
  as.integer(value %% modulus + 1)
}

ensure_directory <- function(path) {
  path <- path.expand(path)
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create directory: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

safe_write_csv <- function(data, path) {
  parent <- ensure_directory(dirname(path))
  temporary <- tempfile(pattern = basename(path), tmpdir = parent, fileext = ".tmp")
  write.csv(data, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) {
    if (!file.copy(temporary, path, overwrite = TRUE)) {
      unlink(temporary)
      stop("Could not write ", path, call. = FALSE)
    }
    unlink(temporary)
  }
  invisible(path)
}

append_csv <- function(data, path) {
  if (!nrow(data)) return(invisible(path))
  ensure_directory(dirname(path))
  exists <- file.exists(path) && file.info(path)$size > 0
  write.table(
    data,
    file = path,
    sep = ",",
    row.names = FALSE,
    col.names = !exists,
    append = exists,
    quote = TRUE,
    qmethod = "double",
    na = ""
  )
  invisible(path)
}

weighted_tabulate <- function(index, weights, nbins) {
  index <- as.integer(index)
  weights <- as.numeric(weights)
  nbins <- as.integer(nbins)
  if (length(index) != length(weights)) {
    stop("index and weights must have the same length.", call. = FALSE)
  }
  out <- numeric(nbins)
  if (!length(index)) return(out)
  summed <- rowsum(weights, group = index, reorder = FALSE)
  out[as.integer(rownames(summed))] <- summed[, 1L]
  out
}

check_binary_factor_matrix <- function(X) {
  X <- as.matrix(X)
  storage.mode(X) <- "integer"
  if (!nrow(X) || !ncol(X) || anyNA(X) || any(!X %in% c(0L, 1L))) {
    stop("All minimization factors must be present and coded 0/1.", call. = FALSE)
  }
  X
}

all_binary_patterns <- function(K) {
  K <- as.integer(K)
  if (length(K) != 1L || K < 1L || K > 20L) {
    stop("K must be an integer between 1 and 20.", call. = FALSE)
  }
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

validate_minimization_inputs <- function(K, pbc, weights) {
  K <- as.integer(K)
  pbc <- as.numeric(pbc)
  weights <- as.numeric(weights)
  if (length(pbc) != 1L || !is.finite(pbc) || pbc <= 0.5 || pbc >= 1) {
    stop("pbc must be strictly between 0.5 and 1.", call. = FALSE)
  }
  if (length(weights) != K + 1L || any(!is.finite(weights)) || any(weights <= 0)) {
    stop("weights must contain K+1 positive finite values.", call. = FALSE)
  }
  invisible(TRUE)
}

profile_prob_independent <- function(factor_prob, patterns = NULL) {
  factor_prob <- as.numeric(factor_prob)
  if (any(!is.finite(factor_prob)) || any(factor_prob <= 0 | factor_prob >= 1)) {
    stop("Each factor probability must be strictly between 0 and 1.", call. = FALSE)
  }
  patterns <- patterns %||% all_binary_patterns(length(factor_prob))
  patterns <- check_binary_factor_matrix(patterns)
  probability <- apply(patterns, 1L, function(x) {
    prod(ifelse(x == 1L, factor_prob, 1 - factor_prob))
  })
  probability / sum(probability)
}

profile_prob_latent_correlated <- function(factor_prob,
                                           latent_strength = 1,
                                           mixing_prob = 0.5,
                                           patterns = NULL) {
  factor_prob <- as.numeric(factor_prob)
  latent_strength <- as.numeric(latent_strength)
  mixing_prob <- as.numeric(mixing_prob)
  if (any(!is.finite(factor_prob)) || any(factor_prob <= 0 | factor_prob >= 1)) {
    stop("Each factor probability must be strictly between 0 and 1.", call. = FALSE)
  }
  if (length(latent_strength) != 1L || !is.finite(latent_strength) || latent_strength < 0) {
    stop("latent_strength must be finite and nonnegative.", call. = FALSE)
  }
  if (length(mixing_prob) != 1L || !is.finite(mixing_prob) ||
      mixing_prob <= 0 || mixing_prob >= 1) {
    stop("mixing_prob must be strictly between 0 and 1.", call. = FALSE)
  }
  K <- length(factor_prob)
  patterns <- patterns %||% all_binary_patterns(K)
  patterns <- check_binary_factor_matrix(patterns)
  intercept <- q0 <- q1 <- numeric(K)
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
  component_probability <- function(q) {
    apply(patterns, 1L, function(x) prod(ifelse(x == 1L, q, 1 - q)))
  }
  probability <- (1 - mixing_prob) * component_probability(q0) +
    mixing_prob * component_probability(q1)
  probability <- probability / sum(probability)
  attr(probability, "latent_intercept") <- intercept
  attr(probability, "conditional_prob_0") <- q0
  attr(probability, "conditional_prob_1") <- q1
  probability
}

profile_marginals <- function(patterns, profile_prob) {
  patterns <- check_binary_factor_matrix(patterns)
  profile_prob <- as.numeric(profile_prob)
  drop(crossprod(profile_prob, patterns))
}

profile_correlations <- function(patterns, profile_prob) {
  patterns <- check_binary_factor_matrix(patterns)
  profile_prob <- as.numeric(profile_prob)
  means <- profile_marginals(patterns, profile_prob)
  second <- crossprod(patterns * sqrt(profile_prob), patterns * sqrt(profile_prob))
  covariance <- second - tcrossprod(means)
  standard_deviation <- sqrt(pmax(diag(covariance), 0))
  denominator <- outer(standard_deviation, standard_deviation)
  correlation <- covariance / denominator
  correlation[!is.finite(correlation)] <- 0
  diag(correlation) <- 1
  correlation
}

nearest_psd <- function(M, tolerance = 1e-12) {
  M <- as.matrix(M)
  M <- (M + t(M)) / 2
  decomposition <- eigen(M, symmetric = TRUE)
  cutoff <- tolerance * max(1, max(abs(decomposition$values)))
  values <- pmax(decomposition$values, 0)
  values[values < cutoff] <- 0
  out <- decomposition$vectors %*% (values * t(decomposition$vectors))
  (out + t(out)) / 2
}

# Generate one treatment path from the historical equal-weight absolute-range
# criterion. The first weight is for overall balance; the remaining K weights
# are for the active levels of the K minimization factors.
absolute_minimization_assign <- function(X,
                                         pbc = 0.80,
                                         weights = NULL,
                                         seed = 1) {
  X <- check_binary_factor_matrix(X)
  n <- nrow(X)
  K <- ncol(X)
  weights <- as.numeric(weights %||% rep(1, K + 1L))
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

# Internal vectorized update used by allocation calibration and RT.
.absolute_minimization_step <- function(overall,
                                        marginal,
                                        X,
                                        pbc,
                                        weights) {
  m <- length(overall)
  K <- ncol(X)
  rows <- seq_len(m)
  score_plus <- weights[1L] * abs(overall + 1L)
  score_minus <- weights[1L] * abs(overall - 1L)
  columns <- matrix(0L, nrow = m, ncol = K)
  for (j in seq_len(K)) {
    columns[, j] <- 2L * (j - 1L) + X[, j] + 1L
    index <- rows + (columns[, j] - 1L) * m
    current <- marginal[index]
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
  overall <- overall + z
  for (j in seq_len(K)) {
    index <- rows + (columns[, j] - 1L) * m
    marginal[index] <- marginal[index] + z
  }
  list(z = z, overall = overall, marginal = marginal)
}

# Fixed-score reference randomization test. All score columns use the same B
# regenerated paths, which is important for equivalence and paired comparisons.
fixed_score_randomization_test <- function(X,
                                           scores,
                                           z_observed,
                                           B = 4999L,
                                           pbc = 0.80,
                                           weights = NULL,
                                           seed = 1,
                                           tolerance = 1e-12) {
  X <- check_binary_factor_matrix(X)
  scores <- as.matrix(scores)
  storage.mode(scores) <- "double"
  z_observed <- as.integer(z_observed)
  n <- nrow(X)
  K <- ncol(X)
  L <- ncol(scores)
  B <- as.integer(B)
  weights <- as.numeric(weights %||% rep(1, K + 1L))
  validate_minimization_inputs(K, pbc, weights)
  if (nrow(scores) != n || length(z_observed) != n || B < 1L) {
    stop("Dimensions or B are invalid in fixed_score_randomization_test().", call. = FALSE)
  }
  set.seed(normalise_seed(seed))
  statistic_observed <- drop(0.5 * crossprod(z_observed, scores))
  statistic_randomized <- matrix(0, nrow = B, ncol = L)
  overall <- integer(B)
  marginal <- matrix(0L, nrow = B, ncol = 2L * K)
  for (i in seq_len(n)) {
    Xi <- matrix(rep(X[i, ], each = B), nrow = B, ncol = K)
    update <- .absolute_minimization_step(overall, marginal, Xi, pbc, weights)
    z <- update$z
    overall <- update$overall
    marginal <- update$marginal
    statistic_randomized <- statistic_randomized + 0.5 * tcrossprod(z, scores[i, ])
  }
  greater <- less <- two_sided <- numeric(L)
  for (j in seq_len(L)) {
    greater[j] <- (1 + sum(statistic_randomized[, j] >=
                              statistic_observed[j] - tolerance)) / (B + 1)
    less[j] <- (1 + sum(statistic_randomized[, j] <=
                           statistic_observed[j] + tolerance)) / (B + 1)
    two_sided[j] <- (1 + sum(abs(statistic_randomized[, j]) >=
                                abs(statistic_observed[j]) - tolerance)) / (B + 1)
  }
  list(
    statistic_observed = statistic_observed,
    greater = greater,
    less = less,
    two_sided = two_sided,
    randomization_mean = colMeans(statistic_randomized),
    randomization_variance = apply(statistic_randomized, 2L, var),
    regenerated_statistics = statistic_randomized,
    B = B
  )
}

# One-path allocation-only calibration for SIGA-S.
calibrate_one_path_design <- function(B0,
                                      n,
                                      patterns,
                                      profile_prob,
                                      pbc = 0.80,
                                      weights = NULL,
                                      seed = 1,
                                      batch_size = 1000L,
                                      progress = interactive()) {
  B0 <- as.integer(B0)
  n <- as.integer(n)
  patterns <- check_binary_factor_matrix(patterns)
  profile_prob <- as.numeric(profile_prob)
  K <- ncol(patterns)
  J <- nrow(patterns)
  weights <- as.numeric(weights %||% rep(1, K + 1L))
  validate_minimization_inputs(K, pbc, weights)
  if (B0 < 2L || n < 2L || J != 2L ^ K || length(profile_prob) != J ||
      any(!is.finite(profile_prob)) || any(profile_prob <= 0)) {
    stop("Invalid one-path calibration inputs.", call. = FALSE)
  }
  profile_prob <- profile_prob / sum(profile_prob)
  batch_size <- max(1L, min(as.integer(batch_size), B0))
  set.seed(normalise_seed(seed))

  sum_u <- numeric(J)
  sum_uu <- matrix(0, J, J)
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
    treated <- integer(m)
    for (i in seq_len(n)) {
      sid <- sample.int(J, size = m, replace = TRUE, prob = profile_prob)
      X <- patterns[sid, , drop = FALSE]
      update <- .absolute_minimization_step(overall, marginal, X, pbc, weights)
      z <- update$z
      overall <- update$overall
      marginal <- update$marginal
      index <- rows + (sid - 1L) * m
      counts[index] <- counts[index] + 1L
      imbalance[index] <- imbalance[index] + z
      treated <- treated + as.integer(z == 1L)
    }
    U <- matrix(0, nrow = m, ncol = J)
    nonzero <- counts > 0L
    U[nonzero] <- imbalance[nonzero] / sqrt(counts[nonzero])
    sum_u <- sum_u + colSums(U)
    sum_uu <- sum_uu + crossprod(U)
    treated_count <- treated_count + tabulate(treated + 1L, nbins = n + 1L)
    completed <- completed + m
    if (progress) message("one-path calibration: ", completed, "/", B0)
  }
  mean_u <- sum_u / B0
  gamma <- (sum_uu - B0 * tcrossprod(mean_u)) / (B0 - 1)
  structure(list(
    gamma = nearest_psd(gamma),
    mean_u = mean_u,
    treated_count_prob = treated_count / B0,
    profile_prob = profile_prob,
    patterns = patterns,
    B0 = B0,
    n = n,
    K = K,
    J = J,
    pbc = pbc,
    weights = weights,
    seed = normalise_seed(seed),
    elapsed_seconds = proc.time()[3L] - start,
    calibration_type = "one_path"
  ), class = c("siga_calibration", "list"))
}

# Three-path calibration. Copy 0 provides Gamma; the pair blocks (0,1) and
# (0,2) provide two exchangeable estimates of Psi, which are averaged.
calibrate_three_path_design <- function(B0,
                                        n,
                                        patterns,
                                        profile_prob,
                                        pbc = 0.80,
                                        weights = NULL,
                                        seed = 1,
                                        batch_size = 1000L,
                                        progress = interactive()) {
  B0 <- as.integer(B0)
  n <- as.integer(n)
  patterns <- check_binary_factor_matrix(patterns)
  profile_prob <- as.numeric(profile_prob)
  K <- ncol(patterns)
  J <- nrow(patterns)
  weights <- as.numeric(weights %||% rep(1, K + 1L))
  validate_minimization_inputs(K, pbc, weights)
  if (B0 < 2L || n < 2L || J != 2L ^ K || length(profile_prob) != J ||
      any(!is.finite(profile_prob)) || any(profile_prob <= 0)) {
    stop("Invalid three-path calibration inputs.", call. = FALSE)
  }
  profile_prob <- profile_prob / sum(profile_prob)
  batch_size <- max(1L, min(as.integer(batch_size), B0))
  set.seed(normalise_seed(seed))

  sum_u <- numeric(J)
  sum_uu <- matrix(0, J, J)
  sum_d <- lapply(seq_len(3L), function(i) numeric(J))
  sum_dd <- lapply(seq_len(3L), function(i) matrix(0, J, J))
  sum_d01 <- sum_d02 <- sum_d12 <- matrix(0, J, J)
  sum_g01 <- sum_g02 <- numeric(J)
  sum_g01g01 <- sum_g02g02 <- sum_g01g02 <- matrix(0, J, J)
  treated_count <- numeric(n + 1L)
  completed <- 0L
  start <- proc.time()[3L]

  while (completed < B0) {
    m <- min(batch_size, B0 - completed)
    rows <- seq_len(m)
    overall <- lapply(seq_len(3L), function(i) integer(m))
    marginal <- lapply(seq_len(3L), function(i) matrix(0L, nrow = m, ncol = 2L * K))
    counts <- matrix(0L, nrow = m, ncol = J)
    imbalance <- lapply(seq_len(3L), function(i) matrix(0L, nrow = m, ncol = J))
    pair01 <- pair02 <- matrix(0L, nrow = m, ncol = J)
    treated <- integer(m)
    for (i in seq_len(n)) {
      sid <- sample.int(J, size = m, replace = TRUE, prob = profile_prob)
      X <- patterns[sid, , drop = FALSE]
      z_copy <- vector("list", 3L)
      for (copy in seq_len(3L)) {
        update <- .absolute_minimization_step(
          overall[[copy]], marginal[[copy]], X, pbc, weights
        )
        z_copy[[copy]] <- update$z
        overall[[copy]] <- update$overall
        marginal[[copy]] <- update$marginal
      }
      index <- rows + (sid - 1L) * m
      counts[index] <- counts[index] + 1L
      for (copy in seq_len(3L)) {
        imbalance[[copy]][index] <- imbalance[[copy]][index] + z_copy[[copy]]
      }
      pair01[index] <- pair01[index] + z_copy[[1L]] * z_copy[[2L]]
      pair02[index] <- pair02[index] + z_copy[[1L]] * z_copy[[3L]]
      treated <- treated + as.integer(z_copy[[1L]] == 1L)
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
    treated_count <- treated_count + tabulate(treated + 1L, nbins = n + 1L)
    completed <- completed + m
    if (progress) message("three-path calibration: ", completed, "/", B0)
  }

  covariance_from_sums <- function(sum_x, sum_xx) {
    mean_x <- sum_x / B0
    (sum_xx - B0 * tcrossprod(mean_x)) / (B0 - 1)
  }
  cross_covariance_from_sums <- function(sum_x, sum_y, sum_xy) {
    (sum_xy - tcrossprod(sum_x, sum_y) / B0) / (B0 - 1)
  }
  mean_u <- sum_u / B0
  gamma <- covariance_from_sums(sum_u, sum_uu)
  sigma_d <- lapply(seq_len(3L), function(i) {
    nearest_psd(covariance_from_sums(sum_d[[i]], sum_dd[[i]]))
  })
  psi01 <- nearest_psd(covariance_from_sums(sum_g01, sum_g01g01))
  psi02 <- nearest_psd(covariance_from_sums(sum_g02, sum_g02g02))
  psi <- nearest_psd((psi01 + psi02) / 2)

  structure(list(
    gamma = nearest_psd(gamma),
    mean_u = mean_u,
    sigma_d = sigma_d,
    psi = psi,
    psi01 = psi01,
    psi02 = psi02,
    mean_g01 = sum_g01 / B0,
    mean_g02 = sum_g02 / B0,
    cross_d01 = cross_covariance_from_sums(sum_d[[1L]], sum_d[[2L]], sum_d01),
    cross_d02 = cross_covariance_from_sums(sum_d[[1L]], sum_d[[3L]], sum_d02),
    cross_d12 = cross_covariance_from_sums(sum_d[[2L]], sum_d[[3L]], sum_d12),
    cross_g0102 = cross_covariance_from_sums(sum_g01, sum_g02, sum_g01g02),
    treated_count_prob = treated_count / B0,
    profile_prob = profile_prob,
    patterns = patterns,
    B0 = B0,
    n = n,
    K = K,
    J = J,
    pbc = pbc,
    weights = weights,
    seed = normalise_seed(seed),
    elapsed_seconds = proc.time()[3L] - start,
    calibration_type = "three_path"
  ), class = c("siga_calibration", "list"))
}

null_space_basis_for_profile_mean <- function(profile_prob) {
  profile_prob <- as.numeric(profile_prob)
  J <- length(profile_prob)
  if (J < 2L || any(!is.finite(profile_prob)) || any(profile_prob <= 0)) {
    stop("profile_prob must contain at least two positive values.", call. = FALSE)
  }
  qfit <- qr(matrix(profile_prob, ncol = 1L), LAPACK = FALSE)
  qr.Q(qfit, complete = TRUE)[, -1L, drop = FALSE]
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
  eig <- eigen((M + t(M)) / 2, symmetric = TRUE)
  ordered <- order(eig$values)
  make_direction <- function(index) {
    y <- eig$vectors[, index]
    x <- drop(Rinv %*% y)
    d <- drop(Q %*% x)
    d <- d - sum(profile_prob * d) / sum(profile_prob)
    d <- d / max(abs(d))
    numerator <- drop(crossprod(d, psi %*% d))
    denominator <- sum(profile_prob * d^2)
    list(direction = d, ratio = numerator / denominator,
         numerator = numerator, denominator = denominator)
  }
  minimum <- make_direction(ordered[1L])
  maximum <- make_direction(ordered[length(ordered)])
  list(
    min_ratio = minimum$ratio,
    min_direction = minimum$direction,
    max_ratio = maximum$ratio,
    max_direction = maximum$direction,
    all_ratios = sort(eig$values),
    basis = Q
  )
}

make_boundary_scores <- function(y, A, X, boundary) {
  y <- as.numeric(y)
  A <- as.numeric(A)
  X <- check_binary_factor_matrix(X)
  if (length(y) != length(A) || length(y) != nrow(X)) {
    stop("y, A and X must contain the same participants.", call. = FALSE)
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
  X <- check_binary_factor_matrix(X)
  id <- joint_stratum_id(X)
  counts <- tabulate(id, nbins = J)
  sums <- weighted_tabulate(id, score, J)
  means <- numeric(J)
  nonempty <- counts > 0L
  means[nonempty] <- sums[nonempty] / counts[nonempty]
  residual <- score - means[id]
  list(id = id, counts = counts, means = means, residual = residual)
}

siga_sampling_variance <- function(score, X, calibration, truncate_kappa = TRUE) {
  if (!inherits(calibration, "siga_calibration")) {
    stop("calibration must be produced by a SIGA calibration function.", call. = FALSE)
  }
  X <- check_binary_factor_matrix(X)
  if (nrow(X) != calibration$n || ncol(X) != calibration$K) {
    stop("Trial dimensions do not match calibration dimensions.", call. = FALSE)
  }
  decomposition <- score_decomposition(score, X, calibration$J)
  sqrt_counts <- sqrt(decomposition$counts)
  Omega <- tcrossprod(sqrt_counts) * calibration$gamma
  nonempty <- decomposition$counts > 0L
  Jplus <- sum(nonempty)
  trace_low <- sum(diag(Omega)[nonempty] / decomposition$counts[nonempty])
  denominator <- calibration$n - Jplus
  kappa_raw <- if (denominator > 0L) {
    (calibration$n - trace_low) / denominator
  } else {
    0
  }
  kappa <- if (truncate_kappa) max(kappa_raw, 0) else kappa_raw
  low <- drop(crossprod(decomposition$means, Omega %*% decomposition$means))
  within <- sum(decomposition$residual^2)
  variance <- 0.25 * (low + kappa * within)
  if (!is.finite(variance) || variance <= 0) variance <- .Machine$double.eps
  list(
    variance = variance,
    kappa = kappa,
    kappa_raw = kappa_raw,
    kappa_truncated = isTRUE(truncate_kappa) && kappa_raw < 0,
    low_component = 0.25 * low,
    within_component = 0.25 * kappa * within,
    decomposition = decomposition
  )
}

siga_randomization_variance <- function(sampling_variance_result,
                                         d,
                                         calibration,
                                         epsilon = NULL) {
  if (is.null(calibration$psi)) {
    stop("SIGA-R requires a three-path calibration containing psi.", call. = FALSE)
  }
  d <- as.numeric(d)
  if (length(d) != calibration$J) {
    stop("d must have one entry per joint profile.", call. = FALSE)
  }
  epsilon <- as.numeric(epsilon %||% (1 / calibration$n))
  if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0) {
    stop("epsilon must be positive and finite.", call. = FALSE)
  }
  counts <- sampling_variance_result$decomposition$counts
  pi_hat <- counts / sum(counts)
  correction_matrix <- calibration$psi - diag(pi_hat, nrow = length(pi_hat))
  scalar_gap <- drop(crossprod(d, correction_matrix %*% d))
  correction <- calibration$n * scalar_gap / 16
  variance_raw <- sampling_variance_result$variance + correction
  lower_bound <- epsilon * sampling_variance_result$variance
  variance <- if (is.finite(variance_raw)) max(variance_raw, lower_bound) else lower_bound
  list(
    variance = variance,
    variance_raw = variance_raw,
    correction = correction,
    scalar_gap = scalar_gap,
    pi_hat = pi_hat,
    lower_bound = lower_bound,
    safeguard_active = !is.finite(variance_raw) || variance_raw < lower_bound,
    variance_ratio = variance_raw / sampling_variance_result$variance
  )
}

gaussian_pvalue <- function(statistic,
                            variance,
                            alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  statistic <- as.numeric(statistic)
  variance <- as.numeric(variance)
  if (length(statistic) != 1L || length(variance) != 1L ||
      !is.finite(statistic) || !is.finite(variance) || variance <= 0) {
    stop("statistic and variance must be finite scalars with variance > 0.", call. = FALSE)
  }
  z <- statistic / sqrt(variance)
  p <- switch(
    alternative,
    two.sided = 2 * pnorm(-abs(z)),
    greater = pnorm(z, lower.tail = FALSE),
    less = pnorm(z)
  )
  list(z = z, p = min(max(p, 0), 1))
}

# Lattice-normal mixture used only for unadjusted binary superiority with
# SIGA-S. It is not used for SIGA-R in the final manuscript.
lattice_normal_mixture_pvalue <- function(y,
                                           A,
                                           statistic,
                                           variance,
                                           treated_count_prob,
                                           alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  y <- as.integer(y)
  A <- as.integer(A)
  n <- length(y)
  if (length(A) != n || any(!y %in% c(0L, 1L)) || any(!A %in% c(0L, 1L))) {
    stop("y and A must be equal-length binary vectors.", call. = FALSE)
  }
  if (length(treated_count_prob) != n + 1L) {
    stop("treated_count_prob must have n+1 entries.", call. = FALSE)
  }
  weights <- as.numeric(treated_count_prob)
  weights <- weights / sum(weights)
  y_total <- sum(y)
  tail_probability <- 0
  for (k in 0:n) {
    wk <- weights[k + 1L]
    if (!is.finite(wk) || wk <= 0) next
    lower_x <- max(0L, k - (n - y_total))
    upper_x <- min(k, y_total)
    x <- seq.int(lower_x, upper_x)
    t_lattice <- x - k * y_total / n
    cell <- pnorm((t_lattice + 0.5) / sqrt(variance)) -
      pnorm((t_lattice - 0.5) / sqrt(variance))
    normalizer <- sum(cell)
    if (!is.finite(normalizer) || normalizer <= 0) next
    keep <- switch(
      alternative,
      two.sided = abs(t_lattice) >= abs(statistic) - 1e-12,
      greater = t_lattice >= statistic - 1e-12,
      less = t_lattice <= statistic + 1e-12
    )
    tail_probability <- tail_probability + wk * sum(cell[keep]) / normalizer
  }
  min(max(tail_probability, 0), 1)
}

wilson_interval <- function(x, n, conf_level = 0.95) {
  x <- as.numeric(x)
  n <- as.numeric(n)
  if (!is.finite(n) || n <= 0) {
    return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  z <- qnorm(1 - (1 - conf_level) / 2)
  phat <- x / n
  denominator <- 1 + z^2 / n
  center <- (phat + z^2 / (2 * n)) / denominator
  half <- z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2)) / denominator
  c(estimate = phat, lower = max(0, center - half), upper = min(1, center + half))
}

paired_difference_interval <- function(x, y, conf_level = 0.95) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y)
  delta <- x[keep] - y[keep]
  if (!length(delta)) {
    return(c(estimate = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  estimate <- mean(delta)
  se <- if (length(delta) > 1L) sd(delta) / sqrt(length(delta)) else NA_real_
  z <- qnorm(1 - (1 - conf_level) / 2)
  c(estimate = estimate, se = se,
    lower = estimate - z * se, upper = estimate + z * se)
}

matrix_max_abs <- function(M) if (length(M)) max(abs(M)) else NA_real_
matrix_frobenius <- function(M) sqrt(sum(M^2))
