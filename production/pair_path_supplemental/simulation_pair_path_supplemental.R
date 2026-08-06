#!/usr/bin/env Rscript

# Standalone full program for the revised supplemental SIGA pair-path sensitivity analysis.
# Semantically reconstructed from the successful R parse output of the 105077-byte
# production source used for the completed multi-machine run.
# Original source SHA-256:
# 35e0adb0da0b7b584fb093cee5f4869c2ad7e79242cd7b02c0c0ea2e19988ca6
#
# The only post-run change is unname() in rejection_summary(), which corrects
# aggregate column names and does not alter trial simulation or RT calculations.

options(stringsAsFactors = FALSE, warn = 1)

`%||%` <- function(x, 
    y) if (is.null(x)) y else x

normalise_seed <- function(seed) {
    seed <- as.double(seed)
    if (!is.finite(seed)) 
        stop("seed must be finite.", call. = FALSE)
    as.integer(abs(seed)%%2147483646 + 1)
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
    if (K < 1L || K > 20L) 
        stop("K must be between 1 and 20.", call. = FALSE)
    out <- as.matrix(expand.grid(rep(list(c(0L, 1L)), K), KEEP.OUT.ATTRS = FALSE, 
        stringsAsFactors = FALSE))
    storage.mode(out) <- "integer"
    colnames(out) <- paste0("X", seq_len(K))
    out
}

joint_stratum_id <- function(X) {
    X <- check_binary_factor_matrix(X)
    K <- ncol(X)
    as.integer(1L + drop(X %*% (2L^(0:(K - 1L)))))
}

nearest_psd <- function(M, tolerance = 1e-12) {
    M <- as.matrix(M)
    M <- (M + t(M))/2
    ee <- eigen(M, symmetric = TRUE)
    cutoff <- tolerance * max(1, max(abs(ee$values)))
    values <- pmax(ee$values, 0)
    values[values < cutoff] <- 0
    out <- ee$vectors %*% (values * t(ee$vectors))
    (out + t(out))/2
}

profile_prob_independent <- function(factor_prob, patterns = NULL) {
    factor_prob <- as.numeric(factor_prob)
    if (any(!is.finite(factor_prob)) || any(factor_prob <= 0 | 
        factor_prob >= 1)) {
        stop("Each factor probability must be strictly between 0 and 1.", 
            call. = FALSE)
    }
    patterns <- patterns %||% all_binary_patterns(length(factor_prob))
    prob <- apply(patterns, 1L, function(x) {
        prod(ifelse(x == 1L, factor_prob, 1 - factor_prob))
    })
    prob/sum(prob)
}

profile_prob_latent_correlated <- function(factor_prob, latent_strength = 1, 
    mixing_prob = 0.5, patterns = NULL) {
    factor_prob <- as.numeric(factor_prob)
    if (any(!is.finite(factor_prob)) || any(factor_prob <= 0 | 
        factor_prob >= 1)) {
        stop("Each factor probability must be strictly between 0 and 1.", 
            call. = FALSE)
    }
    if (!is.finite(latent_strength) || latent_strength < 0) {
        stop("latent_strength must be finite and nonnegative.", 
            call. = FALSE)
    }
    if (!is.finite(mixing_prob) || mixing_prob <= 0 || mixing_prob >= 
        1) {
        stop("mixing_prob must be strictly between 0 and 1.", 
            call. = FALSE)
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
        intercept[j] <- uniroot(objective, interval = c(-30, 
            30), tol = 1e-13)$root
        q0[j] <- plogis(intercept[j] - latent_strength)
        q1[j] <- plogis(intercept[j] + latent_strength)
    }
    component_prob <- function(q) {
        apply(patterns, 1L, function(x) prod(ifelse(x == 1L, 
            q, 1 - q)))
    }
    prob <- (1 - mixing_prob) * component_prob(q0) + mixing_prob * 
        component_prob(q1)
    prob <- prob/sum(prob)
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
    second <- crossprod(patterns * sqrt(profile_prob), patterns * 
        sqrt(profile_prob))
    covariance <- second - tcrossprod(mu)
    sdv <- sqrt(pmax(diag(covariance), 0))
    denominator <- outer(sdv, sdv)
    correlation <- covariance/denominator
    correlation[!is.finite(correlation)] <- 0
    diag(correlation) <- 1
    correlation
}

validate_minimization_inputs <- function(K, pbc, weights) {
    if (!(pbc >= 0.5 && pbc <= 1)) {
        stop("pbc must lie in [0.5, 1].", call. = FALSE)
    }
    if (length(weights) != K + 1L) {
        stop("weights must have length K + 1: overall plus one per factor.", 
            call. = FALSE)
    }
    if (any(!is.finite(weights)) || any(weights < 0)) {
        stop("weights must be finite and nonnegative.", call. = FALSE)
    }
    invisible(TRUE)
}

ps_assign_R <- function(X, pbc = 0.8, weights_ = NULL, seed = 1) {
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
            score_plus <- score_plus + weights[j + 1L] * abs(current + 
                1L)
            score_minus <- score_minus + weights[j + 1L] * abs(current - 
                1L)
        }
        u <- runif(1L)
        zi <- if (score_plus < score_minus) {
            if (u < pbc) 
                1L
            else -1L
        }
        else if (score_plus > score_minus) {
            if (u < 1 - pbc) 
                1L
            else -1L
        }
        else {
            if (u < 0.5) 
                1L
            else -1L
        }
        z[i] <- zi
        overall <- overall + zi
        for (j in seq_len(K)) {
            marginal[j, X[i, j] + 1L] <- marginal[j, X[i, j] + 
                1L] + zi
        }
    }
    z
}

rt_pvalues_R <- function(X, scores, z_obs, B = 1999L, pbc = 0.8, 
    weights_ = NULL, seed = 1, tolerance = 1e-12, return_randomization_moments = TRUE, 
    return_shape_diagnostics = TRUE) {
    X <- check_binary_factor_matrix(X)
    scores <- as.matrix(scores)
    storage.mode(scores) <- "double"
    z_obs <- as.integer(z_obs)
    n <- nrow(X)
    K <- ncol(X)
    L <- ncol(scores)
    B <- as.integer(B)
    if (nrow(scores) != n) 
        stop("scores must have n rows.", call. = FALSE)
    if (length(z_obs) != n) 
        stop("z_obs must have length n.", call. = FALSE)
    if (B < 1L) 
        stop("B must be positive.", call. = FALSE)
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
            score_plus <- score_plus + weights[j + 1L] * abs(current + 
                1L)
            score_minus <- score_minus + weights[j + 1L] * abs(current - 
                1L)
        }
        u <- runif(B)
        z <- integer(B)
        better_plus <- score_plus < score_minus
        better_minus <- score_plus > score_minus
        ties <- !(better_plus | better_minus)
        z[better_plus] <- ifelse(u[better_plus] < pbc, 1L, -1L)
        z[better_minus] <- ifelse(u[better_minus] < 1 - pbc, 
            1L, -1L)
        z[ties] <- ifelse(u[ties] < 0.5, 1L, -1L)
        overall <- overall + z
        for (j in seq_len(K)) marginal[, cols[j]] <- marginal[, 
            cols[j]] + z
        for (ell in seq_len(L)) {
            t_rand[, ell] <- t_rand[, ell] + 0.5 * z * scores[i, 
                ell]
        }
    }
    abs_obs <- matrix(abs(t_obs), nrow = B, ncol = L, byrow = TRUE)
    obs_mat <- matrix(t_obs, nrow = B, ncol = L, byrow = TRUE)
    p_two <- (1 + colSums(abs(t_rand) + tolerance >= abs_obs))/(B + 
        1)
    p_upper <- (1 + colSums(t_rand + tolerance >= obs_mat))/(B + 
        1)
    p_lower <- (1 + colSums(t_rand - tolerance <= obs_mat))/(B + 
        1)
    out <- list(statistic = t_obs, two_sided = p_two, greater = p_upper, 
        less = p_lower, B = B, engine = "pure_R_vectorized")
    if (isTRUE(return_randomization_moments)) {
        randomization_mean <- colMeans(t_rand)
        randomization_variance <- if (B > 1L) {
            apply(t_rand, 2L, var)
        }
        else {
            rep(NA_real_, L)
        }
        out$randomization_mean <- randomization_mean
        out$randomization_variance <- randomization_variance
        if (isTRUE(return_shape_diagnostics)) {
            centered <- sweep(t_rand, 2L, randomization_mean, 
                "-")
            second_moment <- colMeans(centered^2)
            valid <- is.finite(second_moment) & second_moment > 
                0
            skewness <- rep(NA_real_, L)
            excess_kurtosis <- rep(NA_real_, L)
            if (any(valid)) {
                skewness[valid] <- colMeans(centered[, valid, 
                  drop = FALSE]^3)/second_moment[valid]^(3/2)
                excess_kurtosis[valid] <- colMeans(centered[, 
                  valid, drop = FALSE]^4)/second_moment[valid]^2 - 
                  3
            }
            out$randomization_skewness <- skewness
            out$randomization_excess_kurtosis <- excess_kurtosis
        }
    }
    out
}

covariance_from_sums <- function(sum_x, sum_xx, n) {
    if (n < 2L) 
        stop("At least two observations are required.", call. = FALSE)
    (sum_xx - tcrossprod(sum_x)/n)/(n - 1L)
}

cross_covariance_from_sums <- function(sum_x, sum_y, sum_xy, 
    n) {
    if (n < 2L) 
        stop("At least two observations are required.", call. = FALSE)
    (sum_xy - tcrossprod(sum_x, sum_y)/n)/(n - 1L)
}

calibrate_pair_path_design_R <- function(B0, n, patterns, 
    profile_prob, pbc = 0.8, weights_ = NULL, seed = 1, batch_size = 1000L, 
    progress = TRUE) {
    B0 <- as.integer(B0)
    n <- as.integer(n)
    patterns <- check_binary_factor_matrix(patterns)
    profile_prob <- as.numeric(profile_prob)
    K <- ncol(patterns)
    J <- nrow(patterns)
    weights <- as.numeric(weights_ %||% rep(1, K + 1L))
    if (B0 < 2L) 
        stop("B0 must be at least 2.", call. = FALSE)
    if (n < 2L) 
        stop("n must be at least 2.", call. = FALSE)
    if (J != 2L^K) 
        stop("patterns must contain all 2^K binary profiles.", 
            call. = FALSE)
    if (length(profile_prob) != J || any(!is.finite(profile_prob)) || 
        any(profile_prob <= 0)) {
        stop("profile_prob must contain positive finite probabilities for all profiles.", 
            call. = FALSE)
    }
    profile_prob <- profile_prob/sum(profile_prob)
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
                score_plus <- weights[1L] * abs(overall[[copy]] + 
                  1L)
                score_minus <- weights[1L] * abs(overall[[copy]] - 
                  1L)
                for (j in seq_len(K)) {
                  idx_j <- rows + (cols[, j] - 1L) * m
                  current <- marginal[[copy]][idx_j]
                  score_plus <- score_plus + weights[j + 1L] * 
                    abs(current + 1L)
                  score_minus <- score_minus + weights[j + 1L] * 
                    abs(current - 1L)
                }
                u <- runif(m)
                z <- integer(m)
                better_plus <- score_plus < score_minus
                better_minus <- score_plus > score_minus
                ties <- !(better_plus | better_minus)
                z[better_plus] <- ifelse(u[better_plus] < pbc, 
                  1L, -1L)
                z[better_minus] <- ifelse(u[better_minus] < 1 - 
                  pbc, 1L, -1L)
                z[ties] <- ifelse(u[ties] < 0.5, 1L, -1L)
                z_copy[[copy]] <- z
                overall[[copy]] <- overall[[copy]] + z
                for (j in seq_len(K)) {
                  idx_j <- rows + (cols[, j] - 1L) * m
                  marginal[[copy]][idx_j] <- marginal[[copy]][idx_j] + 
                    z
                }
            }
            idx_s <- rows + (sid - 1L) * m
            counts[idx_s] <- counts[idx_s] + 1L
            for (copy in seq_len(3L)) {
                imbalance[[copy]][idx_s] <- imbalance[[copy]][idx_s] + 
                  z_copy[[copy]]
            }
            pair01[idx_s] <- pair01[idx_s] + z_copy[[1L]] * z_copy[[2L]]
            pair02[idx_s] <- pair02[idx_s] + z_copy[[1L]] * z_copy[[3L]]
        }
        U <- matrix(0, nrow = m, ncol = J)
        nonzero <- counts > 0L
        U[nonzero] <- imbalance[[1L]][nonzero]/sqrt(counts[nonzero])
        D <- lapply(imbalance, function(x) x/sqrt(n))
        G01 <- pair01/sqrt(n)
        G02 <- pair02/sqrt(n)
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
        treated <- as.integer((n + overall[[1L]])/2L)
        treated_count <- treated_count + tabulate(treated + 1L, 
            nbins = n + 1L)
        completed <- completed + m
        if (isTRUE(progress)) {
            message("  three-copy calibration: ", completed, 
                "/", B0, " paths; elapsed ", sprintf("%.1f min", 
                  (proc.time()[3L] - start)/60))
        }
    }
    gamma <- covariance_from_sums(sum_u, sum_uu, B0)
    sigma_d <- lapply(seq_len(3L), function(copy) {
        covariance_from_sums(sum_d[[copy]], sum_dd[[copy]], B0)
    })
    psi01 <- covariance_from_sums(sum_g01, sum_g01g01, B0)
    psi02 <- covariance_from_sums(sum_g02, sum_g02g02, B0)
    psi <- nearest_psd((psi01 + psi02)/2)
    cross_d01 <- cross_covariance_from_sums(sum_d[[1L]], sum_d[[2L]], 
        sum_d01, B0)
    cross_d02 <- cross_covariance_from_sums(sum_d[[1L]], sum_d[[3L]], 
        sum_d02, B0)
    cross_d12 <- cross_covariance_from_sums(sum_d[[2L]], sum_d[[3L]], 
        sum_d12, B0)
    cross_g0102 <- cross_covariance_from_sums(sum_g01, sum_g02, 
        sum_g01g02, B0)
    structure(list(gamma = nearest_psd(gamma), mean_u = sum_u/B0, 
        sigma_d = lapply(sigma_d, nearest_psd), psi = psi, psi01 = nearest_psd(psi01), 
        psi02 = nearest_psd(psi02), mean_g01 = sum_g01/B0, mean_g02 = sum_g02/B0, 
        cross_d01 = cross_d01, cross_d02 = cross_d02, cross_d12 = cross_d12, 
        cross_g0102 = cross_g0102, treated_count_prob = treated_count/B0, 
        profile_prob = profile_prob, patterns = patterns, B0 = B0, 
        n = n, K = K, J = J, pbc = pbc, weights = weights, seed = seed, 
        elapsed_seconds = proc.time()[3L] - start, engine = "pure_R_three_copy_vectorized"), 
        class = c("siga_pair_calibration", "list"))
}

null_space_basis_for_profile_mean <- function(profile_prob) {
    profile_prob <- as.numeric(profile_prob)
    J <- length(profile_prob)
    if (J < 2L || any(!is.finite(profile_prob)) || any(profile_prob <= 
        0)) {
        stop("profile_prob must contain at least two positive values.", 
            call. = FALSE)
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
        stop("psi and profile_prob dimensions do not agree.", 
            call. = FALSE)
    }
    Q <- null_space_basis_for_profile_mean(profile_prob)
    Pi <- diag(profile_prob, nrow = J, ncol = J)
    A <- crossprod(Q, psi %*% Q)
    C <- crossprod(Q, Pi %*% Q)
    R <- chol((C + t(C))/2)
    Rinv <- backsolve(R, diag(ncol(R)))
    M <- crossprod(Rinv, A %*% Rinv)
    ee <- eigen((M + t(M))/2, symmetric = TRUE)
    order_index <- order(ee$values)
    make_direction <- function(index) {
        y <- ee$vectors[, index]
        x <- drop(Rinv %*% y)
        d <- drop(Q %*% x)
        d <- d/max(abs(d))
        d <- d - sum(profile_prob * d)/sum(profile_prob)
        d <- d/max(abs(d))
        numerator <- drop(crossprod(d, psi %*% d))
        denominator <- sum(profile_prob * d^2)
        list(direction = d, ratio = numerator/denominator, numerator = numerator, 
            denominator = denominator)
    }
    min_result <- make_direction(order_index[1L])
    max_result <- make_direction(order_index[length(order_index)])
    list(min_ratio = min_result$ratio, min_direction = min_result$direction, 
        max_ratio = max_result$ratio, max_direction = max_result$direction, 
        all_ratios = sort(ee$values), basis = Q)
}

scale_direction <- function(direction, target_max_abs = 1) {
    direction <- as.numeric(direction)
    if (!length(direction) || any(!is.finite(direction))) {
        stop("direction must be a finite nonempty vector.", call. = FALSE)
    }
    m <- max(abs(direction))
    if (m <= 0) 
        return(rep(0, length(direction)))
    direction * (target_max_abs/m)
}

scale_direction_for_binary_delta <- function(direction, boundary, 
    target_max_abs_d = 0.8, max_abs_delta = 0.9, safety = 0.98) {
    direction <- scale_direction(direction, 1)
    upper <- max_abs_delta
    lower <- -max_abs_delta
    candidate <- Inf
    positive <- direction > 0
    negative <- direction < 0
    if (any(positive)) {
        candidate <- min(candidate, min((upper - boundary)/direction[positive]))
    }
    if (any(negative)) {
        candidate <- min(candidate, min((lower - boundary)/direction[negative]))
    }
    if (!is.finite(candidate) || candidate <= 0) {
        stop("No positive heterogeneity scale satisfies the binary bounds.", 
            call. = FALSE)
    }
    scale <- min(target_max_abs_d, safety * candidate)
    d <- scale * direction
    delta <- boundary + d
    if (any(abs(delta) >= 1)) {
        stop("Scaled binary risk differences are outside (-1,1).", 
            call. = FALSE)
    }
    d
}

weighted_center <- function(x, weight) {
    x <- as.numeric(x)
    weight <- as.numeric(weight)
    x - sum(weight * x)/sum(weight)
}

pattern_baseline_mean <- function(patterns, profile_prob, 
    beta_master = c(0.5, 0.4, 0.3, 0.2, 0.1), interaction = 0.25) {
    patterns <- check_binary_factor_matrix(patterns)
    K <- ncol(patterns)
    beta <- as.numeric(beta_master)[seq_len(K)]
    mu <- drop(patterns %*% beta)
    if (K >= 2L) 
        mu <- mu + interaction * patterns[, 1L] * patterns[, 
            2L]
    weighted_center(mu, profile_prob)
}

calibrate_logistic_baseline <- function(patterns, profile_prob, 
    target_control_risk = 0.6, beta_master = c(0.35, -0.25, 0.2, 
        -0.15, 0.1), interaction = 0.2) {
    patterns <- check_binary_factor_matrix(patterns)
    profile_prob <- as.numeric(profile_prob)
    K <- ncol(patterns)
    beta <- as.numeric(beta_master)[seq_len(K)]
    eta_no_intercept <- drop(patterns %*% beta)
    if (K >= 2L) {
        eta_no_intercept <- eta_no_intercept + interaction * 
            patterns[, 1L] * patterns[, 2L]
    }
    objective <- function(alpha0) {
        sum(profile_prob * plogis(alpha0 + eta_no_intercept)) - 
            target_control_risk
    }
    alpha0 <- uniroot(objective, interval = c(-30, 30), tol = 1e-13)$root
    list(intercept = alpha0, eta_no_intercept = eta_no_intercept, 
        p0 = plogis(alpha0 + eta_no_intercept), achieved_control_risk = sum(profile_prob * 
            plogis(alpha0 + eta_no_intercept)))
}

calibrate_common_log_odds_model <- function(patterns, profile_prob, 
    target_control_risk, target_risk_difference, beta_master = c(0.35, 
        -0.25, 0.2, -0.15, 0.1), interaction = 0.2) {
    baseline <- calibrate_logistic_baseline(patterns = patterns, 
        profile_prob = profile_prob, target_control_risk = target_control_risk, 
        beta_master = beta_master, interaction = interaction)
    objective <- function(theta) {
        p1 <- plogis(baseline$intercept + baseline$eta_no_intercept + 
            theta)
        sum(profile_prob * (p1 - baseline$p0)) - target_risk_difference
    }
    theta <- uniroot(objective, interval = c(-30, 30), tol = 1e-13)$root
    p1 <- plogis(baseline$intercept + baseline$eta_no_intercept + 
        theta)
    list(type = "binary_common_log_odds", p0 = baseline$p0, p1 = p1, 
        theta = theta, d = p1 - baseline$p0 - target_risk_difference, 
        achieved_control_risk = sum(profile_prob * baseline$p0), 
        achieved_treatment_risk = sum(profile_prob * p1), achieved_effect = sum(profile_prob * 
            (p1 - baseline$p0)))
}

make_continuous_homogeneous_model <- function(patterns, profile_prob, 
    boundary, outcome_sd = 1) {
    mu0 <- pattern_baseline_mean(patterns, profile_prob)
    list(type = "continuous_realistic", boundary = boundary, 
        d = rep(0, nrow(patterns)), delta = rep(boundary, nrow(patterns)), 
        mu0 = mu0, outcome_sd = outcome_sd, individual_effect_sd = 0, 
        achieved_effect = boundary)
}

make_continuous_pure_model <- function(patterns, profile_prob, 
    boundary, direction, target_max_abs_d = 1, shared_noise_sd = 0.05) {
    d <- scale_direction(direction, target_max_abs_d)
    d <- weighted_center(d, profile_prob)
    list(type = "continuous_symmetric", boundary = boundary, 
        d = d, delta = boundary + d, shared_noise_sd = shared_noise_sd, 
        achieved_effect = sum(profile_prob * (boundary + d)))
}

make_continuous_realistic_model <- function(patterns, profile_prob, 
    boundary, direction, target_max_abs_d = 0.5, outcome_sd = 1, 
    individual_effect_sd = 0.25) {
    d <- scale_direction(direction, target_max_abs_d)
    d <- weighted_center(d, profile_prob)
    mu0 <- pattern_baseline_mean(patterns, profile_prob)
    list(type = "continuous_realistic", boundary = boundary, 
        d = d, delta = boundary + d, mu0 = mu0, outcome_sd = outcome_sd, 
        individual_effect_sd = individual_effect_sd, achieved_effect = sum(profile_prob * 
            (boundary + d)))
}

make_binary_homogeneous_model <- function(patterns, profile_prob, 
    boundary, target_control_risk = 0.6) {
    baseline <- calibrate_logistic_baseline(patterns = patterns, 
        profile_prob = profile_prob, target_control_risk = target_control_risk)
    p1 <- baseline$p0 + boundary
    if (any(p1 <= 0 | p1 >= 1)) {
        stop("The constant risk-difference model produces invalid treatment risks.", 
            call. = FALSE)
    }
    list(type = "binary_direct_probability", boundary = boundary, 
        p0 = baseline$p0, p1 = p1, d = rep(0, nrow(patterns)), 
        achieved_control_risk = sum(profile_prob * baseline$p0), 
        achieved_treatment_risk = sum(profile_prob * p1), achieved_effect = sum(profile_prob * 
            (p1 - baseline$p0)))
}

make_binary_symmetric_model <- function(patterns, profile_prob, 
    boundary, direction, target_max_abs_d = 0.8, max_abs_delta = 0.9) {
    d <- scale_direction_for_binary_delta(direction = direction, 
        boundary = boundary, target_max_abs_d = target_max_abs_d, 
        max_abs_delta = max_abs_delta)
    d <- weighted_center(d, profile_prob)
    delta <- boundary + d
    p0 <- 0.5 - 0.5 * delta
    p1 <- 0.5 + 0.5 * delta
    if (any(p0 <= 0 | p0 >= 1 | p1 <= 0 | p1 >= 1)) {
        stop("The symmetric binary stress model produces invalid probabilities.", 
            call. = FALSE)
    }
    list(type = "binary_direct_probability", boundary = boundary, 
        p0 = p0, p1 = p1, d = p1 - p0 - boundary, achieved_control_risk = sum(profile_prob * 
            p0), achieved_treatment_risk = sum(profile_prob * 
            p1), achieved_effect = sum(profile_prob * (p1 - p0)))
}

generate_profile_sequence <- function(n, patterns, profile_prob) {
    sid <- sample.int(nrow(patterns), size = n, replace = TRUE, 
        prob = profile_prob)
    list(id = sid, X = patterns[sid, , drop = FALSE])
}

generate_outcome_from_model <- function(profile_id, A, z, 
    model) {
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
        }
        else {
            rep(0, n)
        }
        return(y0 + A * (model$delta[profile_id] + eta))
    }
    if (model$type == "binary_direct_probability" || model$type == 
        "binary_common_log_odds") {
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
        stop("y, A and X must contain the same number of participants.", 
            call. = FALSE)
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
    means[nonempty] <- sums[nonempty]/counts[nonempty]
    residual <- score - means[id]
    list(id = id, counts = counts, means = means, residual = residual)
}

siga_sampling_variance <- function(score, X, calibration, 
    truncate_kappa = TRUE) {
    if (!inherits(calibration, "siga_pair_calibration")) {
        stop("calibration must be produced by calibrate_pair_path_design_R().", 
            call. = FALSE)
    }
    X <- check_binary_factor_matrix(X)
    if (nrow(X) != calibration$n || ncol(X) != calibration$K) {
        stop("Trial dimensions do not match the allocation calibration.", 
            call. = FALSE)
    }
    dec <- score_decomposition(score, X, calibration$J)
    sqrt_counts <- sqrt(dec$counts)
    Omega <- tcrossprod(sqrt_counts) * calibration$gamma
    nonempty <- dec$counts > 0L
    Jplus <- sum(nonempty)
    trace_low <- sum(diag(Omega)[nonempty]/dec$counts[nonempty])
    denominator <- calibration$n - Jplus
    kappa_raw <- if (denominator > 0L) {
        (calibration$n - trace_low)/denominator
    }
    else {
        0
    }
    kappa <- if (truncate_kappa) 
        max(kappa_raw, 0)
    else kappa_raw
    low <- drop(crossprod(dec$means, Omega %*% dec$means))
    within <- sum(dec$residual^2)
    variance <- 0.25 * (low + kappa * within)
    if (!is.finite(variance) || variance <= 0) 
        variance <- .Machine$double.eps
    list(variance = variance, kappa = kappa, kappa_raw = kappa_raw, 
        low_component = 0.25 * low, within_component = 0.25 * 
            kappa * within, decomposition = dec)
}

siga_randomization_variance <- function(sampling_variance_result, 
    d, calibration, epsilon_n = NULL) {
    d <- as.numeric(d)
    if (length(d) != calibration$J) {
        stop("d must have one entry per joint profile.", call. = FALSE)
    }
    if (is.null(epsilon_n)) 
        epsilon_n <- 1/calibration$n
    epsilon_n <- as.numeric(epsilon_n)
    if (length(epsilon_n) != 1L || !is.finite(epsilon_n) || epsilon_n <= 
        0) {
        stop("epsilon_n must be a positive finite scalar.", call. = FALSE)
    }
    counts <- sampling_variance_result$decomposition$counts
    pi_hat <- counts/sum(counts)
    correction_matrix <- calibration$psi - diag(pi_hat, nrow = length(pi_hat))
    scalar_gap <- drop(crossprod(d, correction_matrix %*% d))
    correction <- calibration$n * scalar_gap/16
    variance_raw <- sampling_variance_result$variance + correction
    variance_floor <- epsilon_n * sampling_variance_result$variance
    safeguard_active <- !is.finite(variance_raw) || variance_raw < 
        variance_floor
    variance <- if (safeguard_active) 
        variance_floor
    else variance_raw
    list(variance = variance, variance_raw = variance_raw, variance_floor = variance_floor, 
        correction = correction, scalar_gap = scalar_gap, pi_hat = pi_hat, 
        epsilon_n = epsilon_n, safeguard_active = safeguard_active, 
        truncated = safeguard_active)
}

gaussian_pvalue <- function(statistic, variance, alternative = c("two.sided", 
    "greater", "less")) {
    alternative <- match.arg(alternative)
    z <- statistic/sqrt(max(variance, .Machine$double.eps))
    p <- switch(alternative, two.sided = 2 * pnorm(-abs(z)), 
        greater = pnorm(z, lower.tail = FALSE), less = pnorm(z))
    list(z = z, p = min(max(p, 0), 1))
}

wilson_interval <- function(x, n, conf.level = 0.95) {
    if (!is.finite(n) || n <= 0) {
        return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
    }
    z <- qnorm(1 - (1 - conf.level)/2)
    phat <- x/n
    denominator <- 1 + z^2/n
    center <- (phat + z^2/(2 * n))/denominator
    half <- z * sqrt(phat * (1 - phat)/n + z^2/(4 * n^2))/denominator
    c(estimate = phat, lower = max(0, center - half), upper = min(1, 
        center + half))
}

paired_difference_interval <- function(x, y, conf.level = 0.95) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]
    y <- y[keep]
    if (!length(x)) {
        return(c(estimate = NA_real_, se = NA_real_, lower = NA_real_, 
            upper = NA_real_))
    }
    delta <- x - y
    estimate <- mean(delta)
    se <- if (length(delta) > 1L) 
        sd(delta)/sqrt(length(delta))
    else NA_real_
    z <- qnorm(1 - (1 - conf.level)/2)
    c(estimate = estimate, se = se, lower = estimate - z * se, 
        upper = estimate + z * se)
}

matrix_max_abs <- function(M) {
    if (!length(M)) 
        return(NA_real_)
    max(abs(M))
}

matrix_frobenius <- function(M) {
    sqrt(sum(M^2))
}

options(stringsAsFactors = FALSE, warn = 1)

script_directory <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg)) {
        return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]), 
            winslash = "/", mustWork = FALSE)))
    }
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

BASE_DIR <- script_directory()

ensure_directory <- function(path, 
    attempts = 5L, wait_seconds = 1) {
    path <- path.expand(path)
    for (attempt in seq_len(attempts)) {
        if (dir.exists(path)) 
            return(normalizePath(path, winslash = "/"))
        if (dir.create(path, recursive = TRUE, showWarnings = FALSE) && 
            dir.exists(path)) {
            return(normalizePath(path, winslash = "/"))
        }
        if (attempt < attempts) 
            Sys.sleep(wait_seconds)
    }
    stop("Could not create or access directory: ", path, call. = FALSE)
}

env_integer <- function(name, default) {
    value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
    if (is.na(value)) 
        stop(name, " must be an integer.", call. = FALSE)
    value
}

env_numeric <- function(name, default) {
    value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
    if (!is.finite(value)) 
        stop(name, " must be finite.", call. = FALSE)
    value
}

env_character <- function(name, default) {
    value <- trimws(Sys.getenv(name, unset = as.character(default)))
    if (!nzchar(value)) 
        stop(name, " must not be empty.", call. = FALSE)
    value
}

env_flag <- function(name, default = FALSE) {
    default_text <- if (isTRUE(default)) 
        "1"
    else "0"
    value <- tolower(trimws(Sys.getenv(name, unset = default_text)))
    if (!value %in% c("0", "1", "false", "true", "no", "yes")) {
        stop(name, " must be one of 0/1, false/true, or no/yes.", 
            call. = FALSE)
    }
    value %in% c("1", "true", "yes")
}

env_integer_vector <- function(name) {
    text <- trimws(Sys.getenv(name, unset = ""))
    if (!nzchar(text)) 
        return(integer())
    value <- suppressWarnings(as.integer(trimws(strsplit(text, 
        ",", fixed = TRUE)[[1L]])))
    if (anyNA(value)) 
        stop(name, " must contain comma-separated integers.", 
            call. = FALSE)
    unique(value)
}

append_csv <- function(data, path, attempts = 5L, wait_seconds = 1) {
    if (!nrow(data)) 
        return(invisible(NULL))
    parent <- ensure_directory(dirname(path))
    last_error <- NULL
    for (attempt in seq_len(attempts)) {
        tmp <- NULL
        ok <- tryCatch({
            exists <- file.exists(path)
            tmp <- tempfile(pattern = paste0(basename(path), 
                ".batch_"), tmpdir = parent, fileext = ".tmp")
            write.table(data, file = tmp, sep = ",", row.names = FALSE, 
                col.names = !exists, append = FALSE, quote = TRUE, 
                qmethod = "double")
            if (exists) {
                if (!file.append(path, tmp)) 
                  stop("file.append returned FALSE")
                unlink(tmp)
            }
            else if (!file.rename(tmp, path)) {
                if (!file.copy(tmp, path, overwrite = FALSE)) {
                  stop("Could not move the first checkpoint batch into place")
                }
                unlink(tmp)
            }
            TRUE
        }, error = function(e) {
            last_error <<- conditionMessage(e)
            if (!is.null(tmp) && file.exists(tmp)) 
                unlink(tmp)
            FALSE
        })
        if (isTRUE(ok)) 
            return(invisible(NULL))
        if (attempt < attempts) 
            Sys.sleep(wait_seconds)
    }
    stop("Failed to append checkpoint: ", path, if (!is.null(last_error)) 
        paste0("; ", last_error)
    else "", call. = FALSE)
}

safe_write_csv <- function(data, path) {
    parent <- ensure_directory(dirname(path))
    tmp <- tempfile(pattern = basename(path), tmpdir = parent, 
        fileext = ".tmp")
    write.csv(data, tmp, row.names = FALSE, na = "")
    if (!file.rename(tmp, path)) {
        if (!file.copy(tmp, path, overwrite = TRUE)) {
            unlink(tmp)
            stop("Could not write ", path, call. = FALSE)
        }
        unlink(tmp)
    }
    invisible(path)
}

format_elapsed <- function(seconds) {
    if (!is.finite(seconds)) 
        return("NA")
    if (seconds < 60) 
        return(sprintf("%.1f seconds", seconds))
    if (seconds < 3600) 
        return(sprintf("%.1f minutes", seconds/60))
    sprintf("%.2f hours", seconds/3600)
}

PROFILE <- tolower(trimws(Sys.getenv("PWRT_PROFILE", unset = "diagnostic")))

if (!PROFILE %in% c("smoke", "diagnostic", "manuscript")) {
        stop("PWRT_PROFILE must be smoke, diagnostic, or manuscript.", 
            call. = FALSE)
    }

profile_defaults <- switch(PROFILE, smoke = c(outer = 200L, 
        rerand = 199L, calibration = 2000L, batch = 20L), diagnostic = c(outer = 5000L, 
        rerand = 999L, calibration = 20000L, batch = 20L), manuscript = c(outer = 100000L, 
        rerand = 4999L, calibration = 100000L, batch = 10L))

direction_calibration_defaults <- switch(PROFILE, smoke = 5000L, 
        diagnostic = 50000L, manuscript = 200000L)

RUN_MODE <- tolower(trimws(Sys.getenv("PWRT_MODE", 
        unset = "audit")))

if (!RUN_MODE %in% c("audit", "calibrate", 
        "run", "aggregate", "validate")) {
        stop("PWRT_MODE must be audit, calibrate, run, aggregate, or validate.", 
            call. = FALSE)
    }

SCENARIO_SET <- tolower(trimws(Sys.getenv("PWRT_SCENARIO_SET", 
        unset = "supplemental")))

if (!SCENARIO_SET %in% c("supplemental", 
        "extended", "all")) {
        stop("PWRT_SCENARIO_SET must be supplemental, extended, or all.", 
            call. = FALSE)
    }

RUN_VERSION <- env_character("PWRT_RUN_VERSION", "v20260805_supplemental_v2")

MACHINE_ID <- env_character("PWRT_MACHINE_ID", Sys.info()[["nodename"]] %||% 
        "unknown_machine")

REQUIRE_EXISTING_CALIBRATION <- env_flag("PWRT_REQUIRE_EXISTING_CALIBRATION", 
        FALSE)

INCLUDE_ORACLE_N <- env_flag("PWRT_INCLUDE_ORACLE_N", 
        TRUE)

RT_SHAPE_DIAGNOSTICS <- env_flag("PWRT_RT_SHAPE_DIAGNOSTICS", 
        TRUE)

CALIBRATION_SCHEMA_VERSION <- "paircal_v2_main_design_sizes"

N_OUTER <- env_integer("PWRT_N_OUTER", profile_defaults["outer"])

N_RERANDOMIZATIONS <- env_integer("PWRT_N_RERAND", profile_defaults["rerand"])

N_CALIBRATION <- env_integer("PWRT_N_CALIBRATION", profile_defaults["calibration"])

N_DIRECTION_CALIBRATION <- env_integer("PWRT_N_DIRECTION_CALIBRATION", 
        direction_calibration_defaults)

CALIBRATION_BATCH <- env_integer("PWRT_CALIBRATION_BATCH", 
        1000L)

OUTER_BATCH <- env_integer("PWRT_OUTER_BATCH", 
        profile_defaults["batch"])

N_SHARDS <- env_integer("PWRT_N_SHARDS", 
        1L)

SHARD_ID <- env_integer("PWRT_SHARD_ID", 1L)

BASE_SEED <- env_integer("PWRT_SEED", 
        20260805L)

SCENARIO_FILTER <- env_integer_vector("PWRT_SCENARIO_IDS")

ALLOW_PARTIAL <- env_flag("PWRT_ALLOW_PARTIAL", FALSE)

EPSILON_EXPONENT <- env_numeric("PWRT_EPSILON_EXPONENT", 
        1)

STRICT_VALIDATION <- env_flag("PWRT_STRICT_VALIDATION", 
        FALSE)

VALIDATION_ABS_TOLERANCE <- env_numeric("PWRT_VALIDATION_ABS_TOLERANCE", 
        0.003)

NOMINAL_ABS_TOLERANCE <- env_numeric("PWRT_NOMINAL_ABS_TOLERANCE", 
        0.004)

VARIANCE_REL_TOLERANCE <- env_numeric("PWRT_VARIANCE_REL_TOLERANCE", 
        0.12)

physical_cores <- parallel::detectCores(logical = FALSE)

if (is.na(physical_cores)) physical_cores <- parallel::detectCores(logical = TRUE)

if (is.na(physical_cores)) physical_cores <- 1L

DEFAULT_CORES <- if (.Platform$OS.type == 
        "windows") {
        1L
    } else {
        max(1L, min(8L, physical_cores - 1L))
    }

N_CORES <- env_integer("PWRT_CORES", DEFAULT_CORES)

if (.Platform$OS.type == 
        "windows") N_CORES <- 1L

USE_PARALLEL <- .Platform$OS.type == 
        "unix" && N_CORES > 1L

PROJECT_DIR <- ensure_directory(Sys.getenv("PWRT_PROJECT_DIR", 
        unset = file.path(dirname(BASE_DIR), "pair_path_theory_output")))

OUTPUT_ROOT <- ensure_directory(Sys.getenv("PWRT_OUTPUT_DIR", 
        unset = file.path(PROJECT_DIR, "simulation_output")))

RUN_TAG <- paste0("profile_", PROFILE, "_set_", SCENARIO_SET, 
        "_M", N_OUTER, "_B", N_RERANDOMIZATIONS, "_B0", N_CALIBRATION, 
        "_Bdir", N_DIRECTION_CALIBRATION, "_eps", gsub("[.]", 
            "p", formatC(EPSILON_EXPONENT, format = "f", digits = 2)), 
        "_oracle", as.integer(INCLUDE_ORACLE_N), "_shape", as.integer(RT_SHAPE_DIAGNOSTICS), 
        "_seed", BASE_SEED, "_ver_", gsub("[^A-Za-z0-9._-]", 
            "-", RUN_VERSION), "_S", N_SHARDS)

OUTPUT_DIR <- ensure_directory(file.path(OUTPUT_ROOT, 
        RUN_TAG))

CALIBRATION_DIR <- ensure_directory(file.path(PROJECT_DIR, 
        "calibration_cache"))

SCENARIO_DIR <- ensure_directory(file.path(OUTPUT_DIR, 
        "scenario_shards"))

stopifnot(N_OUTER >= 1L, N_RERANDOMIZATIONS >= 
        1L, N_CALIBRATION >= 2L, N_DIRECTION_CALIBRATION >= 2L, 
        CALIBRATION_BATCH >= 1L, OUTER_BATCH >= 1L, N_SHARDS >= 
            1L, SHARD_ID >= 1L, SHARD_ID <= N_SHARDS, N_CORES >= 
            1L, EPSILON_EXPONENT > 0, VALIDATION_ABS_TOLERANCE > 
            0, NOMINAL_ABS_TOLERANCE > 0, VARIANCE_REL_TOLERANCE > 
            0)

if (PROFILE == "manuscript") {
        expected <- c(outer = 100000L, rerand = 4999L, calibration = 100000L, 
            direction_calibration = 200000L)
        observed <- c(outer = N_OUTER, rerand = N_RERANDOMIZATIONS, 
            calibration = N_CALIBRATION, direction_calibration = N_DIRECTION_CALIBRATION)
        if (any(expected != observed)) {
            warning("Manuscript profile was selected but settings differ from the recommended values: ", 
                paste(names(observed), observed, sep = "=", collapse = ", "))
        }
    }

seed_value <- function(scenario_id, replicate_id, stream_id) {
        modulus <- 2147483646
        value <- as.double(BASE_SEED) + 1000003 * as.double(stream_id) + 
            104729 * as.double(scenario_id) + 1009 * as.double(replicate_id)
        as.integer(value%%modulus + 1)
    }

calibration_seed <- function(design_id) {
        modulus <- 2147483646
        value <- as.double(BASE_SEED) + 7000003 + 99991 * as.double(design_id)
        as.integer(value%%modulus + 1)
    }

direction_calibration_seed <- function(design_id) {
        modulus <- 2147483646
        value <- as.double(BASE_SEED) + 17000009 + 199999 * as.double(design_id)
        as.integer(value%%modulus + 1)
    }

make_design_table <- function() {
        target_per_group <- c(100L, 500L, 200L, 1000L, 100L, 
            500L, 200L, 1000L)
        data.frame(design_id = 1:8, design_label = c("K2_npg100_n200_p080_independent", 
            "K2_npg500_n1000_p080_independent", "K5_npg200_n400_p080_independent", 
            "K5_npg1000_n2000_p080_independent", "K2_npg100_n200_p095_independent", 
            "K2_npg500_n1000_p095_independent", "K5_npg200_n400_p080_correlated", 
            "K5_npg1000_n2000_p080_correlated"), factor_count = c(2L, 
            2L, 5L, 5L, 2L, 2L, 5L, 5L), target_per_group = target_per_group, 
            total_n = 2L * target_per_group, pbc = c(0.8, 0.8, 
                0.8, 0.8, 0.95, 0.95, 0.8, 0.8), profile_type = c("independent", 
                "independent", "independent", "independent", 
                "independent", "independent", "latent_correlated", 
                "latent_correlated"), latent_strength = c(0, 
                0, 0, 0, 0, 0, 1, 1), set = c(rep("supplemental", 
                6L), rep("extended", 2L)), stringsAsFactors = FALSE)
    }

FACTOR_PREVALENCE_MASTER <- c(0.5, 0.4, 0.3, 0.2, 0.1)

prepare_design <- function(row) {
        K <- as.integer(row$factor_count)
        patterns <- all_binary_patterns(K)
        factor_prob <- FACTOR_PREVALENCE_MASTER[seq_len(K)]
        profile_prob <- if (row$profile_type == "independent") {
            profile_prob_independent(factor_prob, patterns)
        }
        else {
            profile_prob_latent_correlated(factor_prob = factor_prob, 
                latent_strength = row$latent_strength, mixing_prob = 0.5, 
                patterns = patterns)
        }
        list(design_id = as.integer(row$design_id), design_label = as.character(row$design_label), 
            factor_count = K, target_per_group = as.integer(row$target_per_group), 
            total_n = as.integer(row$total_n), pbc = as.numeric(row$pbc), 
            weights = rep(1, K + 1L), profile_type = as.character(row$profile_type), 
            latent_strength = as.numeric(row$latent_strength), 
            patterns = patterns, factor_prob = factor_prob, profile_prob = profile_prob, 
            marginal_prob = profile_marginals(patterns, profile_prob), 
            correlation = profile_correlations(patterns, profile_prob), 
            set = as.character(row$set))
    }

make_scenario_blueprints <- function() {
        supplemental <- data.frame(design_id = c(1, 1, 2, 2, 
            3, 3, 4, 4, 5, 5, 6, 6), outcome_type = rep(c("continuous", 
            "binary"), 6L), model_type = c("realistic", "common_log_odds", 
            "realistic", "common_log_odds", "realistic", "common_log_odds", 
            "realistic", "common_log_odds", "pure_symmetric", 
            "symmetric_stress", "pure_symmetric", "symmetric_stress"), 
            direction = c("first_factor_contrast", "common_log_odds", 
                "first_factor_contrast", "common_log_odds", "first_factor_contrast", 
                "common_log_odds", "first_factor_contrast", "common_log_odds", 
                "robust_max_ratio", "robust_max_ratio", "robust_max_ratio", 
                "robust_max_ratio"), scenario_class = c(rep("practical_heterogeneity", 
                8L), rep("strong_pair_path", 4L)), set = "supplemental", 
            target_max_abs_d = c(0.5, NA, 0.5, NA, 0.5, NA, 0.5, 
                NA, 1, 0.8, 1, 0.8), shared_noise_sd = c(NA, 
                NA, NA, NA, NA, NA, NA, NA, 0.05, NA, 0.05, NA), 
            outcome_sd = c(1, NA, 1, NA, 1, NA, 1, NA, NA, NA, 
                NA, NA), individual_effect_sd = c(0.25, NA, 0.25, 
                NA, 0.25, NA, 0.25, NA, NA, NA, NA, NA), max_abs_delta = c(rep(NA, 
                9L), 0.9, NA, 0.9), stringsAsFactors = FALSE)
        extended <- data.frame(design_id = c(7, 7, 8, 8), outcome_type = c("continuous", 
            "binary", "continuous", "binary"), model_type = c("pure_symmetric", 
            "symmetric_stress", "pure_symmetric", "symmetric_stress"), 
            direction = rep("basis_max_ratio", 4L), scenario_class = rep("correlated_factor_stress", 
                4L), set = "extended", target_max_abs_d = c(1, 
                0.8, 1, 0.8), shared_noise_sd = c(0.05, NA, 0.05, 
                NA), outcome_sd = rep(NA_real_, 4L), individual_effect_sd = rep(NA_real_, 
                4L), max_abs_delta = c(NA, 0.9, NA, 0.9), stringsAsFactors = FALSE)
        out <- rbind(supplemental, extended)
        out$scenario_id <- seq_len(nrow(out))
        out$scenario_code <- c(sprintf("S%02d", seq_len(nrow(supplemental))), 
            sprintf("E%02d", seq_len(nrow(extended))))
        out$boundary <- ifelse(out$outcome_type == "continuous", 
            -0.2, -0.1)
        out$alpha <- 0.025
        out$alternative <- "greater"
        out$scenario_label <- paste0(out$scenario_code, "_D", 
            out$design_id, "_", out$outcome_type, "_", out$model_type, 
            "_", out$direction)
        out
    }

DESIGN_TABLE <- make_design_table()

DESIGNS <- setNames(lapply(seq_len(nrow(DESIGN_TABLE)), 
        function(i) prepare_design(DESIGN_TABLE[i, ])), DESIGN_TABLE$design_id)

SCENARIOS <- make_scenario_blueprints()

if (SCENARIO_SET != 
        "all") {
        SCENARIOS <- SCENARIOS[SCENARIOS$set == SCENARIO_SET, 
            , drop = FALSE]
    }

if (length(SCENARIO_FILTER)) {
        SCENARIOS <- SCENARIOS[SCENARIOS$scenario_id %in% SCENARIO_FILTER, 
            , drop = FALSE]
    }

if (!nrow(SCENARIOS)) stop("No scenarios remain after filtering.", 
        call. = FALSE)

rownames(SCENARIOS) <- NULL

calibration_cache_file <- function(design) {
        file.path(CALIBRATION_DIR, paste0("pair_calibration_", 
            CALIBRATION_SCHEMA_VERSION, "_", design$design_label, 
            "_B0", N_CALIBRATION, "_seed", calibration_seed(design$design_id), 
            ".rds"))
    }

get_design_calibration <- function(design) {
        cache <- calibration_cache_file(design)
        if (file.exists(cache)) 
            return(readRDS(cache))
        if (isTRUE(REQUIRE_EXISTING_CALIBRATION)) {
            stop("Required analysis calibration is missing on worker machine: ", 
                cache, call. = FALSE)
        }
        lock <- paste0(cache, ".lock")
        have_lock <- dir.create(lock, showWarnings = FALSE)
        if (!have_lock) {
            for (attempt in seq_len(720L)) {
                if (file.exists(cache)) 
                  return(readRDS(cache))
                Sys.sleep(5)
                have_lock <- dir.create(lock, showWarnings = FALSE)
                if (have_lock) 
                  break
            }
        }
        if (!have_lock) 
            stop("Could not acquire calibration lock: ", lock, 
                call. = FALSE)
        on.exit(unlink(lock, recursive = TRUE, force = TRUE), 
            add = TRUE)
        if (file.exists(cache)) 
            return(readRDS(cache))
        message("Creating three-copy allocation calibration for ", 
            design$design_label, ": n=", design$total_n, ", B0=", 
            N_CALIBRATION)
        calibration <- calibrate_pair_path_design_R(B0 = N_CALIBRATION, 
            n = design$total_n, patterns = design$patterns, profile_prob = design$profile_prob, 
            pbc = design$pbc, weights_ = design$weights, seed = calibration_seed(design$design_id), 
            batch_size = CALIBRATION_BATCH, progress = TRUE)
        calibration$directions <- generalized_pair_directions(calibration$psi, 
            calibration$profile_prob)
        calibration$calibration_schema_version <- CALIBRATION_SCHEMA_VERSION
        calibration$calibration_role <- "analysis"
        tmp <- tempfile(pattern = "pair_calibration_", tmpdir = CALIBRATION_DIR, 
            fileext = ".rds")
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

USED_DESIGN_IDS <- sort(unique(SCENARIOS$design_id))

CALIBRATIONS <- setNames(lapply(USED_DESIGN_IDS, function(id) get_design_calibration(DESIGNS[[as.character(id)]])), 
        USED_DESIGN_IDS)

scenario_uses_calibrated_direction <- function(direction) {
        direction %in% c("robust_max_ratio", "basis_max_ratio")
    }

USED_DIRECTION_DESIGN_IDS <- sort(unique(SCENARIOS$design_id[vapply(SCENARIOS$direction, 
        scenario_uses_calibrated_direction, logical(1L))]))

direction_calibration_cache_file <- function(design) {
        file.path(CALIBRATION_DIR, paste0("direction_calibration_", 
            CALIBRATION_SCHEMA_VERSION, "_", design$design_label, 
            "_B0", N_DIRECTION_CALIBRATION, "_seed", direction_calibration_seed(design$design_id), 
            ".rds"))
    }

get_direction_calibration <- function(design) {
        cache <- direction_calibration_cache_file(design)
        if (file.exists(cache)) 
            return(readRDS(cache))
        if (isTRUE(REQUIRE_EXISTING_CALIBRATION)) {
            stop("Required direction-selection calibration is missing on worker machine: ", 
                cache, call. = FALSE)
        }
        lock <- paste0(cache, ".lock")
        have_lock <- dir.create(lock, showWarnings = FALSE)
        if (!have_lock) {
            for (attempt in seq_len(720L)) {
                if (file.exists(cache)) 
                  return(readRDS(cache))
                Sys.sleep(5)
                have_lock <- dir.create(lock, showWarnings = FALSE)
                if (have_lock) 
                  break
            }
        }
        if (!have_lock) 
            stop("Could not acquire direction-calibration lock: ", 
                lock, call. = FALSE)
        on.exit(unlink(lock, recursive = TRUE, force = TRUE), 
            add = TRUE)
        if (file.exists(cache)) 
            return(readRDS(cache))
        message("Creating independent direction-selection calibration for ", 
            design$design_label, ": n=", design$total_n, ", Bdir=", 
            N_DIRECTION_CALIBRATION)
        calibration <- calibrate_pair_path_design_R(B0 = N_DIRECTION_CALIBRATION, 
            n = design$total_n, patterns = design$patterns, profile_prob = design$profile_prob, 
            pbc = design$pbc, weights_ = design$weights, seed = direction_calibration_seed(design$design_id), 
            batch_size = CALIBRATION_BATCH, progress = TRUE)
        calibration$directions <- generalized_pair_directions(calibration$psi, 
            calibration$profile_prob)
        calibration$calibration_role <- "direction_selection"
        calibration$calibration_schema_version <- CALIBRATION_SCHEMA_VERSION
        tmp <- tempfile(pattern = "direction_calibration_", tmpdir = CALIBRATION_DIR, 
            fileext = ".rds")
        saveRDS(calibration, tmp)
        if (!file.rename(tmp, cache)) {
            if (!file.copy(tmp, cache, overwrite = TRUE)) {
                unlink(tmp)
                stop("Could not save direction calibration: ", 
                  cache, call. = FALSE)
            }
            unlink(tmp)
        }
        calibration
    }

DIRECTION_CALIBRATIONS <- setNames(lapply(USED_DIRECTION_DESIGN_IDS, 
        function(id) {
            get_direction_calibration(DESIGNS[[as.character(id)]])
        }), USED_DIRECTION_DESIGN_IDS)

first_factor_contrast_direction <- function(patterns, 
        profile_prob) {
        patterns <- check_binary_factor_matrix(patterns)
        profile_prob <- as.numeric(profile_prob)
        direction <- weighted_center(patterns[, 1L], profile_prob)
        scale_direction(direction, 1)
    }

basis_generalized_pair_direction <- function(psi, profile_prob, 
        patterns, which = c("max", "min")) {
        which <- match.arg(which)
        psi <- nearest_psd(psi)
        patterns <- check_binary_factor_matrix(patterns)
        profile_prob <- as.numeric(profile_prob)
        H <- apply(patterns, 2L, weighted_center, weight = profile_prob)
        H <- as.matrix(H)
        PiH <- profile_prob * H
        C <- crossprod(H, PiH)
        A <- crossprod(H, psi %*% H)
        cc <- eigen((C + t(C))/2, symmetric = TRUE)
        tolerance <- 1e-10 * max(1, max(abs(cc$values)))
        keep <- cc$values > tolerance
        if (!any(keep)) 
            stop("The prespecified direction basis has zero rank.", 
                call. = FALSE)
        whitening <- cc$vectors[, keep, drop = FALSE] %*% diag(1/sqrt(cc$values[keep]), 
            nrow = sum(keep))
        M <- crossprod(whitening, A %*% whitening)
        ee <- eigen((M + t(M))/2, symmetric = TRUE)
        index <- if (which == "max") 
            which.max(ee$values)
        else which.min(ee$values)
        coefficient <- drop(whitening %*% ee$vectors[, index])
        direction <- drop(H %*% coefficient)
        direction <- weighted_center(direction, profile_prob)
        scale_direction(direction, 1)
    }

weighted_direction_cosine <- function(x, y, weight) {
        x <- as.numeric(x)
        y <- as.numeric(y)
        weight <- as.numeric(weight)
        denominator <- sqrt(sum(weight * x^2) * sum(weight * 
            y^2))
        if (!is.finite(denominator) || denominator <= 0) 
            return(NA_real_)
        abs(sum(weight * x * y)/denominator)
    }

direction_for_scenario <- function(scenario, design, analysis_calibration) {
        direction_calibration <- DIRECTION_CALIBRATIONS[[as.character(scenario$design_id)]]
        switch(scenario$direction, zero = rep(0, analysis_calibration$J), 
            first_factor_contrast = first_factor_contrast_direction(design$patterns, 
                design$profile_prob), robust_max_ratio = direction_calibration$directions$max_direction, 
            basis_max_ratio = basis_generalized_pair_direction(direction_calibration$psi, 
                design$profile_prob, design$patterns, which = "max"), 
            common_log_odds = rep(NA_real_, analysis_calibration$J), 
            stop("Unknown direction: ", scenario$direction, call. = FALSE))
    }

build_scenario_model <- function(scenario) {
        design <- DESIGNS[[as.character(scenario$design_id)]]
        calibration <- CALIBRATIONS[[as.character(scenario$design_id)]]
        direction <- direction_for_scenario(scenario, design, 
            calibration)
        boundary <- scenario$boundary
        value_or <- function(x, fallback) {
            x <- suppressWarnings(as.numeric(x))
            if (length(x) == 1L && is.finite(x)) 
                x
            else fallback
        }
        model <- if (scenario$outcome_type == "continuous") {
            switch(scenario$model_type, homogeneous = make_continuous_homogeneous_model(design$patterns, 
                design$profile_prob, boundary, outcome_sd = value_or(scenario$outcome_sd, 
                  1)), pure_symmetric = make_continuous_pure_model(design$patterns, 
                design$profile_prob, boundary, direction, target_max_abs_d = value_or(scenario$target_max_abs_d, 
                  1), shared_noise_sd = value_or(scenario$shared_noise_sd, 
                  0.05)), realistic = make_continuous_realistic_model(design$patterns, 
                design$profile_prob, boundary, direction, target_max_abs_d = value_or(scenario$target_max_abs_d, 
                  0.5), outcome_sd = value_or(scenario$outcome_sd, 
                  1), individual_effect_sd = value_or(scenario$individual_effect_sd, 
                  0.25)), stop("Unknown continuous model type: ", 
                scenario$model_type, call. = FALSE))
        }
        else {
            switch(scenario$model_type, homogeneous = make_binary_homogeneous_model(design$patterns, 
                design$profile_prob, boundary, target_control_risk = 0.6), 
                symmetric_stress = make_binary_symmetric_model(design$patterns, 
                  design$profile_prob, boundary, direction, target_max_abs_d = value_or(scenario$target_max_abs_d, 
                    0.8), max_abs_delta = value_or(scenario$max_abs_delta, 
                    0.9)), common_log_odds = calibrate_common_log_odds_model(design$patterns, 
                  design$profile_prob, target_control_risk = 0.6, 
                  target_risk_difference = boundary), stop("Unknown binary model type: ", 
                  scenario$model_type, call. = FALSE))
        }
        d <- as.numeric(model$d)
        weighted_mean_d <- sum(design$profile_prob * d)
        if (abs(weighted_mean_d) > 1e-10) {
            stop("Scenario ", scenario$scenario_id, " does not satisfy the marginal boundary: weighted mean d=", 
                weighted_mean_d, call. = FALSE)
        }
        d_pi <- sum(design$profile_prob * d^2)
        d_psi <- drop(crossprod(d, calibration$psi %*% d))
        ratio <- if (d_pi > 0) 
            d_psi/d_pi
        else 1
        direction_calibration <- DIRECTION_CALIBRATIONS[[as.character(scenario$design_id)]]
        d_psi_direction <- if (!is.null(direction_calibration)) {
            drop(crossprod(d, direction_calibration$psi %*% d))
        }
        else d_psi
        ratio_direction <- if (d_pi > 0) 
            d_psi_direction/d_pi
        else 1
        model$d_pi <- d_pi
        model$d_psi <- d_psi
        model$d_psi_direction <- d_psi_direction
        model$pair_ratio <- ratio
        model$pair_ratio_direction <- ratio_direction
        model$pair_scalar_gap <- d_psi - d_pi
        model
    }

SCENARIO_MODELS <- setNames(lapply(seq_len(nrow(SCENARIOS)), 
        function(i) build_scenario_model(SCENARIOS[i, ])), SCENARIOS$scenario_id)

build_design_diagnostics <- function() {
        rows <- lapply(USED_DESIGN_IDS, function(id) {
            design <- DESIGNS[[as.character(id)]]
            calibration <- CALIBRATIONS[[as.character(id)]]
            direction_calibration <- DIRECTION_CALIBRATIONS[[as.character(id)]]
            corr <- design$correlation
            offdiag <- corr[row(corr) != col(corr)]
            direction_fields <- if (is.null(direction_calibration)) {
                list(direction_calibration_paths = NA_integer_, 
                  direction_min_ratio = NA_real_, direction_max_ratio = NA_real_, 
                  analysis_min_ratio = calibration$directions$min_ratio, 
                  analysis_max_ratio = calibration$directions$max_ratio, 
                  min_direction_weighted_cosine = NA_real_, max_direction_weighted_cosine = NA_real_, 
                  direction_min_eigengap = NA_real_, direction_max_eigengap = NA_real_)
            }
            else {
                ratios <- direction_calibration$directions$all_ratios
                list(direction_calibration_paths = direction_calibration$B0, 
                  direction_min_ratio = direction_calibration$directions$min_ratio, 
                  direction_max_ratio = direction_calibration$directions$max_ratio, 
                  analysis_min_ratio = calibration$directions$min_ratio, 
                  analysis_max_ratio = calibration$directions$max_ratio, 
                  min_direction_weighted_cosine = weighted_direction_cosine(direction_calibration$directions$min_direction, 
                    calibration$directions$min_direction, design$profile_prob), 
                  max_direction_weighted_cosine = weighted_direction_cosine(direction_calibration$directions$max_direction, 
                    calibration$directions$max_direction, design$profile_prob), 
                  direction_min_eigengap = if (length(ratios) >= 
                    2L) ratios[2L] - ratios[1L] else NA_real_, 
                  direction_max_eigengap = if (length(ratios) >= 
                    2L) {
                    ratios[length(ratios)] - ratios[length(ratios) - 
                      1L]
                  } else NA_real_)
            }
            data.frame(design_id = design$design_id, design_label = design$design_label, 
                factor_count = design$factor_count, target_per_group = design$target_per_group, 
                total_n = design$total_n, pbc = design$pbc, profile_type = design$profile_type, 
                latent_strength = design$latent_strength, analysis_calibration_paths = calibration$B0, 
                analysis_calibration_seconds = calibration$elapsed_seconds, 
                direction_calibration_paths = direction_fields$direction_calibration_paths, 
                direction_min_ratio = direction_fields$direction_min_ratio, 
                direction_max_ratio = direction_fields$direction_max_ratio, 
                analysis_min_ratio = direction_fields$analysis_min_ratio, 
                analysis_max_ratio = direction_fields$analysis_max_ratio, 
                min_direction_weighted_cosine = direction_fields$min_direction_weighted_cosine, 
                max_direction_weighted_cosine = direction_fields$max_direction_weighted_cosine, 
                direction_min_eigengap = direction_fields$direction_min_eigengap, 
                direction_max_eigengap = direction_fields$direction_max_eigengap, 
                max_abs_mean_u = max(abs(calibration$mean_u)), 
                max_abs_mean_g01 = max(abs(calibration$mean_g01)), 
                max_abs_mean_g02 = max(abs(calibration$mean_g02)), 
                max_abs_cross_D01 = matrix_max_abs(calibration$cross_d01), 
                max_abs_cross_D02 = matrix_max_abs(calibration$cross_d02), 
                max_abs_cross_D12 = matrix_max_abs(calibration$cross_d12), 
                max_abs_cross_G0102 = matrix_max_abs(calibration$cross_g0102), 
                frobenius_cross_G0102 = matrix_frobenius(calibration$cross_g0102), 
                max_abs_factor_correlation = if (length(offdiag)) 
                  max(abs(offdiag))
                else 0, stringsAsFactors = FALSE)
        })
        do.call(rbind, rows)
    }

build_scenario_definition_table <- function() {
        rows <- lapply(seq_len(nrow(SCENARIOS)), function(i) {
            scenario <- SCENARIOS[i, ]
            design <- DESIGNS[[as.character(scenario$design_id)]]
            model <- SCENARIO_MODELS[[as.character(scenario$scenario_id)]]
            expected_empty <- sum((1 - design$profile_prob)^design$total_n)
            min_expected_count <- design$total_n * min(design$profile_prob)
            p0_min <- if (!is.null(model$p0)) 
                min(model$p0)
            else NA_real_
            p0_max <- if (!is.null(model$p0)) 
                max(model$p0)
            else NA_real_
            p1_min <- if (!is.null(model$p1)) 
                min(model$p1)
            else NA_real_
            p1_max <- if (!is.null(model$p1)) 
                max(model$p1)
            else NA_real_
            data.frame(scenario_id = scenario$scenario_id, scenario_code = scenario$scenario_code, 
                scenario_label = scenario$scenario_label, scenario_class = scenario$scenario_class, 
                set = scenario$set, design_id = scenario$design_id, 
                design_label = design$design_label, factor_count = design$factor_count, 
                target_per_group = design$target_per_group, joint_strata = nrow(design$patterns), 
                total_n = design$total_n, pbc = design$pbc, profile_type = design$profile_type, 
                min_profile_probability = min(design$profile_prob), 
                min_expected_stratum_count = min_expected_count, 
                expected_empty_strata = expected_empty, outcome_type = scenario$outcome_type, 
                model_type = scenario$model_type, direction = scenario$direction, 
                boundary = scenario$boundary, alpha = scenario$alpha, 
                achieved_effect = model$achieved_effect, min_d = min(model$d), 
                max_d = max(model$d), weighted_mean_d = sum(design$profile_prob * 
                  model$d), d_Pi_d = model$d_pi, d_Psi_d_analysis = model$d_psi, 
                d_Psi_d_direction = model$d_psi_direction, pair_ratio_analysis = model$pair_ratio, 
                pair_ratio_direction = model$pair_ratio_direction, 
                pair_ratio_cross_calibration_difference = model$pair_ratio - 
                  model$pair_ratio_direction, pair_scalar_gap = model$pair_scalar_gap, 
                p0_min = p0_min, p0_max = p0_max, p1_min = p1_min, 
                p1_max = p1_max, target_max_abs_d = scenario$target_max_abs_d, 
                shared_noise_sd = scenario$shared_noise_sd, outcome_sd = scenario$outcome_sd, 
                individual_effect_sd = scenario$individual_effect_sd, 
                max_abs_delta = scenario$max_abs_delta, achieved_control_risk = model$achieved_control_risk %||% 
                  NA_real_, achieved_treatment_risk = model$achieved_treatment_risk %||% 
                  NA_real_, common_log_odds_theta = model$theta %||% 
                  NA_real_, audit_boundary_ok = abs(sum(design$profile_prob * 
                  model$d)) <= 1e-10, audit_binary_probabilities_ok = if (scenario$outcome_type == 
                  "binary") {
                  p0_min > 0 && p0_max < 1 && p1_min > 0 && p1_max < 
                    1
                }
                else TRUE, audit_direction_stability_ok = if (scenario$direction %in% 
                  c("robust_max_ratio", "basis_max_ratio")) {
                  tolerance <- if (scenario$set == "supplemental") 
                    0.02
                  else 0.03
                  abs(model$pair_ratio - model$pair_ratio_direction) <= 
                    tolerance
                }
                else TRUE, stringsAsFactors = FALSE)
        })
        do.call(rbind, rows)
    }

DESIGN_DIAGNOSTICS <- build_design_diagnostics()

SCENARIO_DEFINITIONS <- build_scenario_definition_table()

safe_write_csv(DESIGN_DIAGNOSTICS, file.path(OUTPUT_DIR, 
        "design_pair_path_diagnostics.csv"))

safe_write_csv(SCENARIO_DEFINITIONS, 
        file.path(OUTPUT_DIR, "scenario_definitions.csv"))

if (RUN_MODE == 
        "audit") {
        if (!all(SCENARIO_DEFINITIONS$audit_boundary_ok)) {
            stop("At least one scenario fails the marginal-boundary audit.", 
                call. = FALSE)
        }
        if (!all(SCENARIO_DEFINITIONS$audit_binary_probabilities_ok)) {
            stop("At least one binary scenario has invalid probabilities.", 
                call. = FALSE)
        }
        if (!all(SCENARIO_DEFINITIONS$audit_direction_stability_ok)) {
            stop("A calibrated treatment-effect direction is unstable across independent calibrations.", 
                call. = FALSE)
        }
        message("Design and scenario audits passed. Files were written to: ", 
            OUTPUT_DIR)
        quit(save = "no", status = 0L)
    }

if (RUN_MODE == "calibrate") {
        message("Calibration and scenario-definition files were written to: ", 
            OUTPUT_DIR)
        quit(save = "no", status = 0L)
    }

simulate_trial <- function(scenario, replicate_id, design, 
        model) {
        set.seed(seed_value(scenario$scenario_id, replicate_id, 
            1L))
        profile <- generate_profile_sequence(n = design$total_n, 
            patterns = design$patterns, profile_prob = design$profile_prob)
        z <- ps_assign_R(X = profile$X, pbc = design$pbc, weights_ = design$weights, 
            seed = seed_value(scenario$scenario_id, replicate_id, 
                2L))
        A <- as.integer((z + 1L)/2L)
        set.seed(seed_value(scenario$scenario_id, replicate_id, 
            3L))
        y <- generate_outcome_from_model(profile$id, A, z, model)
        list(profile_id = profile$id, X = profile$X, z = as.integer(z), 
            A = A, y = y)
    }

analyze_one_score <- function(score, trial, scenario, 
        calibration, d) {
        statistic <- 0.5 * sum(trial$z * score)
        sampling <- siga_sampling_variance(score, trial$X, calibration)
        randomization <- siga_randomization_variance(sampling, 
            d, calibration, epsilon_n = calibration$n^(-EPSILON_EXPONENT))
        p_s <- gaussian_pvalue(statistic, sampling$variance, 
            alternative = scenario$alternative)
        p_r <- gaussian_pvalue(statistic, randomization$variance, 
            alternative = scenario$alternative)
        list(statistic = statistic, sampling_variance = sampling$variance, 
            randomization_variance = randomization$variance, 
            randomization_variance_raw = randomization$variance_raw, 
            pair_correction = randomization$correction, pair_scalar_gap_realized = randomization$scalar_gap, 
            p_s = p_s$p, z_s = p_s$z, p_r = p_r$p, z_r = p_r$z, 
            kappa = sampling$kappa, kappa_raw = sampling$kappa_raw, 
            kappa_truncated = sampling$kappa_raw < 0, randomization_variance_floor = randomization$variance_floor, 
            randomization_variance_truncated = randomization$truncated, 
            safeguard_active = randomization$safeguard_active, 
            low_component = sampling$low_component, within_component = sampling$within_component, 
            counts = sampling$decomposition$counts)
    }

run_replicate <- function(scenario, replicate_id, calibration, 
        design, model) {
        generation_start <- proc.time()[3L]
        trial <- simulate_trial(scenario, replicate_id, design, 
            model)
        generation_seconds <- proc.time()[3L] - generation_start
        score_start <- proc.time()[3L]
        scores <- make_boundary_scores(y = trial$y, A = trial$A, 
            X = trial$X, boundary = scenario$boundary)
        score_seconds <- proc.time()[3L] - score_start
        siga_start <- proc.time()[3L]
        unadjusted <- analyze_one_score(score = scores$unadjusted, 
            trial = trial, scenario = scenario, calibration = calibration, 
            d = model$d)
        adjusted <- analyze_one_score(score = scores$adjusted, 
            trial = trial, scenario = scenario, calibration = calibration, 
            d = model$d)
        siga_seconds <- proc.time()[3L] - siga_start
        rt_start <- proc.time()[3L]
        rt <- rt_pvalues_R(X = trial$X, scores = cbind(unadjusted = scores$unadjusted, 
            adjusted = scores$adjusted), z_obs = trial$z, B = N_RERANDOMIZATIONS, 
            pbc = design$pbc, weights_ = design$weights, seed = seed_value(scenario$scenario_id, 
                replicate_id, 4L), return_randomization_moments = TRUE, 
            return_shape_diagnostics = RT_SHAPE_DIAGNOSTICS)
        rt_seconds <- proc.time()[3L] - rt_start
        p_rt_u <- rt$greater[1L]
        p_rt_a <- rt$greater[2L]
        rt_var_u <- rt$randomization_variance[1L]
        rt_var_a <- rt$randomization_variance[2L]
        rt_mean_u <- rt$randomization_mean[1L]
        rt_mean_a <- rt$randomization_mean[2L]
        rt_skew_u <- if (!is.null(rt$randomization_skewness)) 
            rt$randomization_skewness[1L]
        else NA_real_
        rt_skew_a <- if (!is.null(rt$randomization_skewness)) 
            rt$randomization_skewness[2L]
        else NA_real_
        rt_kurt_u <- if (!is.null(rt$randomization_excess_kurtosis)) 
            rt$randomization_excess_kurtosis[1L]
        else NA_real_
        rt_kurt_a <- if (!is.null(rt$randomization_excess_kurtosis)) 
            rt$randomization_excess_kurtosis[2L]
        else NA_real_
        oracle_u <- if (isTRUE(INCLUDE_ORACLE_N) && is.finite(rt_var_u) && 
            rt_var_u > 0) {
            gaussian_pvalue(unadjusted$statistic - rt_mean_u, 
                rt_var_u, alternative = scenario$alternative)
        }
        else list(z = NA_real_, p = NA_real_)
        oracle_a <- if (isTRUE(INCLUDE_ORACLE_N) && is.finite(rt_var_a) && 
            rt_var_a > 0) {
            gaussian_pvalue(adjusted$statistic - rt_mean_a, rt_var_a, 
                alternative = scenario$alternative)
        }
        else list(z = NA_real_, p = NA_real_)
        data.frame(scenario_id = scenario$scenario_id, scenario_code = scenario$scenario_code, 
            scenario_label = scenario$scenario_label, scenario_class = scenario$scenario_class, 
            replicate = replicate_id, set = scenario$set, design_id = scenario$design_id, 
            design_label = design$design_label, factor_count = design$factor_count, 
            target_per_group = design$target_per_group, total_n = design$total_n, 
            pbc = design$pbc, profile_type = design$profile_type, 
            outcome_type = scenario$outcome_type, model_type = scenario$model_type, 
            direction = scenario$direction, boundary = scenario$boundary, 
            alpha = scenario$alpha, achieved_effect = model$achieved_effect, 
            pair_ratio_model = model$pair_ratio, pair_ratio_direction_calibration = model$pair_ratio_direction, 
            pair_scalar_gap_model = model$pair_scalar_gap, configured_outer_trials = N_OUTER, 
            rerandomizations_per_trial = N_RERANDOMIZATIONS, 
            allocation_calibration_paths = N_CALIBRATION, direction_calibration_paths = N_DIRECTION_CALIBRATION, 
            epsilon_exponent = EPSILON_EXPONENT, configured_shards = N_SHARDS, 
            base_seed = BASE_SEED, run_version = RUN_VERSION, 
            machine_id = MACHINE_ID, include_oracle_n = INCLUDE_ORACLE_N, 
            rt_shape_diagnostics = RT_SHAPE_DIAGNOSTICS, calibration_schema_version = CALIBRATION_SCHEMA_VERSION, 
            statistic_unadjusted = unadjusted$statistic, statistic_adjusted = adjusted$statistic, 
            p_siga_s_unadjusted = unadjusted$p_s, p_siga_r_unadjusted = unadjusted$p_r, 
            p_oracle_n_unadjusted = oracle_u$p, p_rt_unadjusted = p_rt_u, 
            p_siga_s_adjusted = adjusted$p_s, p_siga_r_adjusted = adjusted$p_r, 
            p_oracle_n_adjusted = oracle_a$p, p_rt_adjusted = p_rt_a, 
            reject_siga_s_unadjusted = unadjusted$p_s <= scenario$alpha, 
            reject_siga_r_unadjusted = unadjusted$p_r <= scenario$alpha, 
            reject_oracle_n_unadjusted = if (isTRUE(INCLUDE_ORACLE_N)) 
                oracle_u$p <= scenario$alpha
            else NA, reject_rt_unadjusted = p_rt_u <= scenario$alpha, 
            reject_siga_s_adjusted = adjusted$p_s <= scenario$alpha, 
            reject_siga_r_adjusted = adjusted$p_r <= scenario$alpha, 
            reject_oracle_n_adjusted = if (isTRUE(INCLUDE_ORACLE_N)) 
                oracle_a$p <= scenario$alpha
            else NA, reject_rt_adjusted = p_rt_a <= scenario$alpha, 
            variance_s_unadjusted = unadjusted$sampling_variance, 
            variance_r_unadjusted = unadjusted$randomization_variance, 
            variance_rt_unadjusted = rt_var_u, variance_s_adjusted = adjusted$sampling_variance, 
            variance_r_adjusted = adjusted$randomization_variance, 
            variance_rt_adjusted = rt_var_a, pair_correction_unadjusted = unadjusted$pair_correction, 
            pair_correction_adjusted = adjusted$pair_correction, 
            pair_scalar_gap_realized = unadjusted$pair_scalar_gap_realized, 
            randomization_mean_unadjusted = rt_mean_u, randomization_mean_adjusted = rt_mean_a, 
            randomization_skewness_unadjusted = rt_skew_u, randomization_skewness_adjusted = rt_skew_a, 
            randomization_excess_kurtosis_unadjusted = rt_kurt_u, 
            randomization_excess_kurtosis_adjusted = rt_kurt_a, 
            oracle_z_unadjusted = oracle_u$z, oracle_z_adjusted = oracle_a$z, 
            kappa_raw_unadjusted = unadjusted$kappa_raw, kappa_raw_adjusted = adjusted$kappa_raw, 
            kappa_truncated_unadjusted = unadjusted$kappa_truncated, 
            kappa_truncated_adjusted = adjusted$kappa_truncated, 
            variance_r_floor_unadjusted = unadjusted$randomization_variance_floor, 
            variance_r_floor_adjusted = adjusted$randomization_variance_floor, 
            variance_r_truncated_unadjusted = unadjusted$randomization_variance_truncated, 
            variance_r_truncated_adjusted = adjusted$randomization_variance_truncated, 
            safeguard_active_unadjusted = unadjusted$safeguard_active, 
            safeguard_active_adjusted = adjusted$safeguard_active, 
            low_component_unadjusted = unadjusted$low_component, 
            within_component_unadjusted = unadjusted$within_component, 
            low_component_adjusted = adjusted$low_component, 
            within_component_adjusted = adjusted$within_component, 
            empty_joint_strata = sum(unadjusted$counts == 0L), 
            min_positive_joint_stratum_count = if (any(unadjusted$counts > 
                0L)) {
                min(unadjusted$counts[unadjusted$counts > 0L])
            }
            else 0L, max_joint_stratum_count = max(unadjusted$counts), 
            observed_treated = sum(trial$A), data_generation_seconds = generation_seconds, 
            score_construction_seconds = score_seconds, siga_analysis_seconds = siga_seconds, 
            rt_analysis_seconds = rt_seconds, stringsAsFactors = FALSE)
    }

scenario_paths <- function(scenario) {
        directory <- ensure_directory(file.path(SCENARIO_DIR, 
            sprintf("scenario_%02d_%s", scenario$scenario_id, 
                scenario$scenario_label)))
        list(directory = directory, shard = file.path(directory, 
            sprintf("shard_%03d_of_%03d.csv", SHARD_ID, N_SHARDS)))
    }

read_completed_replicates <- function(path) {
        if (!file.exists(path) || file.info(path)$size <= 0) 
            return(integer())
        dat <- tryCatch(read.csv(path, stringsAsFactors = FALSE), 
            error = function(e) NULL)
        if (is.null(dat) || !"replicate" %in% names(dat)) 
            return(integer())
        unique(as.integer(dat$replicate))
    }

run_scenario_shard <- function(scenario) {
        design <- DESIGNS[[as.character(scenario$design_id)]]
        calibration <- CALIBRATIONS[[as.character(scenario$design_id)]]
        model <- SCENARIO_MODELS[[as.character(scenario$scenario_id)]]
        paths <- scenario_paths(scenario)
        assigned <- seq.int(from = SHARD_ID, to = N_OUTER, by = N_SHARDS)
        completed <- read_completed_replicates(paths$shard)
        pending <- setdiff(assigned, completed)
        message("Scenario ", scenario$scenario_id, " [", scenario$scenario_label, 
            "]", ": shard ", SHARD_ID, "/", N_SHARDS, "; assigned=", 
            length(assigned), "; completed=", length(completed), 
            "; pending=", length(pending))
        if (!length(pending)) 
            return(invisible(NULL))
        start <- proc.time()[3L]
        for (from in seq.int(1L, length(pending), by = OUTER_BATCH)) {
            ids <- pending[from:min(length(pending), from + OUTER_BATCH - 
                1L)]
            worker <- function(id) {
                run_replicate(scenario, id, calibration, design, 
                  model)
            }
            rows <- if (USE_PARALLEL && length(ids) > 1L) {
                parallel::mclapply(ids, worker, mc.cores = min(N_CORES, 
                  length(ids)))
            }
            else {
                lapply(ids, worker)
            }
            batch <- do.call(rbind, rows)
            append_csv(batch, paths$shard)
            completed_now <- length(completed) + min(length(pending), 
                from + length(ids) - 1L)
            message("  completed ", completed_now, "/", length(assigned), 
                " assigned trials; elapsed ", format_elapsed(proc.time()[3L] - 
                  start))
        }
        invisible(NULL)
    }

read_scenario_details <- function(scenario) {
        pattern <- sprintf("scenario_%02d_%s", scenario$scenario_id, 
            scenario$scenario_label)
        directory <- file.path(SCENARIO_DIR, pattern)
        if (!dir.exists(directory)) 
            stop("Scenario directory is missing: ", directory, 
                call. = FALSE)
        shard_pattern <- sprintf("^shard_[0-9]+_of_%03d\\.csv$", 
            N_SHARDS)
        files <- list.files(directory, pattern = shard_pattern, 
            full.names = TRUE)
        if (!length(files)) 
            stop("No shard files found for ", scenario$scenario_label, 
                call. = FALSE)
        pieces <- lapply(files, function(path) {
            dat <- read.csv(path, stringsAsFactors = FALSE)
            dat$source_file <- basename(path)
            dat
        })
        dat <- do.call(rbind, pieces)
        required_metadata <- c("scenario_id", "configured_outer_trials", 
            "rerandomizations_per_trial", "allocation_calibration_paths", 
            "direction_calibration_paths", "epsilon_exponent", 
            "configured_shards", "base_seed", "run_version", 
            "include_oracle_n", "rt_shape_diagnostics", "calibration_schema_version", 
            "machine_id")
        missing_metadata <- setdiff(required_metadata, names(dat))
        if (length(missing_metadata)) {
            stop("Missing metadata columns: ", paste(missing_metadata, 
                collapse = ", "), call. = FALSE)
        }
        expected_numeric_metadata <- list(scenario_id = scenario$scenario_id, 
            configured_outer_trials = N_OUTER, rerandomizations_per_trial = N_RERANDOMIZATIONS, 
            allocation_calibration_paths = N_CALIBRATION, direction_calibration_paths = N_DIRECTION_CALIBRATION, 
            epsilon_exponent = EPSILON_EXPONENT, configured_shards = N_SHARDS, 
            base_seed = BASE_SEED, include_oracle_n = as.integer(INCLUDE_ORACLE_N), 
            rt_shape_diagnostics = as.integer(RT_SHAPE_DIAGNOSTICS))
        for (field in names(expected_numeric_metadata)) {
            observed <- unique(as.numeric(dat[[field]]))
            expected <- as.numeric(expected_numeric_metadata[[field]])
            if (length(observed) != 1L || !isTRUE(all.equal(observed, 
                expected))) {
                stop("Incompatible metadata in ", scenario$scenario_label, 
                  ": ", field, "=", paste(observed, collapse = ","), 
                  "; expected ", expected, call. = FALSE)
            }
        }
        expected_character_metadata <- list(run_version = RUN_VERSION, 
            calibration_schema_version = CALIBRATION_SCHEMA_VERSION)
        for (field in names(expected_character_metadata)) {
            observed <- unique(as.character(dat[[field]]))
            expected <- as.character(expected_character_metadata[[field]])
            if (length(observed) != 1L || !identical(observed, 
                expected)) {
                stop("Incompatible metadata in ", scenario$scenario_label, 
                  ": ", field, "=", paste(observed, collapse = ","), 
                  "; expected ", expected, call. = FALSE)
            }
        }
        dat$replicate <- as.integer(dat$replicate)
        dat <- dat[order(dat$replicate), , drop = FALSE]
        duplicate_id <- duplicated(dat$replicate)
        if (any(duplicate_id)) {
            duplicates <- unique(dat$replicate[duplicate_id])
            for (id in duplicates) {
                block <- dat[dat$replicate == id, setdiff(names(dat), 
                  "source_file"), drop = FALSE]
                if (nrow(unique(block)) > 1L) {
                  stop("Conflicting duplicate replicate ", id, 
                    " in ", scenario$scenario_label, call. = FALSE)
                }
            }
            dat <- dat[!duplicated(dat$replicate), , drop = FALSE]
        }
        if (!ALLOW_PARTIAL && nrow(dat) != N_OUTER) {
            stop(scenario$scenario_label, " contains ", nrow(dat), 
                " unique replicates; expected ", N_OUTER, ".", 
                call. = FALSE)
        }
        dat
    }

rejection_summary <- function(x) {
        x <- as.logical(x)
        x <- x[!is.na(x)]
        if (!length(x)) {
            return(c(probability = NA_real_, lower = NA_real_, 
                upper = NA_real_))
        }
        ci <- wilson_interval(sum(x), length(x))
        c(probability = unname(ci["estimate"]), lower = unname(ci["lower"]), 
            upper = unname(ci["upper"]))
    }

safe_ratio <- function(numerator, denominator) {
        out <- numerator/denominator
        out[!is.finite(out)] <- NA_real_
        out
    }

central_moment_summary <- function(x) {
        x <- as.numeric(x)
        x <- x[is.finite(x)]
        if (length(x) < 4L) {
            return(c(mean = NA_real_, variance = NA_real_, skewness = NA_real_, 
                excess_kurtosis = NA_real_))
        }
        center <- x - mean(x)
        variance <- mean(center^2)
        if (!is.finite(variance) || variance <= 0) {
            return(c(mean = mean(x), variance = variance, skewness = NA_real_, 
                excess_kurtosis = NA_real_))
        }
        c(mean = mean(x), variance = var(x), skewness = mean(center^3)/variance^(3/2), 
            excess_kurtosis = mean(center^4)/variance^2 - 3)
    }

summarize_analysis <- function(dat, analysis, alpha) {
        suffix <- if (analysis == "unadjusted") 
            "unadjusted"
        else "adjusted"
        r_s <- rejection_summary(dat[[paste0("reject_siga_s_", 
            suffix)]])
        r_r <- rejection_summary(dat[[paste0("reject_siga_r_", 
            suffix)]])
        r_o <- if (isTRUE(INCLUDE_ORACLE_N)) {
            rejection_summary(dat[[paste0("reject_oracle_n_", 
                suffix)]])
        }
        else c(probability = NA_real_, lower = NA_real_, upper = NA_real_)
        r_rt <- rejection_summary(dat[[paste0("reject_rt_", suffix)]])
        diff_s_rt <- paired_difference_interval(dat[[paste0("reject_siga_s_", 
            suffix)]], dat[[paste0("reject_rt_", suffix)]])
        diff_r_rt <- paired_difference_interval(dat[[paste0("reject_siga_r_", 
            suffix)]], dat[[paste0("reject_rt_", suffix)]])
        diff_o_rt <- if (isTRUE(INCLUDE_ORACLE_N)) {
            paired_difference_interval(dat[[paste0("reject_oracle_n_", 
                suffix)]], dat[[paste0("reject_rt_", suffix)]])
        }
        else c(estimate = NA_real_, se = NA_real_, lower = NA_real_, 
            upper = NA_real_)
        diff_r_o <- if (isTRUE(INCLUDE_ORACLE_N)) {
            paired_difference_interval(dat[[paste0("reject_siga_r_", 
                suffix)]], dat[[paste0("reject_oracle_n_", suffix)]])
        }
        else c(estimate = NA_real_, se = NA_real_, lower = NA_real_, 
            upper = NA_real_)
        statistic <- as.numeric(dat[[paste0("statistic_", suffix)]])
        v_s <- as.numeric(dat[[paste0("variance_s_", suffix)]])
        v_r <- as.numeric(dat[[paste0("variance_r_", suffix)]])
        v_rt <- as.numeric(dat[[paste0("variance_rt_", suffix)]])
        p_s <- as.numeric(dat[[paste0("p_siga_s_", suffix)]])
        p_r <- as.numeric(dat[[paste0("p_siga_r_", suffix)]])
        p_o <- if (isTRUE(INCLUDE_ORACLE_N)) {
            as.numeric(dat[[paste0("p_oracle_n_", suffix)]])
        }
        else rep(NA_real_, nrow(dat))
        p_rt <- as.numeric(dat[[paste0("p_rt_", suffix)]])
        ratio_r_s <- safe_ratio(v_r, v_s)
        ratio_s_rt <- safe_ratio(v_s, v_rt)
        ratio_r_rt <- safe_ratio(v_r, v_rt)
        stat_moments <- central_moment_summary(statistic)
        empirical_var_stat <- unname(stat_moments["variance"])
        mean_v_s <- mean(v_s, na.rm = TRUE)
        mean_v_r <- mean(v_r, na.rm = TRUE)
        mean_v_rt <- mean(v_rt, na.rm = TRUE)
        target_ratio <- mean_v_rt/empirical_var_stat
        estimated_ratio <- mean_v_r/mean_v_s
        zcrit <- qnorm(1 - alpha)
        predicted_rt_from_target_ratio <- if (is.finite(target_ratio) && 
            target_ratio > 0) {
            1 - pnorm(zcrit * sqrt(target_ratio))
        }
        else NA_real_
        trial_prediction <- 1 - pnorm(zcrit * sqrt(ratio_r_s))
        z_s <- statistic/sqrt(v_s)
        z_r <- statistic/sqrt(v_r)
        z_s_moments <- central_moment_summary(z_s)
        z_r_moments <- central_moment_summary(z_r)
        rt_skewness <- as.numeric(dat[[paste0("randomization_skewness_", 
            suffix)]])
        rt_kurtosis <- as.numeric(dat[[paste0("randomization_excess_kurtosis_", 
            suffix)]])
        gap_source <- if (!isTRUE(INCLUDE_ORACLE_N)) {
            "oracle_not_requested"
        }
        else {
            tolerance <- max(VALIDATION_ABS_TOLERANCE, 3 * max(ifelse(is.finite(diff_r_o["se"]), 
                diff_r_o["se"], 0), ifelse(is.finite(diff_o_rt["se"]), 
                diff_o_rt["se"], 0)))
            variance_gap <- abs(diff_r_o["estimate"]) > tolerance
            tail_gap <- abs(diff_o_rt["estimate"]) > tolerance
            if (variance_gap && tail_gap) {
                "mixed_variance_and_tail"
            }
            else if (variance_gap) {
                "variance_approximation"
            }
            else if (tail_gap) {
                "non_gaussian_or_discrete_tail"
            }
            else {
                "first_order_alignment"
            }
        }
        data.frame(analysis = analysis, siga_s = unname(r_s["probability"]), 
            siga_s_lower_95 = unname(r_s["lower"]), siga_s_upper_95 = unname(r_s["upper"]), 
            siga_r = unname(r_r["probability"]), siga_r_lower_95 = unname(r_r["lower"]), 
            siga_r_upper_95 = unname(r_r["upper"]), oracle_n = unname(r_o["probability"]), 
            oracle_n_lower_95 = unname(r_o["lower"]), oracle_n_upper_95 = unname(r_o["upper"]), 
            rt = unname(r_rt["probability"]), rt_lower_95 = unname(r_rt["lower"]), 
            rt_upper_95 = unname(r_rt["upper"]), difference_siga_s_minus_rt = unname(diff_s_rt["estimate"]), 
            difference_siga_s_minus_rt_se = unname(diff_s_rt["se"]), 
            difference_siga_s_minus_rt_lower_95 = unname(diff_s_rt["lower"]), 
            difference_siga_s_minus_rt_upper_95 = unname(diff_s_rt["upper"]), 
            difference_siga_r_minus_rt = unname(diff_r_rt["estimate"]), 
            difference_siga_r_minus_rt_se = unname(diff_r_rt["se"]), 
            difference_siga_r_minus_rt_lower_95 = unname(diff_r_rt["lower"]), 
            difference_siga_r_minus_rt_upper_95 = unname(diff_r_rt["upper"]), 
            difference_oracle_n_minus_rt = unname(diff_o_rt["estimate"]), 
            difference_oracle_n_minus_rt_se = unname(diff_o_rt["se"]), 
            difference_oracle_n_minus_rt_lower_95 = unname(diff_o_rt["lower"]), 
            difference_oracle_n_minus_rt_upper_95 = unname(diff_o_rt["upper"]), 
            difference_siga_r_minus_oracle_n = unname(diff_r_o["estimate"]), 
            difference_siga_r_minus_oracle_n_se = unname(diff_r_o["se"]), 
            difference_siga_r_minus_oracle_n_lower_95 = unname(diff_r_o["lower"]), 
            difference_siga_r_minus_oracle_n_upper_95 = unname(diff_r_o["upper"]), 
            empirical_variance_statistic = empirical_var_stat, 
            mean_variance_s = mean_v_s, mean_variance_r = mean_v_r, 
            mean_variance_rt = mean_v_rt, relative_bias_variance_s = mean_v_s/empirical_var_stat - 
                1, relative_bias_variance_r_vs_rt = mean_v_r/mean_v_rt - 
                1, target_variance_ratio_rt_over_sampling = target_ratio, 
            estimated_variance_ratio_r_over_s = estimated_ratio, 
            mean_variance_ratio_r_over_s = mean(ratio_r_s, na.rm = TRUE), 
            sd_variance_ratio_r_over_s = sd(ratio_r_s, na.rm = TRUE), 
            q025_variance_ratio_r_over_s = unname(quantile(ratio_r_s, 
                0.025, na.rm = TRUE)), median_variance_ratio_r_over_s = median(ratio_r_s, 
                na.rm = TRUE), q975_variance_ratio_r_over_s = unname(quantile(ratio_r_s, 
                0.975, na.rm = TRUE)), mean_variance_ratio_s_over_rt = mean(ratio_s_rt, 
                na.rm = TRUE), mean_variance_ratio_r_over_rt = mean(ratio_r_rt, 
                na.rm = TRUE), predicted_rt_size_from_target_ratio = predicted_rt_from_target_ratio, 
            mean_trial_level_normal_prediction = mean(trial_prediction, 
                na.rm = TRUE), mean_abs_p_difference_siga_s_rt = mean(abs(p_s - 
                p_rt), na.rm = TRUE), mean_abs_p_difference_siga_r_rt = mean(abs(p_r - 
                p_rt), na.rm = TRUE), mean_abs_p_difference_oracle_n_rt = mean(abs(p_o - 
                p_rt), na.rm = TRUE), p95_abs_p_difference_siga_s_rt = unname(quantile(abs(p_s - 
                p_rt), 0.95, na.rm = TRUE)), p95_abs_p_difference_siga_r_rt = unname(quantile(abs(p_r - 
                p_rt), 0.95, na.rm = TRUE)), p95_abs_p_difference_oracle_n_rt = if (isTRUE(INCLUDE_ORACLE_N)) {
                unname(quantile(abs(p_o - p_rt), 0.95, na.rm = TRUE))
            }
            else NA_real_, statistic_mean = unname(stat_moments["mean"]), 
            statistic_skewness = unname(stat_moments["skewness"]), 
            statistic_excess_kurtosis = unname(stat_moments["excess_kurtosis"]), 
            z_s_mean = unname(z_s_moments["mean"]), z_s_variance = unname(z_s_moments["variance"]), 
            z_s_skewness = unname(z_s_moments["skewness"]), z_s_excess_kurtosis = unname(z_s_moments["excess_kurtosis"]), 
            z_r_mean = unname(z_r_moments["mean"]), z_r_variance = unname(z_r_moments["variance"]), 
            z_r_skewness = unname(z_r_moments["skewness"]), z_r_excess_kurtosis = unname(z_r_moments["excess_kurtosis"]), 
            mean_rt_randomization_skewness = mean(rt_skewness, 
                na.rm = TRUE), mean_rt_randomization_excess_kurtosis = mean(rt_kurtosis, 
                na.rm = TRUE), correlation_zs2_rho = suppressWarnings(cor(z_s^2, 
                ratio_r_s, use = "complete.obs")), mean_empty_joint_strata = mean(dat$empty_joint_strata, 
                na.rm = TRUE), probability_any_empty_joint_stratum = mean(dat$empty_joint_strata > 
                0, na.rm = TRUE), mean_min_positive_joint_stratum_count = mean(dat$min_positive_joint_stratum_count, 
                na.rm = TRUE), kappa_raw_mean = mean(dat[[paste0("kappa_raw_", 
                suffix)]], na.rm = TRUE), kappa_truncation_rate = mean(dat[[paste0("kappa_truncated_", 
                suffix)]], na.rm = TRUE), safeguard_rate = mean(dat[[paste0("safeguard_active_", 
                suffix)]], na.rm = TRUE), gap_explanation = gap_source, 
            stringsAsFactors = FALSE)
    }

summarize_scenario <- function(dat, scenario) {
        design <- DESIGNS[[as.character(scenario$design_id)]]
        model <- SCENARIO_MODELS[[as.character(scenario$scenario_id)]]
        by_analysis <- rbind(summarize_analysis(dat, "unadjusted", 
            scenario$alpha), summarize_analysis(dat, "adjusted", 
            scenario$alpha))
        fixed <- data.frame(scenario_id = scenario$scenario_id, 
            scenario_code = scenario$scenario_code, scenario_label = scenario$scenario_label, 
            scenario_class = scenario$scenario_class, set = scenario$set, 
            design_id = scenario$design_id, design_label = design$design_label, 
            factor_count = design$factor_count, target_per_group = design$target_per_group, 
            total_n = design$total_n, pbc = design$pbc, profile_type = design$profile_type, 
            outcome_type = scenario$outcome_type, model_type = scenario$model_type, 
            direction = scenario$direction, boundary = scenario$boundary, 
            alpha = scenario$alpha, achieved_effect = model$achieved_effect, 
            pair_ratio_model = model$pair_ratio, pair_scalar_gap_model = model$pair_scalar_gap, 
            n_outer = nrow(dat), rerandomizations_per_trial = N_RERANDOMIZATIONS, 
            allocation_calibration_paths = N_CALIBRATION, direction_calibration_paths = N_DIRECTION_CALIBRATION, 
            base_seed = BASE_SEED, run_version = RUN_VERSION, 
            include_oracle_n = INCLUDE_ORACLE_N, rt_shape_diagnostics = RT_SHAPE_DIAGNOSTICS, 
            machine_count = length(unique(dat$machine_id)), machine_ids = paste(sort(unique(dat$machine_id)), 
                collapse = ";"), mean_realized_pair_scalar_gap = mean(dat$pair_scalar_gap_realized, 
                na.rm = TRUE), mean_observed_treated = mean(dat$observed_treated, 
                na.rm = TRUE), data_generation_minutes = sum(dat$data_generation_seconds, 
                na.rm = TRUE)/60, score_construction_minutes = sum(dat$score_construction_seconds, 
                na.rm = TRUE)/60, siga_analysis_minutes = sum(dat$siga_analysis_seconds, 
                na.rm = TRUE)/60, rt_analysis_minutes = sum(dat$rt_analysis_seconds, 
                na.rm = TRUE)/60, stringsAsFactors = FALSE)
        cbind(fixed[rep(1L, nrow(by_analysis)), , drop = FALSE], 
            by_analysis)
    }

write_latex_tables <- function(summary, output_directory) {
        fmt <- function(x, digits = 2) {
            ifelse(is.finite(x), formatC(x, format = "f", digits = digits), 
                "--")
        }
        bs <- intToUtf8(92L)
        tex <- function(x) paste0(bs, x)
        row_end <- paste0(bs, bs)
        operating <- c(tex("begin{landscape}"), tex("begin{table}[!htbp]"), 
            tex("centering"), tex("caption{Supplemental pair-path sensitivity analysis: rejection probabilities.}"), 
            tex("label{tab:supp-pair-path-operating}"), tex("scriptsize"), 
            paste0(tex("setlength{"), tex("tabcolsep"), "}{3.2pt}"), 
            paste0(tex("renewcommand{"), tex("arraystretch"), 
                "}{1.06}"), tex("begin{tabular}{rcrllllrrrrr}"), 
            tex("toprule"), paste0("Code & $F$ & $n$/group & Outcome & Model & Analysis & Class & ", 
                "SIGA-S & SIGA-R & Oracle-N & RT & R$-$RT ", 
                row_end), tex("midrule"))
        for (i in seq_len(nrow(summary))) {
            row <- summary[i, ]
            operating <- c(operating, paste0(row$scenario_code, 
                " & ", row$factor_count, " & ", row$target_per_group, 
                " & ", row$outcome_type, " & ", gsub("_", "-", 
                  row$model_type, fixed = TRUE), " & ", row$analysis, 
                " & ", gsub("_", "-", row$scenario_class, fixed = TRUE), 
                " & ", fmt(100 * row$siga_s, 2), " & ", fmt(100 * 
                  row$siga_r, 2), " & ", fmt(100 * row$oracle_n, 
                  2), " & ", fmt(100 * row$rt, 2), " & ", fmt(100 * 
                  row$difference_siga_r_minus_rt, 2), " ", row_end))
        }
        operating <- c(operating, tex("bottomrule"), tex("end{tabular}"), 
            paste0(tex("par"), tex("vspace{2pt}"), tex("footnotesize "), 
                "Rejection probabilities and R$-$RT are percentages and percentage points, respectively. ", 
                "Oracle-N is the Gaussian tail probability using the Monte Carlo RT conditional mean and variance ", 
                "already calculated from the same regenerated paths; it requires no additional rerandomization."), 
            tex("end{table}"), tex("end{landscape}"))
        diagnostic <- c(tex("begin{landscape}"), tex("begin{table}[!htbp]"), 
            tex("centering"), tex("caption{Supplemental diagnostics for the pair-path sensitivity analysis.}"), 
            tex("label{tab:supp-pair-path-diagnostics}"), tex("scriptsize"), 
            paste0(tex("setlength{"), tex("tabcolsep"), "}{2.8pt}"), 
            paste0(tex("renewcommand{"), tex("arraystretch"), 
                "}{1.06}"), tex("begin{tabular}{rlrrrrrrrrl}"), 
            tex("toprule"), paste0("Code & Analysis & Mean $", 
                tex("widehat{"), tex("rho"), "_n}$ & Mean $", 
                tex("widehat{V}_R/V_{RT}"), "$ & RT skew. & RT ex. kurt. & $P$(empty) & ", 
                "Mean empty & R$-$Oracle & Oracle$-$RT & Diagnostic ", 
                row_end), tex("midrule"))
        for (i in seq_len(nrow(summary))) {
            row <- summary[i, ]
            diagnostic <- c(diagnostic, paste0(row$scenario_code, 
                " & ", row$analysis, " & ", fmt(row$mean_variance_ratio_r_over_s, 
                  3), " & ", fmt(row$mean_variance_ratio_r_over_rt, 
                  3), " & ", fmt(row$mean_rt_randomization_skewness, 
                  3), " & ", fmt(row$mean_rt_randomization_excess_kurtosis, 
                  3), " & ", fmt(100 * row$probability_any_empty_joint_stratum, 
                  1), " & ", fmt(row$mean_empty_joint_strata, 
                  2), " & ", fmt(100 * row$difference_siga_r_minus_oracle_n, 
                  2), " & ", fmt(100 * row$difference_oracle_n_minus_rt, 
                  2), " & ", gsub("_", "-", row$gap_explanation, 
                  fixed = TRUE), " ", row_end))
        }
        diagnostic <- c(diagnostic, tex("bottomrule"), tex("end{tabular}"), 
            paste0(tex("par"), tex("vspace{2pt}"), tex("footnotesize "), 
                "$V_{RT}$ denotes the Monte Carlo conditional variance of the regenerated fixed-score statistics. ", 
                "$P$(empty) is the percentage of outer trials with at least one empty joint stratum. ", 
                "R$-$Oracle and Oracle$-$RT are paired rejection-probability differences in percentage points. ", 
                "The final column distinguishes variance-approximation error from non-Gaussian or discrete-tail error."), 
            tex("end{table}"), tex("end{landscape}"))
        operating_path <- file.path(output_directory, "pair_path_supplemental_operating_table.tex")
        diagnostic_path <- file.path(output_directory, "pair_path_supplemental_diagnostic_table.tex")
        combined_path <- file.path(output_directory, "pair_path_supplemental_tables.tex")
        writeLines(operating, operating_path)
        writeLines(diagnostic, diagnostic_path)
        writeLines(c(operating, "", diagnostic), combined_path)
        invisible(c(operating_path, diagnostic_path, combined_path))
    }

apply_validation_gates <- function(summary) {
        rows <- lapply(seq_len(nrow(summary)), function(i) {
            row <- summary[i, ]
            se_candidates <- c(row$difference_siga_r_minus_rt_se, 
                row$difference_siga_r_minus_oracle_n_se, row$difference_oracle_n_minus_rt_se)
            se_candidates <- se_candidates[is.finite(se_candidates)]
            paired_tolerance <- max(VALIDATION_ABS_TOLERANCE, 
                if (length(se_candidates)) 3 * max(se_candidates) else 0)
            r_close_rt <- abs(row$difference_siga_r_minus_rt) <= 
                paired_tolerance
            r_close_oracle <- if (isTRUE(INCLUDE_ORACLE_N)) {
                abs(row$difference_siga_r_minus_oracle_n) <= 
                  paired_tolerance
            }
            else NA
            oracle_close_rt <- if (isTRUE(INCLUDE_ORACLE_N)) {
                abs(row$difference_oracle_n_minus_rt) <= paired_tolerance
            }
            else NA
            rt_agreement_or_explained <- if (isTRUE(INCLUDE_ORACLE_N)) {
                r_close_rt || r_close_oracle
            }
            else r_close_rt
            s_nominal <- abs(row$siga_s - row$alpha) <= NOMINAL_ABS_TOLERANCE
            variance_r_close <- abs(row$relative_bias_variance_r_vs_rt) <= 
                VARIANCE_REL_TOLERANCE
            variance_s_close <- abs(row$relative_bias_variance_s) <= 
                VARIANCE_REL_TOLERANCE
            safeguard_ok <- row$safeguard_rate <= 0.001
            direction_ok <- TRUE
            if (row$scenario_class %in% c("strong_pair_path", 
                "correlated_factor_stress") && row$pair_ratio_model > 
                1.05) {
                direction_ok <- row$siga_r <= row$siga_s && row$rt <= 
                  row$siga_s
            }
            primary_pass <- if (row$set == "supplemental") {
                s_nominal && variance_r_close && variance_s_close && 
                  safeguard_ok && direction_ok && rt_agreement_or_explained
            }
            else NA
            data.frame(scenario_id = row$scenario_id, scenario_code = row$scenario_code, 
                scenario_class = row$scenario_class, analysis = row$analysis, 
                paired_tolerance = paired_tolerance, r_close_rt = r_close_rt, 
                r_close_oracle = r_close_oracle, oracle_close_rt = oracle_close_rt, 
                rt_agreement_or_explained = rt_agreement_or_explained, 
                s_nominal = s_nominal, variance_r_close_to_rt = variance_r_close, 
                variance_s_close_to_sampling = variance_s_close, 
                safeguard_ok = safeguard_ok, direction_ok = direction_ok, 
                gap_explanation = row$gap_explanation, primary_pass = primary_pass, 
                stringsAsFactors = FALSE)
        })
        do.call(rbind, rows)
    }

aggregate_results <- function(validate = FALSE) {
        all_summary <- list()
        all_details <- list()
        for (i in seq_len(nrow(SCENARIOS))) {
            scenario <- SCENARIOS[i, ]
            message("Aggregating ", scenario$scenario_label)
            dat <- read_scenario_details(scenario)
            all_details[[i]] <- dat
            all_summary[[i]] <- summarize_scenario(dat, scenario)
        }
        summary <- do.call(rbind, all_summary)
        details <- do.call(rbind, all_details)
        summary <- summary[order(summary$scenario_id, summary$analysis), 
            , drop = FALSE]
        details <- details[order(details$scenario_id, details$replicate), 
            , drop = FALSE]
        validation <- apply_validation_gates(summary)
        safe_write_csv(summary, file.path(OUTPUT_DIR, "pair_path_supplemental_summary.csv"))
        safe_write_csv(details, file.path(OUTPUT_DIR, "pair_path_supplemental_details.csv"))
        safe_write_csv(validation, file.path(OUTPUT_DIR, "pair_path_supplemental_gates.csv"))
        write_latex_tables(summary, OUTPUT_DIR)
        message("Summary written to: ", file.path(OUTPUT_DIR, 
            "pair_path_supplemental_summary.csv"))
        if (validate) {
            failed <- validation$primary_pass %in% FALSE
            if (any(failed)) {
                failed_labels <- paste(paste0(validation$scenario_code[failed], 
                  "/", validation$analysis[failed]), collapse = ", ")
                message("Supplemental validation failures: ", 
                  failed_labels)
                if (STRICT_VALIDATION) {
                  stop("One or more supplemental validation gates failed.", 
                    call. = FALSE)
                }
            }
            else {
                message("All supplemental validation gates passed.")
            }
        }
        invisible(list(summary = summary, validation = validation))
    }

message("Final audit-first SIGA pair-path supplemental sensitivity analysis")

message("Profile: ", PROFILE)

message("Mode: ", RUN_MODE)

message("Scenario set: ", SCENARIO_SET)

message("Run version: ", 
        RUN_VERSION)

message("Machine ID: ", MACHINE_ID)

message("Oracle-N included: ", 
        INCLUDE_ORACLE_N)

message("RT shape diagnostics included: ", 
        RT_SHAPE_DIAGNOSTICS)

message("Outer trials per scenario: ", 
        N_OUTER)

message("Rerandomizations per trial: ", N_RERANDOMIZATIONS)

message("Analysis allocation calibration paths: ", N_CALIBRATION)

message("Independent direction-selection paths: ", N_DIRECTION_CALIBRATION)

message("SIGA-R safeguard exponent: ", EPSILON_EXPONENT)

message("Selected scenarios: ", paste(SCENARIOS$scenario_id, 
        collapse = ","))

message("Output directory: ", OUTPUT_DIR)

if (RUN_MODE == "run") {
        for (i in seq_len(nrow(SCENARIOS))) run_scenario_shard(SCENARIOS[i, 
            ])
    } else if (RUN_MODE == "aggregate") {
        aggregate_results(validate = FALSE)
    } else if (RUN_MODE == "validate") {
        aggregate_results(validate = TRUE)
    }
