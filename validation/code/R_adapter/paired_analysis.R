# Data-only paired comparison using the user's existing R allocation engine.
# Source siga_pair_path_engine.R, r3_engine_additions.R and
# covariance_refinement.R before this file.

corrected_lattice_tail <- function(T, VR, VS, n, Yplus, weights,
                                   corrected = FALSE, tolerance = 1e-12) {
  if (!is.finite(VR) || !is.finite(VS) || min(VR, VS) <= 0) return(1)
  lambda <- if (corrected) sqrt(VR / VS) else 1
  ans <- 0
  for (k in which(weights > 0) - 1L) {
    x <- seq.int(max(0L, k - (n - Yplus)), min(k, Yplus))
    t <- x - k * Yplus / n
    a <- (t - .5) / sqrt(VR); b <- (t + .5) / sqrt(VR)
    masses <- ifelse(a >= 0,
                      pnorm(a, lower.tail = FALSE) - pnorm(b, lower.tail = FALSE),
                      pnorm(b) - pnorm(a))
    if (sum(masses) <= 0) stop("Zero normal mass on feasible lattice.")
    use <- abs(t) + tolerance >= lambda * abs(T)
    if (corrected) use <- use | abs(abs(t) - abs(T)) <= tolerance
    ans <- ans + weights[k+1L] * sum(masses[use]) / sum(masses)
  }
  min(max(ans, 0), 1)
}

compare_observed_trial <- function(y, A, X, boundaries, tails, calibration,
                                  alpha = .05, B = 4999L, reference_seed = 1,
                                  kernel = "auto") {
  # Inputs do not include the DGP, true effect vector, or true covariance.
  if (length(boundaries) != length(tails) || !length(boundaries)) stop("Boundary/tail mismatch.")
  if (any(!tails %in% c("two.sided", "greater", "less"))) stop("Invalid tails.")
  if (is.null(calibration$xi)) calibration <- add_balancing_covariance(calibration)
  scores <- list(); details <- list(); lambdas <- list(); row_tails <- character()
  id <- joint_stratum_id(X)
  for (j in seq_along(boundaries)) {
    b <- boundaries[j]
    residuals <- make_boundary_scores(y, A, X, b)
    d <- estimate_stratum_effect_deviations(y, A, id, calibration$J, b)$d
    for (score_name in c("unadjusted", "adjusted")) {
      r <- residuals[[score_name]]
      old <- siga_sampling_variance(r, X, calibration)
      new <- siga_sampling_variance_refined(r, X, calibration)
      ro <- corrected_variance(old, d = d, calibration = calibration)
      rn <- corrected_variance(new, d = d, calibration = calibration)
      scores[[length(scores)+1L]] <- r
      details[[length(details)+1L]] <- list(b = b, score_name = score_name,
          S = c(old$variance, new$variance), R = c(ro$variance, rn$variance),
          rank_fallback = new$rank_fallback)
      lambdas[[length(lambdas)+1L]] <- c(1, ro$lambda, rn$lambda)
      row_tails <- c(row_tails, tails[j])
    }
  }
  scores <- do.call(cbind, scores); lambda <- do.call(rbind, lambdas)
  colnames(lambda) <- c("RT", "original_CRT_data", "refined_CRT_data")
  rr <- rt_corrected_pvalues(X = X, scores = scores, z_obs = 2L*A-1L,
       lambda = lambda, B = B, pbc = calibration$pbc, weights_ = calibration$weights,
       seed = reference_seed, kernel = kernel)
  methods <- c("RT", "original_CRT_data", "refined_CRT_data",
               "original_S_normal", "refined_S_normal", "original_R_normal", "refined_R_normal",
               "original_R_lattice", "refined_R_lattice", "original_CRT_lattice", "refined_CRT_lattice")
  allp <- matrix(NA_real_, length(details), length(methods), dimnames = list(NULL, methods))
  for (i in seq_along(details)) {
    d <- details[[i]]; t <- rr$statistic[i]; tail <- row_tails[i]
    tail_key <- switch(tail, two.sided = "two_sided", greater = "greater", less = "less")
    allp[i, 1:3] <- rr[[tail_key]][i, ]
    np <- function(v) {
      z <- t / sqrt(v)
      switch(tail, two.sided = 2*pnorm(-abs(z)), greater = pnorm(z, lower.tail = FALSE), less = pnorm(z))
    }
    allp[i, 4:5] <- vapply(d$S, np, numeric(1))
    allp[i, 6:7] <- vapply(d$R, np, numeric(1))
    allp[i, 8:9] <- allp[i, 6:7]
    allp[i, 10:11] <- allp[i, 4:5]
    if (all(y %in% c(0,1)) && d$b == 0 && d$score_name == "unadjusted" && tail == "two.sided") {
      for (v in 1:2) {
        allp[i, 7+v] <- corrected_lattice_tail(t, d$R[v], d$S[v], length(y), sum(y), calibration$treated_count_prob)
        allp[i, 9+v] <- corrected_lattice_tail(t, d$R[v], d$S[v], length(y), sum(y), calibration$treated_count_prob, corrected = TRUE)
      }
    }
  }
  # One boundary for superiority/NI; max of both boundaries for equivalence.
  final <- rbind(unadjusted = apply(allp[seq(1, nrow(allp), by=2), , drop=FALSE], 2, max),
                 adjusted = apply(allp[seq(2, nrow(allp), by=2), , drop=FALSE], 2, max))
  list(pvalues = final, reject = final <= alpha, boundary_pvalues = allp,
       variance_details = details, reference_variance = rr$randomization_variance)
}
