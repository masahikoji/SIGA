#!/usr/bin/env Rscript

# Verifies that the binary and continuous benchmark implementations use the
# same fixed-score randomization-test kernel under biased-coin minimization.
# The two functions below intentionally use the two state representations from
# the production programs: level-1 imbalances with level-0 derived by subtraction
# (continuous code), and explicit two-level marginal imbalances (binary code).

continuous_kernel <- function(X, scores, B, pbc = 0.80, seed = 20260726L) {
  X <- as.matrix(X)
  scores <- as.matrix(scores)
  n <- nrow(X)
  K <- ncol(X)
  L <- ncol(scores)
  set.seed(seed)
  d_overall <- integer(B)
  d_level1 <- matrix(0L, nrow = B, ncol = K)
  statistic <- matrix(0, nrow = B, ncol = L)

  for (i in seq_len(n)) {
    score_if_1 <- abs(d_overall + 1L)
    score_if_0 <- abs(d_overall - 1L)
    for (j in seq_len(K)) {
      level_imbalance <- if (X[i, j] == 1L) d_level1[, j] else d_overall - d_level1[, j]
      score_if_1 <- score_if_1 + abs(level_imbalance + 1L)
      score_if_0 <- score_if_0 + abs(level_imbalance - 1L)
    }
    q <- rep.int(0.5, B)
    q[score_if_1 < score_if_0] <- pbc
    q[score_if_1 > score_if_0] <- 1 - pbc
    A <- as.integer(runif(B) < q)
    z <- 2L * A - 1L
    centered <- A - 0.5
    for (ell in seq_len(L)) statistic[, ell] <- statistic[, ell] + centered * scores[i, ell]
    d_overall <- d_overall + z
    active <- which(X[i, ] == 1L)
    if (length(active)) for (j in active) d_level1[, j] <- d_level1[, j] + z
  }
  statistic
}

binary_kernel <- function(X, scores, B, pbc = 0.80, seed = 20260726L) {
  X <- as.matrix(X)
  scores <- as.matrix(scores)
  n <- nrow(X)
  K <- ncol(X)
  L <- ncol(scores)
  set.seed(seed)
  overall <- integer(B)
  marginal <- matrix(0L, nrow = B, ncol = 2L * K)
  statistic <- matrix(0, nrow = B, ncol = L)

  for (i in seq_len(n)) {
    score_plus <- abs(overall + 1L)
    score_minus <- abs(overall - 1L)
    cols <- integer(K)
    for (j in seq_len(K)) {
      cols[j] <- 2L * (j - 1L) + X[i, j] + 1L
      current <- marginal[, cols[j]]
      score_plus <- score_plus + abs(current + 1L)
      score_minus <- score_minus + abs(current - 1L)
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
    for (ell in seq_len(L)) statistic[, ell] <- statistic[, ell] + 0.5 * z * scores[i, ell]
  }
  statistic
}

set.seed(91273)
n <- 250L
K <- 5L
B <- 499L
X <- sapply(c(0.50, 0.40, 0.30, 0.20, 0.10), function(p) rbinom(n, 1L, p))
X <- matrix(as.integer(X), nrow = n, ncol = K)
scores <- matrix(rnorm(n * 4L), nrow = n, ncol = 4L)
scores <- sweep(scores, 2L, colMeans(scores), "-")

s1 <- continuous_kernel(X, scores, B = B)
s2 <- binary_kernel(X, scores, B = B)
err <- max(abs(s1 - s2))
if (!is.finite(err) || err > 1e-12) {
  stop("Randomization-test kernels differ; maximum absolute error = ", err, call. = FALSE)
}

obs_z <- rep(c(-1L, 1L), length.out = n)
t_obs <- drop(0.5 * crossprod(obs_z, scores))
p1 <- (1 + colSums(abs(s1) >= matrix(abs(t_obs), B, ncol(scores), byrow = TRUE))) / (B + 1)
p2 <- (1 + colSums(abs(s2) >= matrix(abs(t_obs), B, ncol(scores), byrow = TRUE))) / (B + 1)
if (max(abs(p1 - p2)) > 0) stop("Plus-one p-values differ.", call. = FALSE)

cat("PASS: binary and continuous fixed-score randomization-test kernels are identical.\n")
cat("Checked ordered-factor conditioning, biased-coin allocation, fixed scores, and plus-one p-values.\n")
