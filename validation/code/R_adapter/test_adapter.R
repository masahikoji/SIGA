# Run from this directory: Rscript test_adapter.R
source("covariance_refinement.R")
set.seed(901)
for (F in c(2L, 5L)) {
  patterns <- as.matrix(expand.grid(rep(list(0:1), F)))
  L <- t(cbind(1, 2 * patterns - 1)); J <- nrow(patterns)
  x <- matrix(rnorm(J * J), J); gamma <- tcrossprod(x) / J
  x <- matrix(rnorm((F+1)^2), F+1); xi <- tcrossprod(x)
  for (empty in c(FALSE, TRUE)) {
    counts <- sample(3:30, J, replace = TRUE)
    if (empty) counts[1] <- 0
    z <- refined_imbalance_covariance(counts, gamma, xi, L)
    stopifnot(!z$fallback, z$balance_error < 1e-8,
              min(eigen(z$Omega, symmetric = TRUE, only.values = TRUE)$values) > -1e-8)
  }
  counts <- c(100, rep(0, J-1))
  z <- refined_imbalance_covariance(counts, gamma, xi, L)
  stopifnot(z$fallback, identical(z$Omega, z$original))
}
cat("R adapter algebra checks passed.\n")
