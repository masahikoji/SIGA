#!/usr/bin/env Rscript
# Reduced internal checks for the R3 additions.  Runs in a few seconds.
options(stringsAsFactors = FALSE, warn = 2)

script_directory <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE)))
  normalizePath(getwd(), mustWork = FALSE)
}
BASE_DIR <- script_directory()
source(file.path(BASE_DIR, "siga_pair_path_engine.R"), local = FALSE)
source(file.path(BASE_DIR, "r3_engine_additions.R"), local = FALSE)

set.seed(1)
K <- 2L; n <- 60L
patterns <- all_binary_patterns(K)
pi0 <- profile_prob_independent(c(0.5, 0.4), patterns)
cal <- calibrate_pair_path_design_R(B0 = 400L, n = n, patterns = patterns, profile_prob = pi0,
                                    pbc = 0.8, weights_ = rep(1, K + 1L), seed = 7, batch_size = 200L, progress = FALSE)
cal$directions <- generalized_pair_directions(cal$psi, cal$profile_prob)

# 1. Heterogeneous models satisfy pi'd0 = 0 and the marginal effect
v <- make_direction("first_factor", patterns, pi0)
d0 <- scale_centred_direction(v, pi0, 0.5)
stopifnot(abs(sum(pi0 * d0)) < 1e-12, abs(max(abs(d0)) - 0.5) < 1e-12)
mc <- make_continuous_model_r3(patterns, pi0, true_effect = -0.2, d0 = d0)
stopifnot(abs(mc$achieved_effect + 0.2) < 1e-12)
mb <- make_binary_direct_model_r3(patterns, pi0, true_effect = -0.1, d0 = scale_centred_direction(v, pi0, 0.15))
stopifnot(abs(mb$achieved_effect + 0.1) < 1e-12, all(mb$p1 > 0 & mb$p1 < 1))
stopifnot(max(abs(model_based_deviation(mc, -0.2) - d0)) < 1e-12)
vi <- make_direction("interaction", patterns, pi0)
stopifnot(abs(sum(pi0 * vi)) < 1e-12)

# 2. Data-based deviations agree with a direct computation
prof <- generate_profile_sequence(n, patterns, pi0)
z <- ps_assign_R(prof$X, pbc = 0.8, seed = 11); A <- as.integer((z + 1L) / 2L)
y <- generate_outcome_from_model(prof$id, A, z, mc)
dh <- estimate_stratum_effect_deviations(y, A, prof$id, nrow(patterns), -0.2)
for (s in seq_len(nrow(patterns))) {
  i1 <- prof$id == s & A == 1L; i0 <- prof$id == s & A == 0L
  if (any(i1) && any(i0)) {
    stopifnot(abs(dh$d[s] - (mean(y[i1]) - mean(y[i0]) + 0.2)) < 1e-12)
    if (sum(i1) > 1 && sum(i0) > 1) stopifnot(abs(dh$var_d[s] - (var(y[i1]) / sum(i1) + var(y[i0]) / sum(i0))) < 1e-12)
  } else stopifnot(dh$d[s] == 0)
}

# 3. Corrected RT with lambda = 1 reproduces the reference test on identical paths
scores <- make_boundary_scores(y, A, prof$X, -0.2)
S <- cbind(unadjusted = scores$unadjusted, adjusted = scores$adjusted)
ref <- rt_pvalues_R(prof$X, S, z, B = 199L, pbc = 0.8, seed = 99)
lam <- cbind(rt = c(1, 1), crt = c(1.2, 0.8))
cor <- rt_corrected_pvalues_R(prof$X, S, z, lam, B = 199L, pbc = 0.8, seed = 99)
stopifnot(max(abs(cor$greater[, "rt"] - ref$greater)) < 1e-12,
          max(abs(cor$two_sided[, "rt"] - ref$two_sided)) < 1e-12,
          max(abs(cor$less[, "rt"] - ref$less)) < 1e-12,
          max(abs(cor$randomization_variance - ref$randomization_variance)) < 1e-10)
# lambda > 1 with positive statistic cannot increase the upper-tail p-value's complement
stopifnot(all(cor$two_sided[, "crt"] >= 0 & cor$two_sided[, "crt"] <= 1))

# 3b. Tie convention: with lambda > 1 the corrected p-value never falls below the count of
#     regenerated statistics tied with the observed one; on a lattice-valued score with
#     lambda - 1 below the lattice spacing it coincides with the reference test.
yb <- as.numeric(runif(n) < 0.5); sb <- make_boundary_scores(yb, A, prof$X, 0)
Sb <- cbind(unadjusted = sb$unadjusted)
refb <- rt_pvalues_R(prof$X, Sb, z, B = 399L, pbc = 0.8, seed = 5)
corb <- rt_corrected_pvalues_R(prof$X, Sb, z, cbind(rt = 1, crt = 1.0005), B = 399L, pbc = 0.8, seed = 5)
stopifnot(abs(corb$two_sided[1, "crt"] - refb$two_sided[1]) < 1e-12, abs(corb$greater[1, "crt"] - refb$greater[1]) < 1e-12)

# 4. corrected_variance: d = 0 gives lambda = 1; floor behaves
samp <- siga_sampling_variance(S[, 1L], prof$X, cal)
cv0 <- corrected_variance(samp, rep(0, nrow(patterns)), calibration = cal)
stopifnot(abs(cv0$lambda - 1) < 1e-12, !cv0$floored)
cvm <- corrected_variance(samp, model_based_deviation(mc, -0.2), calibration = cal)
stopifnot(is.finite(cvm$lambda), cvm$lambda > 0)
cvd <- corrected_variance(samp, dh$d, var_d = dh$var_d, calibration = cal, debias = TRUE)
stopifnot(is.finite(cvd$lambda), cvd$variance >= cvd$variance_raw - 1e-12)

# 5. Lattice mixture p-value equals one at an extreme statistic and lies in [0,1]
pl <- lattice_mixture_pvalue(-1e6, samp$variance, n, as.integer(sum(y > 0)), cal$treated_count_prob, "greater")
stopifnot(abs(pl - 1) < 1e-9)
pl2 <- lattice_mixture_pvalue(2, samp$variance, n, 30L, cal$treated_count_prob, "two.sided")
stopifnot(pl2 >= 0, pl2 <= 1)

# 6. Compiled kernel (if a C compiler is available) is bit-identical to the pure-R kernel
if (load_rt_kernel(BASE_DIR, quiet = TRUE)) {
  cc <- rt_corrected_pvalues_C(prof$X, S, z, lam, B = 199L, pbc = 0.8, seed = 99)
  stopifnot(identical(cc$greater, cor$greater), identical(cc$two_sided, cor$two_sided), identical(cc$less, cor$less),
            isTRUE(all.equal(cc$randomization_variance, cor$randomization_variance)))
  cat("Compiled kernel: available and identical to the pure-R kernel.\n")
} else {
  cat("Compiled kernel: not available (no C compiler); the pure-R kernel will be used.\n")
}

cat("PASS: R3 engine checks completed.\n")
