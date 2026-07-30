#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 2)

script_directory <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

BASE_DIR <- script_directory()
source(file.path(BASE_DIR, "siga_pair_path_engine.R"), local = FALSE)

patterns <- all_binary_patterns(2L)
stopifnot(identical(joint_stratum_id(patterns), seq_len(4L)))

pi_ind <- profile_prob_independent(c(0.5, 0.4), patterns)
stopifnot(abs(sum(pi_ind) - 1) < 1e-12)
stopifnot(max(abs(profile_marginals(patterns, pi_ind) - c(0.5, 0.4))) < 1e-12)

pi_cor <- profile_prob_latent_correlated(c(0.5, 0.4), 0.8, patterns = patterns)
stopifnot(abs(sum(pi_cor) - 1) < 1e-12)
stopifnot(max(abs(profile_marginals(patterns, pi_cor) - c(0.5, 0.4))) < 1e-10)

calibration <- calibrate_pair_path_design_R(
  B0 = 300L,
  n = 40L,
  patterns = patterns,
  profile_prob = pi_ind,
  pbc = 0.8,
  weights_ = rep(1, 3L),
  seed = 12345,
  batch_size = 100L,
  progress = FALSE
)
stopifnot(all(dim(calibration$gamma) == c(4L, 4L)))
stopifnot(all(dim(calibration$psi) == c(4L, 4L)))
stopifnot(min(eigen(calibration$gamma, symmetric = TRUE, only.values = TRUE)$values) > -1e-9)
stopifnot(min(eigen(calibration$psi, symmetric = TRUE, only.values = TRUE)$values) > -1e-9)

directions <- generalized_pair_directions(calibration$psi, pi_ind)
stopifnot(abs(sum(pi_ind * directions$min_direction)) < 1e-9)
stopifnot(abs(sum(pi_ind * directions$max_direction)) < 1e-9)
stopifnot(directions$min_ratio <= directions$max_ratio + 1e-10)

model <- make_continuous_pure_model(
  patterns,
  pi_ind,
  boundary = -0.2,
  direction = directions$max_direction,
  target_max_abs_d = 0.5,
  shared_noise_sd = 0.05
)
stopifnot(abs(sum(pi_ind * model$d)) < 1e-9)

set.seed(10)
profile <- generate_profile_sequence(40L, patterns, pi_ind)
z <- ps_assign_R(profile$X, pbc = 0.8, weights_ = rep(1, 3L), seed = 11)
A <- as.integer((z + 1L) / 2L)
set.seed(12)
y <- generate_outcome_from_model(profile$id, A, z, model)
scores <- make_boundary_scores(y, A, profile$X, boundary = -0.2)
vs <- siga_sampling_variance(scores$unadjusted, profile$X, calibration)
vr <- siga_randomization_variance(vs, model$d, calibration)
stopifnot(is.finite(vs$variance), vs$variance > 0)
stopifnot(is.finite(vr$variance), vr$variance > 0)

rt <- rt_pvalues_R(
  X = profile$X,
  scores = cbind(scores$unadjusted, scores$adjusted),
  z_obs = z,
  B = 49L,
  pbc = 0.8,
  weights_ = rep(1, 3L),
  seed = 13
)
stopifnot(length(rt$greater) == 2L)
stopifnot(all(rt$greater >= 0 & rt$greater <= 1))
stopifnot(all(is.finite(rt$randomization_variance)))

message("PASS: pair-path engine checks completed.")
