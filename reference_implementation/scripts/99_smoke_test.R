#!/usr/bin/env Rscript
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])))), "R", "load_all.R"))

patterns <- all_binary_patterns(2L)
stopifnot(identical(joint_stratum_id(patterns), 1:4))
probability <- profile_prob_independent(c(0.5, 0.4), patterns)
stopifnot(abs(sum(probability) - 1) < 1e-12)
stopifnot(max(abs(profile_marginals(patterns, probability) - c(0.5, 0.4))) < 1e-12)
correlated <- profile_prob_latent_correlated(c(0.5, 0.4), 1, 0.5, patterns)
stopifnot(abs(sum(correlated) - 1) < 1e-12)
stopifnot(max(abs(profile_marginals(patterns, correlated) - c(0.5, 0.4))) < 1e-10)

one <- calibrate_one_path_design(
  B0 = 100L, n = 40L, patterns = patterns, profile_prob = probability,
  pbc = 0.8, seed = 1234, batch_size = 50L, progress = FALSE
)
three <- calibrate_three_path_design(
  B0 = 100L, n = 40L, patterns = patterns, profile_prob = probability,
  pbc = 0.8, seed = 5678, batch_size = 50L, progress = FALSE
)
stopifnot(all(dim(one$gamma) == c(4L, 4L)))
stopifnot(all(dim(three$gamma) == c(4L, 4L)))
stopifnot(all(dim(three$psi) == c(4L, 4L)))
stopifnot(min(eigen(one$gamma, symmetric = TRUE, only.values = TRUE)$values) > -1e-8)
stopifnot(min(eigen(three$psi, symmetric = TRUE, only.values = TRUE)$values) > -1e-8)

design <- make_independent_design(2L, 20L)
model <- make_continuous_benchmark_model(design, 0.2)
trial <- generate_benchmark_trial(design, model, 11, 12, 13)
scores <- make_boundary_scores(trial$y, trial$A, trial$X, 0)
vs <- siga_sampling_variance(scores$unadjusted, trial$X, three)
vr <- siga_randomization_variance(vs, model$d_at(0), three, 1 / design$total_n)
stopifnot(is.finite(vs$variance), vs$variance > 0)
stopifnot(is.finite(vr$variance), vr$variance > 0)
rt <- fixed_score_randomization_test(
  trial$X, cbind(scores$unadjusted, scores$adjusted), trial$z,
  B = 19L, pbc = design$pbc, weights = design$weights, seed = 14
)
stopifnot(all(rt$greater >= 0 & rt$greater <= 1))
stopifnot(all(rt$two_sided >= 0 & rt$two_sided <= 1))
stopifnot(nrow(make_full_grid_scenarios()) == 56L)
stopifnot(nrow(make_stress_scenarios()) == 20L)

cat("SIGA smoke test passed.\n")
