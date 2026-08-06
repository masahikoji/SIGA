#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
base <- if (length(args)) args[1L] else file.path("data", "production_results", "raw", "pair_path_supplemental")
summary <- read.csv(file.path(base, "pair_path_supplemental_summary.csv"), check.names = FALSE)
gates <- read.csv(file.path(base, "pair_path_supplemental_gates.csv"), check.names = FALSE)
stopifnot(
  nrow(summary) == 24L,
  length(unique(summary$scenario_id)) == 12L,
  all(summary$n_outer == 100000L),
  all(summary$rerandomizations_per_trial == 4999L),
  all(summary$allocation_calibration_paths == 100000L),
  all(summary$direction_calibration_paths == 200000L),
  all(summary$base_seed == 20260805L),
  all(summary$run_version == "v20260805_supplemental_v2"),
  all(summary$safeguard_rate == 0),
  nrow(gates) == 24L,
  all(gates$primary_pass)
)
assert_close <- function(actual, expected, tolerance, label) {
  if (!is.finite(actual) || abs(actual - expected) > tolerance) {
    stop(label, ": actual=", actual, ", expected=", expected, call. = FALSE)
  }
}
assert_close(100 * max(abs(summary$difference_siga_s_minus_rt)), 0.900, 1e-12, "max |SIGA-S - RT| (pp)")
assert_close(100 * max(abs(summary$difference_siga_r_minus_rt)), 0.244, 1e-12, "max |SIGA-R - RT| (pp)")
assert_close(max(summary$mean_variance_ratio_r_over_s), 1.19359594637313, 1e-13, "max mean variance ratio")
cat("PASS: 12 scenarios, 24 analyses, all primary gates passed; maxima are 0.900 and 0.244 pp.\n")
