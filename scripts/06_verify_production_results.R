#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/")
root <- normalizePath(file.path(dirname(script), ".."), winslash = "/")
source(file.path(root, "R", "production_result_tools.R"))

s <- siga_build_full_grid_s(root)
stopifnot(nrow(s) == 112L, all(s$n_outer == 100000L),
          all(s$rerandomizations_per_trial == 4999L),
          all(s$calibration_paths == 100000L))
siga_assert_close(max(abs(s$difference_pp)), 0.422, label = "SIGA-S maximum absolute difference (pp)")
siga_assert_close(max(abs(s$difference_pp[s$scenario_role != "power"])), 0.185,
                  label = "SIGA-S null-boundary maximum absolute difference (pp)")
siga_assert_close(max(abs(s$difference_pp[s$scenario_role == "power"])), 0.422,
                  label = "SIGA-S power maximum absolute difference (pp)")

r <- siga_build_full_grid_r(root)
stopifnot(nrow(r) == 112L, length(unique(r$scenario_id)) == 56L,
          all(r$completed_outer_trials == 100000L))
siga_assert_close(max(abs(r$difference_pp)), 0.322, label = "SIGA-R maximum absolute difference (pp)")
siga_assert_close(max(abs(r$difference_pp[r$scenario_role != "power"])), 0.205,
                  label = "SIGA-R null-boundary maximum absolute difference (pp)")
siga_assert_close(max(abs(r$difference_pp[r$scenario_role == "power"])), 0.322,
                  label = "SIGA-R power maximum absolute difference (pp)")
ratios <- c(r$mean_ratio_r_over_s_1, r$mean_ratio_r_over_s_2)
ratios <- ratios[is.finite(ratios)]
siga_assert_close(min(ratios), 1.0, label = "minimum mean variance ratio")
siga_assert_close(max(ratios), 1.04157041225235, tolerance = 1e-13,
                  label = "maximum mean variance ratio")
stopifnot(max(c(r$safeguard_rate_1, r$safeguard_rate_2), na.rm = TRUE) == 0)

pair <- siga_build_pair_path_supplemental(root)
stopifnot(
  nrow(pair) == 24L,
  length(unique(pair$scenario_id)) == 12L,
  all(pair$n_outer == 100000L),
  all(pair$rerandomizations_per_trial == 4999L),
  all(pair$allocation_calibration_paths == 100000L),
  all(pair$direction_calibration_paths == 200000L),
  all(pair$base_seed == 20260805L),
  all(pair$run_version == "v20260805_supplemental_v2"),
  all(pair$safeguard_rate == 0)
)
siga_assert_close(
  max(abs(pair$difference_siga_s_minus_rt_pp)), 0.900,
  label = "supplemental pair-path SIGA-S maximum absolute difference (pp)"
)
siga_assert_close(
  max(abs(pair$difference_siga_r_minus_rt_pp)), 0.244,
  label = "supplemental pair-path SIGA-R maximum absolute difference (pp)"
)
siga_assert_close(
  max(pair$mean_variance_ratio_r_over_s), 1.19359594637313,
  tolerance = 1e-13,
  label = "supplemental pair-path maximum mean variance ratio"
)
gates <- read.csv(
  file.path(
    siga_raw_dir(root), "pair_path_supplemental",
    "pair_path_supplemental_gates.csv"
  ),
  check.names = FALSE
)
stopifnot(nrow(gates) == 24L, all(gates$primary_pass))

swift <- siga_build_swift(root)
stopifnot(nrow(swift) == 2L, all(swift$n_outer == 100000L))
siga_assert_close(swift$siga_s_percent, c(81.483, 81.535), label = "SWIFT SIGA-S power (%)")
siga_assert_close(swift$siga_r_percent, c(81.440, 81.477), label = "SWIFT SIGA-R power (%)")
siga_assert_close(swift$rt_percent, c(81.283, 81.497), label = "SWIFT RT power (%)")
stopifnot(all(swift$safeguard_rate == 0),
          max(swift$max_abs_siga_s_pvalue_reproduction_error) < 2e-15)

audit <- read.csv(file.path(siga_raw_dir(root), "swift_direct", "swift_direct_siga_r_audit.csv"), check.names = FALSE)
stopifnot(audit$value[audit$item == "completed outer trials"] == 100000,
          audit$value[audit$item == "pair-calibration replicates"] == 100000)

timing <- read.csv(file.path(siga_raw_dir(root), "timing", "full_grid_timing.csv"), check.names = FALSE)
cal_timing <- read.csv(file.path(siga_raw_dir(root), "timing", "allocation_calibration_timing.csv"), check.names = FALSE)
stopifnot(nrow(timing) == 48L, nrow(cal_timing) == 4L,
          all(timing$benchmark_repeats == 3L), all(cal_timing$allocation_calibration_paths == 100000L))

std_dir <- file.path(root, "data", "standardized_results")
siga_compare_data_frames(s, read.csv(file.path(std_dir, "full_grid_siga_s_standardized.csv"), check.names = FALSE),
                         c("outcome", "factor_count", "target_per_group", "scenario_id", "analysis"))
siga_compare_data_frames(r, read.csv(file.path(std_dir, "full_grid_siga_r_standardized.csv"), check.names = FALSE),
                         c("scenario_id", "analysis"))
siga_compare_data_frames(pair, read.csv(file.path(std_dir, "pair_path_supplemental_standardized.csv"), check.names = FALSE),
                         c("scenario_id", "analysis"))
siga_compare_data_frames(swift, read.csv(file.path(std_dir, "swift_direct_standardized.csv"), check.names = FALSE),
                         c("analysis"))

cat("PASS: production aggregate outputs reproduce all headline manuscript values.\n")
cat("  SIGA-S full grid: 112 comparisons; maxima 0.422, 0.185, 0.422 pp.\n")
cat("  SIGA-R full grid: 112 comparisons; maxima 0.322, 0.205, 0.322 pp.\n")
cat("  Supplemental pair-path sensitivity: 12 scenarios; maxima 0.900 and 0.244 pp; all gates passed.\n")
cat("  SWIFT DIRECT-inspired: 100,000 outer trials and exact reported powers.\n")
