#!/usr/bin/env Rscript

# Build the manuscript-level result files from the frozen production summaries.
# This script performs no simulation. It validates the full-run settings and
# transforms the simulation summaries into one common data structure.

options(stringsAsFactors = FALSE, warn = 1)

script_directory <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

root_dir <- normalizePath(
  path.expand(Sys.getenv("SIGA_PROJECT_DIR", unset = file.path(script_directory(), ".."))),
  winslash = "/", mustWork = TRUE
)
output_dir <- file.path(root_dir, "manuscript_outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop("Could not create output directory: ", output_dir, call. = FALSE)

first_existing <- function(paths, label) {
  paths <- path.expand(paths[nzchar(paths)])
  hit <- paths[file.exists(paths)]
  if (!length(hit)) stop(label, " not found. Checked: ", paste(paths, collapse = "; "), call. = FALSE)
  normalizePath(hit[1L], winslash = "/", mustWork = TRUE)
}

read_checked <- function(path, required, label) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  x
}

coerce_logical_strict <- function(x, label) {
  if (is.logical(x)) {
    if (anyNA(x)) stop(label, " contains NA.", call. = FALSE)
    return(x)
  }
  if (is.numeric(x) || is.integer(x)) {
    if (anyNA(x) || any(!x %in% c(0, 1))) stop(label, " must contain 0/1.", call. = FALSE)
    return(x == 1)
  }
  z <- tolower(trimws(as.character(x)))
  ans <- rep(NA, length(z))
  ans[z %in% c("true", "t", "1", "yes")] <- TRUE
  ans[z %in% c("false", "f", "0", "no")] <- FALSE
  if (anyNA(ans)) stop(label, " cannot be interpreted as logical.", call. = FALSE)
  ans
}

expected_outer <- 100000L
expected_rerandomizations <- 4999L
expected_calibration <- 100000L

continuous_path <- first_existing(c(
  Sys.getenv("SIGA_CONTINUOUS_SUMMARY", unset = ""),
  file.path(root_dir, "all_siga_fixed_score_power85_summary.csv")
), "Continuous summary")

binary_large_path <- first_existing(c(
  Sys.getenv("SIGA_BINARY_LARGE_SUMMARY", unset = ""),
  file.path(root_dir, "binary_siga_weak_null_pureR_100k_4999_output", "binary_siga_weak_null_summary.csv"),
  file.path(root_dir, "binary_siga_weak_null_pureR_100K_4999_output", "binary_siga_weak_null_summary.csv")
), "Large-sample binary summary")

binary_small_path <- first_existing(c(
  Sys.getenv("SIGA_BINARY_SMALL_SUMMARY", unset = ""),
  file.path(root_dir, "binary_siga_weak_null_small_n_100K_4999_output", "binary_siga_weak_null_small_n_summary.csv")
), "Small-sample binary summary")

timing_path <- first_existing(c(
  Sys.getenv("SIGA_TIMING_SUMMARY", unset = ""),
  file.path(root_dir, "standardized_timing_rerun_output", "clean_timing_benchmark_summary_split_complete.csv")
), "Complete timing summary")

calibration_timing_path <- first_existing(c(
  Sys.getenv("SIGA_CALIBRATION_TIMING", unset = ""),
  file.path(root_dir, "standardized_timing_rerun_output", "allocation_calibration_repeat_times.csv")
), "Calibration timing repeats")

case_path <- first_existing(c(
  Sys.getenv("SIGA_CASE_SUMMARY", unset = ""),
  file.path(root_dir, "swift_direct_siga_case_study_pureR_output", "swift_direct_siga_case_study_summary.csv")
), "SWIFT DIRECT-inspired summary")

# Continuous results: two separately run analysis scripts, paired by scenario.
continuous <- read_checked(
  continuous_path,
  c("adjusted", "reference_statistic", "scenario_id", "factor_count", "target_per_group",
    "objective", "scenario_role", "true_effect", "alpha", "completed_outer_trials",
    "rerandomizations_per_trial", "proposed_rejection", "randomization_rejection"),
  "Continuous summary"
)
continuous$adjusted <- coerce_logical_strict(continuous$adjusted, "continuous adjusted")
if (nrow(continuous) != 56L) stop("Expected 56 continuous rows; found ", nrow(continuous), call. = FALSE)
if (any(continuous$reference_statistic != "fixed treatment-score statistic")) {
  stop("Continuous comparator is not uniformly the fixed-score randomization test.", call. = FALSE)
}
if (any(as.integer(continuous$completed_outer_trials) != expected_outer)) stop("Incomplete continuous results.", call. = FALSE)
if (any(as.integer(continuous$rerandomizations_per_trial) != expected_rerandomizations)) stop("Continuous rerandomization count mismatch.", call. = FALSE)
keys <- c("scenario_id", "factor_count", "target_per_group", "objective", "scenario_role", "true_effect", "alpha")
cu <- continuous[!continuous$adjusted, c(keys, "proposed_rejection", "randomization_rejection")]
ca <- continuous[ continuous$adjusted, c(keys, "proposed_rejection", "randomization_rejection")]
names(cu)[names(cu) == "proposed_rejection"] <- "proposed_unadjusted"
names(cu)[names(cu) == "randomization_rejection"] <- "reference_unadjusted"
names(ca)[names(ca) == "proposed_rejection"] <- "proposed_adjusted"
names(ca)[names(ca) == "randomization_rejection"] <- "reference_adjusted"
continuous_wide <- merge(cu, ca, by = keys, sort = FALSE)
if (nrow(continuous_wide) != 28L) stop("Expected 28 paired continuous rows.", call. = FALSE)
continuous_wide$outcome <- "Continuous"

# Binary large-sample results are already jointly summarized by the production script.
binary_large <- read_checked(
  binary_large_path,
  c("scenario_id", "factor_count", "target_per_group", "objective", "scenario_role",
    "true_risk_difference", "n_outer", "rerandomizations_per_trial",
    "allocation_calibration_paths", "proposed_unadjusted", "reference_unadjusted",
    "proposed_adjusted", "reference_adjusted"),
  "Large-sample binary summary"
)
if (nrow(binary_large) != 14L) stop("Expected 14 large-sample binary rows.", call. = FALSE)
if (any(as.integer(binary_large$n_outer) != expected_outer)) stop("Incomplete large-sample binary results.", call. = FALSE)
if (any(as.integer(binary_large$rerandomizations_per_trial) != expected_rerandomizations)) stop("Large-sample binary rerandomization count mismatch.", call. = FALSE)
if (any(as.integer(binary_large$allocation_calibration_paths) != expected_calibration)) stop("Large-sample binary calibration count mismatch.", call. = FALSE)
binary_large$true_effect <- binary_large$true_risk_difference

# Binary small-sample results use the manuscript-frozen margins selected by the pilot.
binary_small <- read_checked(
  binary_small_path,
  c("scenario_id", "factor_count", "target_per_group", "objective", "scenario_role",
    "true_risk_difference", "n_outer", "rerandomizations_per_trial",
    "allocation_calibration_paths", "proposed_unadjusted", "reference_unadjusted",
    "proposed_adjusted", "reference_adjusted"),
  "Small-sample binary summary"
)
if (nrow(binary_small) != 14L) stop("Expected 14 small-sample binary rows.", call. = FALSE)
if (any(as.integer(binary_small$n_outer) != expected_outer)) stop("Incomplete small-sample binary results.", call. = FALSE)
if (any(as.integer(binary_small$rerandomizations_per_trial) != expected_rerandomizations)) stop("Small-sample binary rerandomization count mismatch.", call. = FALSE)
if (any(as.integer(binary_small$allocation_calibration_paths) != expected_calibration)) stop("Small-sample binary calibration count mismatch.", call. = FALSE)
binary_small$true_effect <- binary_small$true_risk_difference

binary_columns <- c("factor_count", "target_per_group", "objective", "scenario_role", "true_effect",
                    "proposed_unadjusted", "reference_unadjusted", "proposed_adjusted", "reference_adjusted")
binary <- rbind(binary_small[, binary_columns], binary_large[, binary_columns])
binary$outcome <- "Binary"
if (nrow(binary) != 28L) stop("Expected 28 binary rows after combining sample sizes.", call. = FALSE)

final_columns <- c("outcome", "factor_count", "target_per_group", "objective", "scenario_role", "true_effect",
                   "proposed_unadjusted", "reference_unadjusted", "proposed_adjusted", "reference_adjusted")
combined <- rbind(continuous_wide[, final_columns], binary[, final_columns])
objective_order <- c(superiority = 1L, noninferiority = 2L, equivalence = 3L)
role_order <- c(type1 = 1L, type1_lower = 1L, type1_upper = 2L, power = 3L)
outcome_order <- c(Continuous = 1L, Binary = 2L)
combined <- combined[order(combined$factor_count, outcome_order[combined$outcome], combined$target_per_group,
                           objective_order[combined$objective], role_order[combined$scenario_role]), ]
rownames(combined) <- NULL
write.csv(combined, file.path(output_dir, "combined_operating_characteristics.csv"), row.names = FALSE)

# Timing: the complete from-scratch run contains 48 rows (all four designs for
# both outcomes, two analyses, and three objectives).
timing <- read_checked(
  timing_path,
  c("outcome", "analysis", "factor_count", "target_per_group", "objective",
    "benchmark_repeats", "trials_per_repeat", "siga_seconds_per_trial", "rt_seconds_per_trial",
    "siga_projected_minutes_100k", "rt_projected_minutes_100k"),
  "Timing summary"
)
if (nrow(timing) != 48L) stop("Expected 48 complete timing rows; found ", nrow(timing), call. = FALSE)
if (any(as.integer(timing$benchmark_repeats) != 3L)) stop("Timing repeats are not uniformly three.", call. = FALSE)
if (any(as.integer(timing$trials_per_repeat) != 30L)) stop("Timing trials per repeat are not uniformly 30.", call. = FALSE)
write.csv(timing, file.path(output_dir, "timing_summary_complete.csv"), row.names = FALSE)

calibration_timing <- read_checked(
  calibration_timing_path,
  c("factor_count", "target_per_group", "total_n", "calibration_repeat",
    "allocation_calibration_paths", "calibration_seconds"),
  "Calibration timing repeats"
)
if (nrow(calibration_timing) != 12L) stop("Expected 12 calibration timing rows (four designs x three repeats).", call. = FALSE)
write.csv(calibration_timing, file.path(output_dir, "calibration_timing_repeats.csv"), row.names = FALSE)

case_summary <- read_checked(
  case_path,
  c("factor_count", "target_per_group", "total_n", "objective", "scenario_role", "n_outer",
    "proposed_unadjusted", "reference_unadjusted", "proposed_adjusted", "reference_adjusted",
    "calibration_minutes", "proposed_analysis_minutes", "proposed_total_minutes", "reference_analysis_minutes"),
  "SWIFT DIRECT-inspired summary"
)
if (nrow(case_summary) != 1L) stop("Expected one SWIFT DIRECT-inspired scenario.", call. = FALSE)
if (as.integer(case_summary$n_outer) != expected_outer) stop("Case-study outer-trial count mismatch.", call. = FALSE)
write.csv(case_summary, file.path(output_dir, "swift_direct_case_study_summary.csv"), row.names = FALSE)

# Summary metrics quoted in the manuscript.
max_abs_diff <- function(dat) {
  max(c(abs(dat$proposed_unadjusted - dat$reference_unadjusted),
        abs(dat$proposed_adjusted - dat$reference_adjusted)))
}
summary_metrics <- data.frame(
  outcome = c("Continuous", "Binary"),
  maximum_absolute_difference_percentage_points = c(
    100 * max_abs_diff(combined[combined$outcome == "Continuous", ]),
    100 * max_abs_diff(combined[combined$outcome == "Binary", ])
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_metrics, file.path(output_dir, "manuscript_summary_metrics.csv"), row.names = FALSE)

message("Manuscript-level outputs written to: ", output_dir)
message("Operating-characteristic rows: ", nrow(combined))
message("Timing rows: ", nrow(timing))
