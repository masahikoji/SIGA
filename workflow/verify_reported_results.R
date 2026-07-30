#!/usr/bin/env Rscript

# Verify that full-run outputs reproduce the values printed in the manuscript.
# Operating characteristics are compared after the same rounding used in the
# manuscript. Timing is hardware-dependent and is checked only when
# SIGA_VERIFY_TIMING=1.

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
expected_dir <- file.path(root_dir, "expected")
actual_dir <- file.path(root_dir, "manuscript_outputs")

read_required <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path, call. = FALSE)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

failures <- character()
record_failure <- function(text) failures <<- c(failures, text)

compare_joined <- function(actual, expected, keys, numeric_rules, label) {
  if (anyDuplicated(actual[keys])) record_failure(paste0(label, ": duplicate keys in actual output."))
  if (anyDuplicated(expected[keys])) record_failure(paste0(label, ": duplicate keys in expected output."))
  merged <- merge(expected, actual, by = keys, all = TRUE, suffixes = c(".expected", ".actual"), sort = FALSE)
  if (nrow(merged) != nrow(expected) || nrow(merged) != nrow(actual)) {
    record_failure(sprintf("%s: key mismatch (expected %d rows, actual %d rows, merged %d rows).",
                           label, nrow(expected), nrow(actual), nrow(merged)))
    return(invisible(NULL))
  }
  for (column in names(numeric_rules)) {
    digits <- numeric_rules[[column]]
    e <- merged[[paste0(column, ".expected")]]
    a <- merged[[paste0(column, ".actual")]]
    bad <- is.na(e) != is.na(a) | (!is.na(e) & round(as.numeric(e), digits) != round(as.numeric(a), digits))
    if (any(bad)) {
      rows <- which(bad)
      preview <- head(rows, 8L)
      detail <- paste(vapply(preview, function(i) {
        key_values <- unlist(merged[i, keys, drop = FALSE], use.names = FALSE)
        key_text <- paste(paste0(keys, "=", key_values), collapse = ", ")
        sprintf("[%s: expected %s, actual %s]", key_text,
                format(e[i], digits = 12), format(a[i], digits = 12))
      }, character(1L)), collapse = "; ")
      record_failure(sprintf("%s: %s differs after rounding to %d decimals. %s",
                             label, column, digits, detail))
    }
  }
  invisible(NULL)
}

# Operating characteristics: manuscript reports percentages to two decimals.
actual_oc <- read_required(file.path(actual_dir, "combined_operating_characteristics.csv"), "Actual operating characteristics")
expected_oc <- read_required(file.path(expected_dir, "reported_operating_characteristics.csv"), "Expected operating characteristics")
keys_oc <- c("outcome", "factor_count", "target_per_group", "objective", "scenario_role", "true_effect")
actual_oc_pct <- actual_oc
expected_oc_pct <- expected_oc
probability_columns <- c("proposed_unadjusted", "reference_unadjusted", "proposed_adjusted", "reference_adjusted")
actual_oc_pct[probability_columns] <- lapply(actual_oc_pct[probability_columns], function(x) 100 * as.numeric(x))
expected_oc_pct[probability_columns] <- lapply(expected_oc_pct[probability_columns], function(x) 100 * as.numeric(x))
compare_joined(
  actual_oc_pct, expected_oc_pct, keys_oc,
  setNames(as.list(rep(2L, length(probability_columns))), probability_columns),
  "Operating characteristics"
)

# Check the two headline maximum absolute differences after manuscript rounding.
metric_actual <- read_required(file.path(actual_dir, "manuscript_summary_metrics.csv"), "Actual summary metrics")
metric_expected <- data.frame(
  outcome = c("Continuous", "Binary"),
  maximum_absolute_difference_percentage_points = c(0.14, 0.42),
  stringsAsFactors = FALSE
)
compare_joined(
  metric_actual, metric_expected, "outcome",
  list(maximum_absolute_difference_percentage_points = 2L),
  "Headline maximum differences"
)

# SWIFT DIRECT-inspired power and timing.
actual_case <- read_required(file.path(actual_dir, "swift_direct_case_study_summary.csv"), "Actual SWIFT DIRECT-inspired summary")
expected_case_power <- read_required(file.path(expected_dir, "reported_swift_direct_power.csv"), "Expected SWIFT DIRECT-inspired power")
actual_case_power <- data.frame(
  analysis = c("unadjusted", "adjusted"),
  siga_power = c(actual_case$proposed_unadjusted[1L], actual_case$proposed_adjusted[1L]),
  rt_power = c(actual_case$reference_unadjusted[1L], actual_case$reference_adjusted[1L]),
  stringsAsFactors = FALSE
)
compare_joined(actual_case_power, expected_case_power, "analysis",
               list(siga_power = 4L, rt_power = 4L), "SWIFT DIRECT-inspired power")

expected_case_timing <- read_required(file.path(expected_dir, "reported_swift_direct_timing.csv"), "Expected SWIFT DIRECT-inspired timing")
actual_case_timing <- data.frame(
  component = c("allocation_calibration", "complete_unadjusted_adjusted_analyses", "total_method_specific"),
  siga_minutes = c(actual_case$calibration_minutes[1L], actual_case$proposed_analysis_minutes[1L], actual_case$proposed_total_minutes[1L]),
  rt_minutes = c(NA_real_, actual_case$reference_analysis_minutes[1L], actual_case$reference_analysis_minutes[1L]),
  stringsAsFactors = FALSE
)
# This benchmark is machine-dependent; compare it only on request.
verify_timing <- identical(Sys.getenv("SIGA_VERIFY_TIMING", unset = "0"), "1")
if (verify_timing) {
  compare_joined(actual_case_timing, expected_case_timing, "component",
                 list(siga_minutes = 2L, rt_minutes = 2L), "SWIFT DIRECT-inspired timing")
}

# General timing tables are also hardware-dependent.
actual_timing <- read_required(file.path(actual_dir, "timing_summary_complete.csv"), "Actual timing summary")
expected_timing <- read_required(file.path(expected_dir, "reported_timing_m3_ultra.csv"), "Expected M3 Ultra timing summary")
keys_timing <- c("outcome", "analysis", "factor_count", "target_per_group", "objective")
if (nrow(actual_timing) != 48L) record_failure(paste0("Timing summary: expected 48 rows, found ", nrow(actual_timing), "."))
if (verify_timing) {
  compare_joined(actual_timing, expected_timing, keys_timing,
                 list(siga_projected_minutes_100k = 3L, rt_projected_minutes_100k = 2L),
                 "M3 Ultra timing")

  actual_calibration <- read_required(file.path(actual_dir, "calibration_timing_repeats.csv"), "Actual calibration timing")
  expected_calibration <- read_required(file.path(expected_dir, "reported_calibration_timing_m3_ultra.csv"), "Expected calibration timing")
  actual_calibration_wide <- reshape(
    actual_calibration[, c("factor_count", "target_per_group", "total_n", "allocation_calibration_paths", "calibration_repeat", "calibration_seconds")],
    idvar = c("factor_count", "target_per_group", "total_n", "allocation_calibration_paths"),
    timevar = "calibration_repeat", direction = "wide"
  )
  names(actual_calibration_wide) <- sub("calibration_seconds\\.", "run", names(actual_calibration_wide))
  names(actual_calibration_wide)[names(actual_calibration_wide) == "allocation_calibration_paths"] <- "B0"
  run_columns <- intersect(c("run1", "run2", "run3"), names(actual_calibration_wide))
  actual_calibration_wide$mean_seconds <- rowMeans(actual_calibration_wide[run_columns])
  compare_joined(actual_calibration_wide, expected_calibration,
                 c("factor_count", "target_per_group", "total_n", "B0"),
                 list(run1 = 3L, run2 = 3L, run3 = 3L, mean_seconds = 3L),
                 "M3 Ultra allocation calibration timing")
}

if (length(failures)) {
  cat("REPRODUCIBILITY CHECK FAILED\n")
  cat(paste0("- ", failures, collapse = "\n"), "\n")
  quit(status = 1L)
}

cat("PASS: operating characteristics and SWIFT DIRECT-inspired power reproduce the manuscript-rounded values.\n")
if (verify_timing) {
  cat("PASS: M3 Ultra timing values also reproduce the manuscript-rounded values.\n")
} else {
  cat("NOTE: machine-dependent timing values were structurally checked but not compared. Set SIGA_VERIFY_TIMING=1 on the Apple M3 Ultra benchmark machine to compare them.\n")
}
