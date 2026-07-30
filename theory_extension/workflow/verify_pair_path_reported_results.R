#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L) {
  stop("Usage: Rscript verify_pair_path_reported_results.R SUMMARY.csv [EXPECTED.csv]", call. = FALSE)
}
summary_path <- normalizePath(args[1L], mustWork = TRUE)
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1L]))) else getwd()
project_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
expected_path <- if (length(args) == 2L) normalizePath(args[2L], mustWork = TRUE) else {
  normalizePath(file.path(project_root, "expected", "reported_pair_path_validation.csv"), mustWork = TRUE)
}

observed <- read.csv(summary_path, stringsAsFactors = FALSE, check.names = FALSE)
expected <- read.csv(expected_path, stringsAsFactors = FALSE, check.names = FALSE)
key <- c("scenario_id", "analysis")
required <- c(
  key, "mean_variance_ratio_r_over_s", "predicted_rt_size_from_mean_ratio",
  "siga_s", "siga_r", "rt", "difference_siga_s_minus_rt", "difference_siga_r_minus_rt"
)
missing <- setdiff(required, names(observed))
if (length(missing)) stop("Missing summary columns: ", paste(missing, collapse = ", "), call. = FALSE)

merged <- merge(expected, observed, by = key, all.x = TRUE, suffixes = c("_expected", "_observed"))
if (anyNA(merged$scenario_id) || nrow(merged) != nrow(expected)) stop("Could not match every expected row.", call. = FALSE)

checks <- list(
  mean_variance_ratio_r_over_s = c("mean_variance_ratio_r_over_s", 1, 0.0011),
  predicted_rt_percent = c("predicted_rt_size_from_mean_ratio", 100, 0.011),
  siga_s_percent = c("siga_s", 100, 0.011),
  siga_r_percent = c("siga_r", 100, 0.011),
  rt_percent = c("rt", 100, 0.011),
  difference_siga_s_minus_rt_percent = c("difference_siga_s_minus_rt", 100, 0.011),
  difference_siga_r_minus_rt_percent = c("difference_siga_r_minus_rt", 100, 0.011)
)

failures <- character()
for (expected_name in names(checks)) {
  observed_name <- checks[[expected_name]][1L]
  scale <- as.numeric(checks[[expected_name]][2L])
  tolerance <- as.numeric(checks[[expected_name]][3L])
  expected_column <- if (expected_name %in% names(merged)) expected_name else paste0(expected_name, "_expected")
  observed_column <- if (observed_name %in% names(merged)) observed_name else paste0(observed_name, "_observed")
  if (!expected_column %in% names(merged) || !observed_column %in% names(merged)) {
    stop("Could not resolve verification columns for ", expected_name, call. = FALSE)
  }
  delta <- abs(as.numeric(merged[[expected_column]]) - scale * as.numeric(merged[[observed_column]]))
  cat(sprintf("%-43s max abs difference = %.6f (tolerance %.6f)\n", expected_name, max(delta, na.rm = TRUE), tolerance))
  if (any(!is.finite(delta)) || any(delta > tolerance)) failures <- c(failures, expected_name)
}
if (length(failures)) stop("Reported-value verification failed for: ", paste(failures, collapse = ", "), call. = FALSE)
cat("PASS: aggregated pair-path results reproduce the reported values at displayed precision.\n")
