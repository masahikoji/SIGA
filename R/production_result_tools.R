siga_repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    script <- normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = FALSE)
    return(normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

siga_raw_dir <- function(root) file.path(root, "data", "production_results", "raw")

siga_build_full_grid_s <- function(root) {
  d <- file.path(siga_raw_dir(root), "full_grid_siga_s")
  cu <- read.csv(file.path(d, "continuous_unadjusted_summary.csv"), check.names = FALSE)
  ca <- read.csv(file.path(d, "continuous_adjusted_summary.csv"), check.names = FALSE)
  bs <- read.csv(file.path(d, "binary_small_n_summary.csv"), check.names = FALSE)
  bl <- read.csv(file.path(d, "binary_large_n_summary.csv"), check.names = FALSE)

  continuous <- function(x, analysis, source_file) data.frame(
    benchmark = "sampling_targeted", outcome = "continuous", analysis = analysis,
    scenario_id = x$scenario_id, scenario_label = x$scenario_label,
    factor_count = x$factor_count, target_per_group = x$target_per_group,
    total_n = x$total_n, objective = x$objective, scenario_role = x$scenario_role,
    true_effect = x$true_effect, alpha = x$alpha,
    n_outer = x$completed_outer_trials,
    rerandomizations_per_trial = x$rerandomizations_per_trial,
    calibration_paths = x$calibration_paths,
    siga_rejection = x$proposed_rejection, rt_rejection = x$randomization_rejection,
    difference = x$proposed_rejection - x$randomization_rejection,
    difference_pp = 100 * (x$proposed_rejection - x$randomization_rejection),
    source_file = source_file, stringsAsFactors = FALSE
  )
  binary <- function(x, source_file) {
    u <- data.frame(
      benchmark = "sampling_targeted", outcome = "binary", analysis = "Unadjusted",
      scenario_id = x$scenario_id, scenario_label = x$scenario_label,
      factor_count = x$factor_count, target_per_group = x$target_per_group,
      total_n = x$total_n, objective = x$objective, scenario_role = x$scenario_role,
      true_effect = x$true_risk_difference, alpha = x$alpha, n_outer = x$n_outer,
      rerandomizations_per_trial = x$rerandomizations_per_trial,
      calibration_paths = if ("allocation_calibration_paths" %in% names(x)) x$allocation_calibration_paths else 100000L,
      siga_rejection = x$proposed_unadjusted, rt_rejection = x$reference_unadjusted,
      difference = x$proposed_unadjusted - x$reference_unadjusted,
      difference_pp = 100 * (x$proposed_unadjusted - x$reference_unadjusted),
      source_file = source_file, stringsAsFactors = FALSE
    )
    a <- u
    a$analysis <- "Adjusted"
    a$siga_rejection <- x$proposed_adjusted
    a$rt_rejection <- x$reference_adjusted
    a$difference <- x$proposed_adjusted - x$reference_adjusted
    a$difference_pp <- 100 * a$difference
    rbind(u, a)
  }
  out <- rbind(
    continuous(cu, "Unadjusted", "continuous_unadjusted_summary.csv"),
    continuous(ca, "Adjusted", "continuous_adjusted_summary.csv"),
    binary(bl, "binary_large_n_summary.csv"),
    binary(bs, "binary_small_n_summary.csv")
  )
  out[order(out$outcome, out$factor_count, out$target_per_group, out$scenario_id, out$analysis), , drop = FALSE]
}

siga_build_full_grid_r <- function(root) {
  x <- read.csv(file.path(siga_raw_dir(root), "full_grid_siga_r", "siga_r_full_benchmark_summary.csv"), check.names = FALSE)
  x$benchmark <- "randomization_targeted"
  x$difference_pp <- 100 * x$siga_r_minus_rt
  x
}

siga_build_pair_stress <- function(root) {
  x <- read.csv(file.path(siga_raw_dir(root), "pair_path_stress", "pair_path_stress_summary.csv"), check.names = FALSE)
  x$difference_siga_s_minus_rt_pp <- 100 * x$difference_siga_s_minus_rt
  x$difference_siga_r_minus_rt_pp <- 100 * x$difference_siga_r_minus_rt
  x
}

siga_build_swift <- function(root) {
  x <- read.csv(file.path(siga_raw_dir(root), "swift_direct", "swift_direct_siga_r_summary.csv"), check.names = FALSE)
  x$siga_s_percent <- 100 * x$siga_s_power
  x$siga_r_percent <- 100 * x$siga_r_power
  x$rt_percent <- 100 * x$rt_power
  x$siga_r_minus_rt_pp <- 100 * x$siga_r_minus_rt
  x
}

siga_assert_close <- function(actual, expected, tolerance = 1e-12, label = "value") {
  if (length(actual) != length(expected) || any(!is.finite(actual)) ||
      any(abs(actual - expected) > tolerance)) {
    stop(label, " mismatch. Actual: ", paste(actual, collapse = ", "),
         "; expected: ", paste(expected, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

siga_compare_data_frames <- function(a, b, key, tolerance = 1e-12) {
  a <- a[do.call(order, a[key]), , drop = FALSE]
  b <- b[do.call(order, b[key]), , drop = FALSE]
  rownames(a) <- rownames(b) <- NULL
  if (!identical(names(a), names(b))) stop("Column names differ.", call. = FALSE)
  if (nrow(a) != nrow(b)) stop("Row counts differ.", call. = FALSE)
  for (nm in names(a)) {
    if (is.numeric(a[[nm]]) || is.integer(a[[nm]])) {
      xa <- as.numeric(a[[nm]]); xb <- as.numeric(b[[nm]])
      ok_na <- is.na(xa) & is.na(xb)
      bad <- !(ok_na | (!is.na(xa) & !is.na(xb) & abs(xa - xb) <= tolerance))
      if (any(bad)) stop("Numeric mismatch in column ", nm, call. = FALSE)
    } else if (!identical(as.character(a[[nm]]), as.character(b[[nm]]))) {
      stop("Character mismatch in column ", nm, call. = FALSE)
    }
  }
  invisible(TRUE)
}
