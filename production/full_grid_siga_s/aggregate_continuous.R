#!/usr/bin/env Rscript

# =============================================================================
# Standalone aggregation for the final fixed-score SIGA simulations
# =============================================================================
# This script never reruns the simulation. It reads the completed scenario CSV
# files, verifies 100,000 unique outer replicates for scenarios 1--28, and
# recreates the final summary CSV and manuscript-ready TeX tables.
#
# Usage from any directory:
#   Rscript aggregate_all_siga_fixed_score_power85.R
#
# Optional override:
#   PWRT_PROJECT_DIR=/path/to/02_program Rscript aggregate_all_siga_fixed_score_power85.R
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

script_directory <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

PROJECT_DIR <- path.expand(Sys.getenv("PWRT_PROJECT_DIR", unset = script_directory()))
PROJECT_DIR <- normalizePath(PROJECT_DIR, winslash = "/", mustWork = TRUE)
EXPECTED_SCENARIOS <- 1:28
EXPECTED_OUTER <- 100000L
N_RERANDOMIZATIONS <- 4999L
N_CALIBRATION <- 100000L
POWER_LOWER <- 0.80
POWER_UPPER <- 0.90

methods <- list(
  list(
    method_id = "unadjusted_siga_vs_fixed_score_rt_power85",
    adjusted = FALSE,
    output_dir = file.path(PROJECT_DIR, "unadjusted_siga_vs_fixed_score_rt_power85_output")
  ),
  list(
    method_id = "covariate_adjusted_siga_vs_fixed_score_rt_power85",
    adjusted = TRUE,
    output_dir = file.path(PROJECT_DIR, "covariate_adjusted_siga_vs_fixed_score_rt_power85_output")
  )
)

make_scenarios <- function() {
  designs <- data.frame(
    factor_count = c(2L, 2L, 5L, 5L),
    target_per_group = c(100L, 500L, 200L, 1000L),
    superiority_power_effect = c(0.430, 0.190, 0.300, 0.135),
    noninferiority_power_effect = c(0.230, -0.010, 0.100, -0.065),
    equivalence_power_effect = c(0.000, 0.280, 0.180, 0.330),
    stringsAsFactors = FALSE
  )

  roles <- data.frame(
    objective = c(
      "superiority", "superiority",
      "noninferiority", "noninferiority",
      "equivalence", "equivalence", "equivalence"
    ),
    scenario_role = c(
      "type1", "power",
      "type1", "power",
      "type1_lower", "type1_upper", "power"
    ),
    effect_code = c(
      "superiority_type1", "superiority_power",
      "noninferiority_type1", "noninferiority_power",
      "equivalence_type1_lower", "equivalence_type1_upper", "equivalence_power"
    ),
    alpha = c(0.05, 0.05, 0.025, 0.025, 0.05, 0.05, 0.05),
    stringsAsFactors = FALSE
  )

  rows <- vector("list", nrow(designs) * nrow(roles))
  pos <- 0L
  for (d in seq_len(nrow(designs))) {
    for (r in seq_len(nrow(roles))) {
      pos <- pos + 1L
      effect <- switch(
        roles$effect_code[r],
        superiority_type1 = 0,
        superiority_power = designs$superiority_power_effect[d],
        noninferiority_type1 = -0.20,
        noninferiority_power = designs$noninferiority_power_effect[d],
        equivalence_type1_lower = -0.45,
        equivalence_type1_upper = 0.45,
        equivalence_power = designs$equivalence_power_effect[d]
      )
      rows[[pos]] <- data.frame(
        factor_count = designs$factor_count[d],
        target_per_group = designs$target_per_group[d],
        total_n = 2L * designs$target_per_group[d],
        objective = roles$objective[r],
        scenario_role = roles$scenario_role[r],
        true_effect = effect,
        alpha = roles$alpha[r],
        stringsAsFactors = FALSE
      )
    }
  }

  out <- do.call(rbind, rows)
  out$scenario_id <- seq_len(nrow(out))
  out$scenario_label <- paste0(
    "K", out$factor_count,
    "_npg", out$target_per_group,
    "_", out$objective,
    "_", out$scenario_role
  )
  out
}

SCENARIOS <- make_scenarios()

wilson_interval <- function(x, n, conf.level = 0.95) {
  z <- qnorm(1 - (1 - conf.level) / 2)
  p <- x / n
  den <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(estimate = p, lower = max(0, center - half), upper = min(1, center + half))
}

role_label <- function(x) {
  switch(
    x,
    type1 = "Type I error",
    power = "Power",
    type1_lower = "Type I error at lower limit",
    type1_upper = "Type I error at upper limit",
    x
  )
}

objective_label <- function(x) {
  switch(
    x,
    superiority = "Superiority",
    noninferiority = "Non-inferiority",
    equivalence = "Equivalence",
    x
  )
}

read_csv_safely <- function(path) {
  tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) {
    stop("Could not read: ", path, "\n", conditionMessage(e), call. = FALSE)
  })
}

get_calibration_info <- function(factor_count, total_n) {
  cache_dir <- file.path(PROJECT_DIR, "siga_allocation_calibration_cache")
  if (!dir.exists(cache_dir)) return(list(paths = N_CALIBRATION, minutes = NA_real_))
  pattern <- sprintf("^siga_cal_K%d_n%d_.*_B%d_seed[0-9]+[.]rds$", factor_count, total_n, N_CALIBRATION)
  files <- list.files(cache_dir, pattern = pattern, full.names = TRUE)
  if (!length(files)) return(list(paths = N_CALIBRATION, minutes = NA_real_))
  info <- file.info(files)
  file <- files[which.max(info$mtime)]
  calibration <- tryCatch(readRDS(file), error = function(e) NULL)
  if (is.null(calibration)) return(list(paths = N_CALIBRATION, minutes = NA_real_))
  list(
    paths = if (!is.null(calibration$B0)) calibration$B0 else N_CALIBRATION,
    minutes = if (!is.null(calibration$elapsed_seconds)) calibration$elapsed_seconds / 60 else NA_real_
  )
}

summarize_one_scenario <- function(method, scenario) {
  chunk_dir <- file.path(method$output_dir, "scenario_shards")
  prefix <- sprintf("scenario_%02d_", scenario$scenario_id)
  dirs <- list.dirs(chunk_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[startsWith(basename(dirs), prefix)]
  if (!length(dirs)) {
    stop("Missing scenario directory for ", method$method_id,
         ", scenario ", scenario$scenario_id, call. = FALSE)
  }

  detail_files <- unlist(lapply(dirs, function(d) {
    list.files(d, pattern = "^shard_[0-9]+_of_[0-9]+[.]csv$", full.names = TRUE)
  }))
  if (!length(detail_files)) {
    stop("No detail CSV for ", method$method_id,
         ", scenario ", scenario$scenario_id, call. = FALSE)
  }

  detail <- do.call(rbind, lapply(detail_files, read_csv_safely))
  detail <- detail[detail$scenario_id == scenario$scenario_id, , drop = FALSE]
  detail <- detail[!duplicated(detail$replicate), , drop = FALSE]
  detail <- detail[order(detail$replicate), , drop = FALSE]

  if (nrow(detail) != EXPECTED_OUTER) {
    stop(method$method_id, ", scenario ", scenario$scenario_id,
         " has ", nrow(detail), " unique replicates; expected ", EXPECTED_OUTER,
         call. = FALSE)
  }
  if (!identical(as.integer(detail$replicate), seq_len(EXPECTED_OUTER))) {
    missing_ids <- setdiff(seq_len(EXPECTED_OUTER), detail$replicate)
    stop(method$method_id, ", scenario ", scenario$scenario_id,
         " has nonconsecutive replicate IDs. First missing IDs: ",
         paste(head(missing_ids, 10L), collapse = ","), call. = FALSE)
  }

  required <- c(
    "reject_proposed", "reject_randomization",
    "p_proposed", "p_randomization",
    "data_generation_seconds", "proposed_analysis_seconds",
    "randomization_test_seconds"
  )
  missing_columns <- setdiff(required, names(detail))
  if (length(missing_columns)) {
    stop("Missing columns in ", method$method_id, ", scenario ", scenario$scenario_id,
         ": ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }

  observed_effects <- unique(round(as.numeric(detail$true_effect), 12))
  if (length(observed_effects) != 1L || abs(observed_effects - scenario$true_effect) > 1e-10) {
    stop("True-effect mismatch for ", method$method_id,
         ", scenario ", scenario$scenario_id, call. = FALSE)
  }

  n <- nrow(detail)
  proposed_ci <- wilson_interval(sum(detail$reject_proposed), n)
  rt_ci <- wilson_interval(sum(detail$reject_randomization), n)

  timing_files <- unlist(lapply(dirs, function(d) {
    list.files(d, pattern = "_timing[.]csv$", full.names = TRUE)
  }))
  timing <- NULL
  if (length(timing_files)) {
    timing <- do.call(rbind, lapply(timing_files, read_csv_safely))
    timing <- timing[!duplicated(
      timing[c("scenario_id", "shard_id", "n_shards")], fromLast = TRUE
    ), , drop = FALSE]
  }

  calibration <- get_calibration_info(scenario$factor_count, scenario$total_n)
  reference_estimate <- unname(rt_ci["estimate"])

  data.frame(
    method_id = method$method_id,
    adjusted = method$adjusted,
    reference_statistic = "fixed treatment-score statistic",
    scenario_id = scenario$scenario_id,
    scenario_label = scenario$scenario_label,
    factor_count = scenario$factor_count,
    target_per_group = scenario$target_per_group,
    total_n = scenario$total_n,
    objective = scenario$objective,
    scenario_role = scenario$scenario_role,
    true_effect = scenario$true_effect,
    alpha = scenario$alpha,
    requested_outer_trials = EXPECTED_OUTER,
    completed_outer_trials = n,
    rerandomizations_per_trial = N_RERANDOMIZATIONS,
    calibration_paths = calibration$paths,
    calibration_minutes = calibration$minutes,
    proposed_rejection = unname(proposed_ci["estimate"]),
    proposed_lower_95 = unname(proposed_ci["lower"]),
    proposed_upper_95 = unname(proposed_ci["upper"]),
    randomization_rejection = reference_estimate,
    randomization_lower_95 = unname(rt_ci["lower"]),
    randomization_upper_95 = unname(rt_ci["upper"]),
    reference_power_in_target_band = if (scenario$scenario_role == "power") {
      reference_estimate >= POWER_LOWER && reference_estimate <= POWER_UPPER
    } else {
      NA
    },
    pvalue_mean_absolute_error = mean(abs(detail$p_proposed - detail$p_randomization)),
    pvalue_root_mean_squared_error = sqrt(mean((detail$p_proposed - detail$p_randomization)^2)),
    pvalue_correlation = cor(detail$p_proposed, detail$p_randomization),
    data_generation_minutes_sum = sum(detail$data_generation_seconds) / 60,
    proposed_analysis_minutes_sum = sum(detail$proposed_analysis_seconds) / 60,
    randomization_test_minutes_sum = sum(detail$randomization_test_seconds) / 60,
    proposed_minutes_mean_per_trial = mean(detail$proposed_analysis_seconds) / 60,
    randomization_minutes_mean_per_trial = mean(detail$randomization_test_seconds) / 60,
    shard_wall_clock_minutes_sum = if (is.null(timing)) NA_real_ else sum(timing$shard_wall_clock_seconds) / 60,
    shard_wall_clock_minutes_max = if (is.null(timing)) NA_real_ else max(timing$shard_wall_clock_seconds) / 60,
    stringsAsFactors = FALSE
  )
}

write_tex_results <- function(summary, method, path) {
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\small",
    "\\renewcommand{\\arraystretch}{1.08}",
    "\\begin{tabular}{rrllrrr}",
    "\\toprule",
    "Factors & $n$/group & Objective & Scenario & True effect & Proposed (\\%) & Randomization test (\\%) \\\\",
    "\\midrule"
  )

  for (i in seq_len(nrow(summary))) {
    row <- summary[i, ]
    lines <- c(lines, paste0(
      row$factor_count, " & ",
      row$target_per_group, " & ",
      objective_label(row$objective), " & ",
      role_label(row$scenario_role), " & ",
      sprintf("%+.3f", row$true_effect), " & ",
      sprintf("%.2f", 100 * row$proposed_rejection), " & ",
      sprintf("%.2f", 100 * row$randomization_rejection), " \\\\"
    ))
  }

  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    paste0(
      "\\caption{Operating characteristics for the proposed SIGA approximation and the fixed-score conditional randomization test. Both procedures use the same boundary-specific statistic $T_n(r)=\\sum_i(A_i-1/2)r_i$. Each reference test uses ",
      N_RERANDOMIZATIONS, " regenerated allocation paths.}"
    ),
    paste0("\\label{tab:", method$method_id, "-operating-characteristics}"),
    "\\end{table}",
    "",
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\small",
    "\\renewcommand{\\arraystretch}{1.08}",
    "\\begin{tabular}{rrllrr}",
    "\\toprule",
    "Factors & $n$/group & Objective & Scenario & Proposed time (min) & Randomization-test time (min) \\\\",
    "\\midrule"
  )

  for (i in seq_len(nrow(summary))) {
    row <- summary[i, ]
    lines <- c(lines, paste0(
      row$factor_count, " & ",
      row$target_per_group, " & ",
      objective_label(row$objective), " & ",
      role_label(row$scenario_role), " & ",
      sprintf("%.2f", row$proposed_analysis_minutes_sum), " & ",
      sprintf("%.2f", row$randomization_test_minutes_sum), " \\\\"
    ))
  }

  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    paste0(
      "\\caption{Cumulative analysis time across all outer trials, reported in minutes. The one-time allocation-only calibration used ",
      N_CALIBRATION, " paths; its time is reported separately in minutes in the CSV output.}"
    ),
    paste0("\\label{tab:", method$method_id, "-timing}"),
    "\\end{table}"
  )

  writeLines(lines, con = path, useBytes = TRUE)
}

aggregate_method <- function(method) {
  if (!dir.exists(method$output_dir)) {
    stop("Output directory not found: ", method$output_dir, call. = FALSE)
  }

  message("Aggregating ", method$method_id, " ...")
  rows <- lapply(seq_len(nrow(SCENARIOS)), function(i) {
    summarize_one_scenario(method, SCENARIOS[i, ])
  })
  summary <- do.call(rbind, rows)
  summary <- summary[order(summary$scenario_id), , drop = FALSE]

  if (!identical(as.integer(summary$scenario_id), EXPECTED_SCENARIOS)) {
    stop("Scenario IDs are not exactly 1--28 for ", method$method_id, call. = FALSE)
  }
  if (any(summary$completed_outer_trials != EXPECTED_OUTER)) {
    stop("At least one scenario is incomplete for ", method$method_id, call. = FALSE)
  }

  csv_path <- file.path(method$output_dir, paste0(method$method_id, "_summary.csv"))
  tex_path <- file.path(method$output_dir, paste0(method$method_id, "_results.tex"))
  write.csv(summary, csv_path, row.names = FALSE)
  write_tex_results(summary, method, tex_path)

  power_rows <- summary$scenario_role == "power"
  outside <- summary[power_rows & !summary$reference_power_in_target_band, , drop = FALSE]
  if (nrow(outside)) {
    warning(
      method$method_id, ": reference power outside 80%--90% in scenario(s): ",
      paste(outside$scenario_id, collapse = ", "),
      ". Inspect the pilot before changing effects.",
      call. = FALSE
    )
  }

  message("  CSV: ", csv_path)
  message("  TeX: ", tex_path)
  message("  Scenarios: ", nrow(summary))
  message("  Completed range: ", min(summary$completed_outer_trials),
          "--", max(summary$completed_outer_trials))
  summary
}

all_summaries <- lapply(methods, aggregate_method)
combined <- do.call(rbind, all_summaries)
combined_path <- file.path(PROJECT_DIR, "all_siga_fixed_score_power85_summary.csv")
write.csv(combined, combined_path, row.names = FALSE)
message("Combined CSV: ", combined_path)
message("All fixed-score SIGA aggregation completed successfully.")
