#!/usr/bin/env Rscript

# =============================================================================
# Standardized single-core timing benchmark for SIGA and the fixed-score
# randomization test under biased-coin Pocock--Simon minimization.
#
# Purpose
#   1. Re-measure all continuous-outcome timings after an explicit warm-up,
#      avoiding first-call overhead in the first benchmarked scenario.
#   2. Add binary-outcome timings for 2 factors with 100 participants/group
#      and 5 factors with 200 participants/group.
#   3. Exclude allocation-only calibration, data generation, and score
#      construction from the measured analysis time.
#   4. Measure the one-time allocation-only calibration cost for each unique
#      minimization design using fresh, uncached B0-path calibrations.
#   5. Produce summaries compatible with clean_timing_benchmark_summary_split.csv.
#
# The script reuses only function definitions from
# simulation_binary_SIGA_small_n_margin_tuning_v3.R; it does not execute that
# script's simulation main block.
#
# Recommended run on the M3 Ultra:
#   caffeinate -dimsu Rscript rerun_standardized_timing_benchmark_v3_with_calibration_fixed.R
#
# Important environment variables
#   PWRT_PROJECT_DIR          Project directory.
#   PWRT_ENGINE_SCRIPT       Path to simulation_binary_SIGA_small_n_margin_tuning_v3.R.
#   PWRT_EXISTING_TIMING_CSV Existing clean timing summary to merge with.
#   PWRT_TIMING_SCOPE        "needed" (default), "continuous", "binary_small", or "all".
#                            needed      = all continuous + binary n/group 100 and 200.
#                            continuous  = all continuous designs only.
#                            binary_small = binary n/group 100 and 200 only.
#                            all         = all continuous and all binary designs.
#   PWRT_TIMING_REPEATS      Independent timing repeats; default 3.
#   PWRT_TIMING_TRIALS       Prepared trials per repeat; default 30.
#   PWRT_TIMING_RERAND       Regenerated paths per randomization test; default 4999.
#   PWRT_TIMING_CALIBRATION  Allocation-only paths; default 100000.
#   PWRT_SIGA_MIN_SECONDS    Minimum cumulative SIGA timing per repeat; default 1.0.
#   PWRT_SIGA_MIN_CALLS      Minimum SIGA calls per repeat; default 1000.
#   PWRT_TIME_CALIBRATION    1 to time fresh allocation calibrations; default 1.
#   PWRT_CALIBRATION_REPEATS Number of fresh full-B0 calibration repeats; default 3.
#   PWRT_CALIBRATION_WARMUP_PATHS Small unmeasured warm-up calibration; default 1000.
#   PWRT_TIMING_OUTPUT_DIR   Output directory.
#
# Output
#   clean_timing_repeat_means_rerun.csv
#   clean_timing_benchmark_summary_rerun.csv
#   clean_timing_benchmark_summary_split_complete.csv
#   allocation_calibration_repeat_times.csv
#   allocation_calibration_timing_summary.csv
#   allocation_calibration_timing_table.tex
#   clean_timing_benchmark_with_calibration_complete.csv
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

# Keep the benchmark single-threaded, including BLAS where supported.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  RCPP_PARALLEL_NUM_THREADS = "1"
)

script_directory <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(
      sub("^--file=", "", file_arg[1L]),
      winslash = "/", mustWork = FALSE
    )))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

env_integer <- function(name, default) {
  value <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) return(as.integer(default))
  out <- suppressWarnings(as.integer(value))
  if (is.na(out)) stop(name, " must be an integer.", call. = FALSE)
  out
}

env_numeric <- function(name, default) {
  value <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) return(as.numeric(default))
  out <- suppressWarnings(as.numeric(value))
  if (!is.finite(out)) stop(name, " must be numeric.", call. = FALSE)
  out
}

ensure_directory <- function(path) {
  path <- path.expand(path)
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop("Could not create directory: ", path, call. = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

PROJECT_DIR <- path.expand(Sys.getenv(
  "PWRT_PROJECT_DIR",
  unset = "~/Dropbox/Tex/03_Post_PhD/029_exact_minimization/02_program"
))
PROJECT_DIR <- ensure_directory(PROJECT_DIR)

ENGINE_SCRIPT <- path.expand(Sys.getenv(
  "PWRT_ENGINE_SCRIPT",
  unset = file.path(PROJECT_DIR, "simulation_binary_SIGA_small_n_margin_tuning_v3.R")
))
if (!file.exists(ENGINE_SCRIPT)) {
  fallback <- file.path(script_directory(), "simulation_binary_SIGA_small_n_margin_tuning_v3.R")
  if (file.exists(fallback)) ENGINE_SCRIPT <- fallback
}
if (!file.exists(ENGINE_SCRIPT)) {
  stop(
    "Could not find simulation_binary_SIGA_small_n_margin_tuning_v3.R. ",
    "Set PWRT_ENGINE_SCRIPT to its full path.", call. = FALSE
  )
}

OUTPUT_DIR <- ensure_directory(path.expand(Sys.getenv(
  "PWRT_TIMING_OUTPUT_DIR",
  unset = file.path(PROJECT_DIR, "standardized_timing_rerun_output")
)))

EXISTING_TIMING_CSV <- path.expand(Sys.getenv(
  "PWRT_EXISTING_TIMING_CSV",
  unset = file.path(PROJECT_DIR, "clean_timing_benchmark_summary_split.csv")
))
if (!file.exists(EXISTING_TIMING_CSV)) {
  alternative <- file.path(
    PROJECT_DIR,
    "final_operating_characteristics_tables_100K_4999",
    "clean_timing_benchmark_summary_split.csv"
  )
  if (file.exists(alternative)) EXISTING_TIMING_CSV <- alternative
}

TIMING_SCOPE <- tolower(trimws(Sys.getenv("PWRT_TIMING_SCOPE", unset = "needed")))
if (!TIMING_SCOPE %in% c("needed", "continuous", "binary_small", "all")) {
  stop(
    "PWRT_TIMING_SCOPE must be 'needed', 'continuous', 'binary_small', or 'all'.",
    call. = FALSE
  )
}

N_REPEATS <- env_integer("PWRT_TIMING_REPEATS", 3L)
TRIALS_PER_REPEAT <- env_integer("PWRT_TIMING_TRIALS", 30L)
N_RERAND <- env_integer("PWRT_TIMING_RERAND", 4999L)
N_CALIBRATION <- env_integer("PWRT_TIMING_CALIBRATION", 100000L)
SIGA_MIN_SECONDS <- env_numeric("PWRT_SIGA_MIN_SECONDS", 1.0)
SIGA_MIN_CALLS <- env_integer("PWRT_SIGA_MIN_CALLS", 1000L)
BASE_SEED <- env_integer("PWRT_TIMING_SEED", 20260729L)
P_BIASED_COIN <- env_numeric("PWRT_P_BIASED_COIN", 0.80)
TIME_CALIBRATION <- identical(Sys.getenv("PWRT_TIME_CALIBRATION", unset = "1"), "1")
CALIBRATION_REPEATS <- env_integer("PWRT_CALIBRATION_REPEATS", 3L)
CALIBRATION_WARMUP_PATHS <- env_integer("PWRT_CALIBRATION_WARMUP_PATHS", 1000L)

stopifnot(
  N_REPEATS >= 2L,
  TRIALS_PER_REPEAT >= 2L,
  N_RERAND >= 1L,
  N_CALIBRATION >= 2L,
  SIGA_MIN_SECONDS > 0,
  SIGA_MIN_CALLS >= 1L,
  P_BIASED_COIN > 0.5,
  P_BIASED_COIN < 1,
  CALIBRATION_REPEATS >= 1L,
  CALIBRATION_WARMUP_PATHS >= 2L
)

# -----------------------------------------------------------------------------
# Load only the required function definitions from the existing binary engine.
# This avoids executing its main simulation block and keeps the timing analysis
# identical to the existing SIGA and randomization-test implementations.
# -----------------------------------------------------------------------------

REQUIRED_ENGINE_FUNCTIONS <- c(
  "%||%",
  "check_binary_factor_matrix",
  "joint_stratum_id",
  "nearest_psd",
  "normalise_seed",
  "validate_minimization_inputs",
  "ps_assign_R",
  "precompute_design_R",
  "rt_pvalues_R",
  "calibrate_siga_binary",
  "all_binary_patterns",
  "pattern_probabilities",
  "calibrate_binary_outcome_model",
  "generate_factors",
  "binary_event_probabilities",
  "generate_binary_outcome",
  "make_boundary_scores",
  "score_decomposition",
  "siga_variance",
  "gaussian_pvalue",
  "lattice_normal_pvalue",
  "siga_component_pvalue",
  "reference_component_pvalue",
  "combine_objective_pvalues"
)

load_named_function_definitions <- function(path, function_names, envir = .GlobalEnv) {
  expressions <- parse(file = path, keep.source = FALSE)
  found <- character()
  for (expr in expressions) {
    if (!is.call(expr) || !identical(expr[[1L]], as.name("<-"))) next
    lhs <- expr[[2L]]
    rhs <- expr[[3L]]
    lhs_name <- if (is.symbol(lhs)) as.character(lhs) else NA_character_
    is_function <- is.call(rhs) && identical(rhs[[1L]], as.name("function"))
    if (!is.na(lhs_name) && lhs_name %in% function_names && is_function) {
      eval(expr, envir = envir)
      found <- c(found, lhs_name)
    }
  }
  missing <- setdiff(function_names, unique(found))
  if (length(missing)) {
    stop(
      "The engine script is missing required function definitions: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  invisible(unique(found))
}

load_named_function_definitions(ENGINE_SCRIPT, REQUIRED_ENGINE_FUNCTIONS)

# -----------------------------------------------------------------------------
# Design and outcome settings.
# -----------------------------------------------------------------------------

FACTOR_PREVALENCE_MASTER <- c(0.50, 0.40, 0.30, 0.20, 0.10)
CONTINUOUS_BETA_MASTER <- c(0.50, 0.40, 0.30, 0.20, 0.10)
CONTINUOUS_INTERACTION <- 0.25
BINARY_BETA_MASTER <- c(0.35, -0.25, 0.20, -0.15, 0.10)
BINARY_INTERACTION <- 0.20
TARGET_CONTROL_RISK <- 0.60

continuous_effects <- data.frame(
  factor_count = c(2L, 2L, 5L, 5L),
  target_per_group = c(100L, 500L, 200L, 1000L),
  superiority = c(0.430, 0.190, 0.300, 0.135),
  noninferiority = c(0.230, -0.010, 0.100, -0.065),
  equivalence = c(0.000, 0.280, 0.180, 0.330)
)

binary_effects <- data.frame(
  factor_count = c(2L, 2L, 5L, 5L),
  target_per_group = c(100L, 500L, 200L, 1000L),
  superiority = c(0.200, 0.100, 0.140, 0.070),
  noninferiority = c(0.100, 0.000, 0.050, -0.030),
  equivalence = c(0.000, 0.000, 0.000, 0.040)
)

all_designs <- rbind(
  data.frame(outcome = "Continuous", continuous_effects),
  data.frame(outcome = "Binary", binary_effects)
)
if (TIMING_SCOPE == "needed") {
  benchmark_designs <- all_designs[
    all_designs$outcome == "Continuous" |
      (all_designs$outcome == "Binary" &
         all_designs$target_per_group %in% c(100L, 200L)),
    , drop = FALSE
  ]
} else if (TIMING_SCOPE == "continuous") {
  benchmark_designs <- all_designs[
    all_designs$outcome == "Continuous",
    , drop = FALSE
  ]
} else if (TIMING_SCOPE == "binary_small") {
  benchmark_designs <- all_designs[
    all_designs$outcome == "Binary" &
      all_designs$target_per_group %in% c(100L, 200L),
    , drop = FALSE
  ]
} else {
  benchmark_designs <- all_designs
}
rownames(benchmark_designs) <- NULL

factor_probabilities <- function(factor_count) {
  FACTOR_PREVALENCE_MASTER[seq_len(factor_count)]
}

get_effect <- function(design, objective) {
  as.numeric(design[[objective]])
}

# Use the selected small-n equivalence margins when they can be identified.
# The values do not materially affect timing, but using them keeps the timing
# benchmark aligned with the final analysis specification.
selected_parameter_candidates <- c(
  file.path(PROJECT_DIR, "binary_small_n_selected_parameters.csv"),
  file.path(PROJECT_DIR, "binary_small_n_selected_alternatives.csv")
)
selected_parameter_path <- selected_parameter_candidates[file.exists(selected_parameter_candidates)][1L]
selected_parameters <- NULL
if (length(selected_parameter_path) && !is.na(selected_parameter_path)) {
  selected_parameters <- tryCatch(
    read.csv(selected_parameter_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
}

first_existing_column <- function(data, candidates) {
  hits <- candidates[candidates %in% names(data)]
  if (!length(hits)) return(NULL)
  hits[[1L]]
}

finite_numeric_scalar <- function(data, column_name) {
  if (is.null(column_name) || length(column_name) != 1L ||
      is.na(column_name) || !nzchar(column_name) ||
      !column_name %in% names(data)) {
    return(NA_real_)
  }
  value <- suppressWarnings(as.numeric(data[[column_name]][1L]))
  if (length(value) != 1L || !is.finite(value)) return(NA_real_)
  value
}

extract_equivalence_limits <- function(factor_count, target_per_group) {
  default_limits <- c(-0.10, 0.10)
  if (is.null(selected_parameters)) return(default_limits)

  required_keys <- c("factor_count", "target_per_group")
  if (!all(required_keys %in% names(selected_parameters))) return(default_limits)

  factor_values <- suppressWarnings(as.integer(selected_parameters$factor_count))
  sample_values <- suppressWarnings(as.integer(selected_parameters$target_per_group))
  matching_rows <- which(
    !is.na(factor_values) & !is.na(sample_values) &
      factor_values == as.integer(factor_count) &
      sample_values == as.integer(target_per_group)
  )
  if (!length(matching_rows)) return(default_limits)
  row <- selected_parameters[matching_rows[[1L]], , drop = FALSE]

  lower_name <- first_existing_column(
    row,
    c("equivalence_lower", "equivalence_lower_limit", "selected_equivalence_lower")
  )
  upper_name <- first_existing_column(
    row,
    c("equivalence_upper", "equivalence_upper_limit", "selected_equivalence_upper")
  )
  lower_value <- finite_numeric_scalar(row, lower_name)
  upper_value <- finite_numeric_scalar(row, upper_name)
  if (is.finite(lower_value) && is.finite(upper_value) &&
      lower_value < upper_value) {
    return(c(lower_value, upper_value))
  }

  margin_name <- first_existing_column(
    row,
    c(
      "equivalence_margin", "selected_equivalence_margin",
      "equivalence_half_width", "equivalence_limit"
    )
  )
  margin <- abs(finite_numeric_scalar(row, margin_name))
  if (is.finite(margin) && margin > 0) return(c(-margin, margin))

  default_limits
}

objective_spec <- function(outcome, objective, factor_count, target_per_group) {
  if (outcome == "Continuous") {
    return(switch(
      objective,
      superiority = list(boundaries = 0, alternatives = "two.sided"),
      noninferiority = list(boundaries = -0.20, alternatives = "greater"),
      equivalence = list(
        boundaries = c(-0.45, 0.45),
        alternatives = c("greater", "less")
      ),
      stop("Unknown continuous objective: ", objective, call. = FALSE)
    ))
  }

  switch(
    objective,
    superiority = list(boundaries = 0, alternatives = "two.sided"),
    noninferiority = list(boundaries = -0.10, alternatives = "greater"),
    equivalence = list(
      boundaries = extract_equivalence_limits(factor_count, target_per_group),
      alternatives = c("greater", "less")
    ),
    stop("Unknown binary objective: ", objective, call. = FALSE)
  )
}

# Fail fast before any long timing run. This specifically validates all binary
# objective specifications, including the optional selected-margin CSV.
for (test_row in seq_len(nrow(binary_effects))) {
  for (test_objective in c("superiority", "noninferiority", "equivalence")) {
    test_spec <- objective_spec(
      outcome = "Binary",
      objective = test_objective,
      factor_count = binary_effects$factor_count[test_row],
      target_per_group = binary_effects$target_per_group[test_row]
    )
    if (!is.list(test_spec) ||
        !length(test_spec$boundaries) ||
        !length(test_spec$alternatives) ||
        length(test_spec$boundaries) != length(test_spec$alternatives) ||
        any(!is.finite(as.numeric(test_spec$boundaries)))) {
      stop(
        "Binary objective preflight check failed for K=",
        binary_effects$factor_count[test_row],
        ", n/group=", binary_effects$target_per_group[test_row],
        ", objective=", test_objective,
        call. = FALSE
      )
    }
  }
}
message("Preflight checks passed for all binary objective specifications.")

# -----------------------------------------------------------------------------
# Allocation-only calibration cache. This is performed outside the timed code.
# The same cache can be reused by continuous and binary outcomes because the
# calibration depends only on the allocation design.
# -----------------------------------------------------------------------------

CALIBRATION_DIR <- ensure_directory(file.path(
  PROJECT_DIR, "binary_siga_allocation_calibration_cache_pureR"
))

calibration_cache_file <- function(factor_count, total_n) {
  probs <- factor_probabilities(factor_count)
  prob_tag <- paste(formatC(probs, format = "f", digits = 3), collapse = "-")
  file.path(
    CALIBRATION_DIR,
    sprintf(
      "binary_siga_pureR_K%d_n%d_p%s_pbc%.3f_B%d.rds",
      factor_count, total_n, prob_tag, P_BIASED_COIN, N_CALIBRATION
    )
  )
}

calibration_seed <- function(total_n, factor_count) {
  modulus <- 2147483646
  value <- as.double(BASE_SEED) + 7000003 +
    1009 * as.double(total_n) + 104729 * as.double(factor_count)
  as.integer(value %% modulus + 1)
}

get_design_calibration <- function(factor_count, total_n) {
  cache <- calibration_cache_file(factor_count, total_n)
  if (file.exists(cache)) {
    message("Using cached allocation calibration: ", basename(cache))
    return(readRDS(cache))
  }
  message(
    "Creating allocation-only calibration outside timed code: K=", factor_count,
    ", total n=", total_n, ", B0=", N_CALIBRATION
  )
  calibration <- calibrate_siga_binary(
    n = total_n,
    factor_prob = factor_probabilities(factor_count),
    pbc = P_BIASED_COIN,
    weights = rep(1, factor_count + 1L),
    B0 = N_CALIBRATION,
    seed = calibration_seed(total_n, factor_count)
  )
  saveRDS(calibration, cache)
  calibration
}

# -----------------------------------------------------------------------------
# Fresh allocation-only calibration timing. These measurements deliberately do
# not read or write the calibration cache. They quantify the one-time offline
# cost of estimating the allocation-law covariance (and treated-count
# probabilities used by the binary lattice refinement) for a design.
# -----------------------------------------------------------------------------

fresh_calibration_seed <- function(total_n, factor_count, repeat_id) {
  modulus <- 2147483646
  value <- as.double(BASE_SEED) + 170000033 +
    1009 * as.double(total_n) +
    104729 * as.double(factor_count) +
    1000003 * as.double(repeat_id)
  as.integer(value %% modulus + 1)
}

time_fresh_calibration <- function(factor_count, total_n, repeat_id) {
  gc(verbose = FALSE)
  start <- proc.time()[["elapsed"]]
  calibration <- calibrate_siga_binary(
    n = total_n,
    factor_prob = factor_probabilities(factor_count),
    pbc = P_BIASED_COIN,
    weights = rep(1, factor_count + 1L),
    B0 = N_CALIBRATION,
    seed = fresh_calibration_seed(total_n, factor_count, repeat_id)
  )
  elapsed <- proc.time()[["elapsed"]] - start
  rm(calibration)
  gc(verbose = FALSE)
  elapsed
}

warm_up_calibration <- function(factor_count, total_n) {
  message(
    "  Warm-up calibration (unmeasured): K=", factor_count,
    ", total n=", total_n, ", B0=", CALIBRATION_WARMUP_PATHS
  )
  invisible(calibrate_siga_binary(
    n = total_n,
    factor_prob = factor_probabilities(factor_count),
    pbc = P_BIASED_COIN,
    weights = rep(1, factor_count + 1L),
    B0 = CALIBRATION_WARMUP_PATHS,
    seed = fresh_calibration_seed(total_n, factor_count, 0L)
  ))
}

# -----------------------------------------------------------------------------
# Trial generation and score construction. These are completed before timing.
# -----------------------------------------------------------------------------

seed_value <- function(outcome, factor_count, target_per_group, objective,
                       repeat_id, trial_id, stream_id) {
  modulus <- 2147483646
  outcome_code <- if (outcome == "Continuous") 11 else 29
  objective_code <- match(objective, c("superiority", "noninferiority", "equivalence"))
  value <- as.double(BASE_SEED) +
    1000003 * stream_id +
    104729 * outcome_code +
    1009 * factor_count +
    10007 * target_per_group +
    100003 * objective_code +
    10000019 * repeat_id +
    97 * trial_id
  as.integer(value %% modulus + 1)
}

generate_continuous_outcome <- function(X, A, factor_prob, tau) {
  factor_count <- ncol(X)
  beta <- CONTINUOUS_BETA_MASTER[seq_len(factor_count)]
  centered_main <- sweep(X, 2L, factor_prob, FUN = "-")
  mean_interaction <- factor_prob[1L] * factor_prob[2L]
  baseline <- drop(centered_main %*% beta) +
    CONTINUOUS_INTERACTION * (X[, 1L] * X[, 2L] - mean_interaction) +
    rnorm(nrow(X))
  baseline + tau * A
}

prepare_trials <- function(design, objective, repeat_id) {
  factor_count <- design$factor_count
  target_per_group <- design$target_per_group
  total_n <- 2L * target_per_group
  factor_prob <- factor_probabilities(factor_count)
  tau_or_rd <- get_effect(design, objective)
  spec <- objective_spec(
    design$outcome, objective, factor_count, target_per_group
  )

  binary_model <- NULL
  if (design$outcome == "Binary") {
    binary_model <- calibrate_binary_outcome_model(
      factor_prob = factor_prob,
      target_control_risk = TARGET_CONTROL_RISK,
      target_risk_difference = tau_or_rd,
      beta = BINARY_BETA_MASTER,
      interaction = BINARY_INTERACTION
    )
  }

  trials <- vector("list", TRIALS_PER_REPEAT)
  for (trial_id in seq_len(TRIALS_PER_REPEAT)) {
    set.seed(seed_value(
      design$outcome, factor_count, target_per_group, objective,
      repeat_id, trial_id, 1L
    ))
    X <- generate_factors(total_n, factor_prob)
    z <- ps_assign_R(
      X = X,
      pbc = P_BIASED_COIN,
      weights_ = rep(1, factor_count + 1L),
      seed = seed_value(
        design$outcome, factor_count, target_per_group, objective,
        repeat_id, trial_id, 2L
      )
    )
    A <- as.integer((z + 1L) / 2L)
    set.seed(seed_value(
      design$outcome, factor_count, target_per_group, objective,
      repeat_id, trial_id, 3L
    ))
    y <- if (design$outcome == "Continuous") {
      generate_continuous_outcome(X, A, factor_prob, tau_or_rd)
    } else {
      generate_binary_outcome(X, A, binary_model)
    }
    scores <- make_boundary_scores(
      y = y,
      A = A,
      X = X,
      boundaries = spec$boundaries
    )
    trials[[trial_id]] <- list(
      X = X,
      z = as.integer(z),
      A = A,
      y = y,
      scores = scores,
      spec = spec
    )
  }
  trials
}

# -----------------------------------------------------------------------------
# Timed analysis functions.
# -----------------------------------------------------------------------------

analyze_siga_fixed_score <- function(trial, outcome, objective, analysis,
                                     calibration) {
  score_matrix <- if (analysis == "Unadjusted") {
    trial$scores$unadjusted
  } else {
    trial$scores$adjusted
  }
  component_p <- numeric(ncol(score_matrix))
  for (j in seq_len(ncol(score_matrix))) {
    use_lattice <- outcome == "Binary" &&
      objective == "superiority" &&
      analysis == "Unadjusted" &&
      abs(trial$spec$boundaries[j]) < 1e-14
    result <- siga_component_pvalue(
      score = score_matrix[, j],
      z = trial$z,
      X = trial$X,
      calibration = calibration,
      alternative = trial$spec$alternatives[j],
      use_binary_lattice = use_lattice,
      y = if (use_lattice) trial$y else NULL
    )
    component_p[j] <- result$p
  }
  combine_objective_pvalues(component_p, objective)
}

analyze_rt_fixed_score <- function(trial, objective, analysis,
                                   factor_count, seed) {
  score_matrix <- if (analysis == "Unadjusted") {
    trial$scores$unadjusted
  } else {
    trial$scores$adjusted
  }
  score_matrix <- as.matrix(score_matrix)
  rt <- rt_pvalues_R(
    X = trial$X,
    scores = score_matrix,
    z_obs = trial$z,
    B = N_RERAND,
    pbc = P_BIASED_COIN,
    weights_ = rep(1, factor_count + 1L),
    seed = seed
  )
  component_p <- numeric(ncol(score_matrix))
  for (j in seq_len(ncol(score_matrix))) {
    component_p[j] <- reference_component_pvalue(
      rt, j, trial$spec$alternatives[j]
    )
  }
  combine_objective_pvalues(component_p, objective)
}

time_siga_repeat <- function(trials, outcome, objective, analysis, calibration) {
  # Explicit warm-up removes first-call allocation, dispatch, and cache effects.
  warmup_n <- min(10L, length(trials))
  for (i in seq_len(warmup_n)) {
    invisible(analyze_siga_fixed_score(
      trials[[i]], outcome, objective, analysis, calibration
    ))
  }

  calls <- 0L
  start <- proc.time()[["elapsed"]]
  elapsed <- 0
  while (elapsed < SIGA_MIN_SECONDS || calls < SIGA_MIN_CALLS) {
    for (trial in trials) {
      invisible(analyze_siga_fixed_score(
        trial, outcome, objective, analysis, calibration
      ))
      calls <- calls + 1L
    }
    elapsed <- proc.time()[["elapsed"]] - start
  }
  list(
    total_seconds = elapsed,
    calls = calls,
    seconds_per_trial = elapsed / calls
  )
}

time_rt_repeat <- function(trials, design, objective, analysis, repeat_id) {
  # One unmeasured full call is sufficient to remove first-call overhead.
  invisible(analyze_rt_fixed_score(
    trials[[1L]], objective, analysis, design$factor_count,
    seed = seed_value(
      design$outcome, design$factor_count, design$target_per_group,
      objective, repeat_id, 0L, 91L
    )
  ))

  start <- proc.time()[["elapsed"]]
  for (trial_id in seq_along(trials)) {
    invisible(analyze_rt_fixed_score(
      trials[[trial_id]], objective, analysis, design$factor_count,
      seed = seed_value(
        design$outcome, design$factor_count, design$target_per_group,
        objective, repeat_id, trial_id, 92L
      )
    ))
  }
  elapsed <- proc.time()[["elapsed"]] - start
  list(
    total_seconds = elapsed,
    calls = length(trials),
    seconds_per_trial = elapsed / length(trials)
  )
}

# -----------------------------------------------------------------------------
# Run benchmark.
# -----------------------------------------------------------------------------

message("Engine script: ", ENGINE_SCRIPT)
message("Output directory: ", OUTPUT_DIR)
message("Timing scope: ", TIMING_SCOPE)
message(
  "Protocol: ", N_REPEATS, " repeats x ", TRIALS_PER_REPEAT,
  " prepared trials; B=", N_RERAND,
  "; SIGA timed for >=", SIGA_MIN_SECONDS, " s and >=", SIGA_MIN_CALLS,
  " calls per repeat."
)
message("Calibration, data generation, and score construction are outside per-trial analysis timing.")
if (TIME_CALIBRATION) {
  message(
    "Fresh allocation-only calibration will also be timed separately: ",
    CALIBRATION_REPEATS, " repeat(s), B0=", N_CALIBRATION,
    " paths per design."
  )
}

calibration_cache <- new.env(parent = emptyenv())
get_cached_calibration <- function(factor_count, total_n) {
  key <- paste(factor_count, total_n, sep = ":")
  if (!exists(key, envir = calibration_cache, inherits = FALSE)) {
    assign(
      key,
      get_design_calibration(factor_count, total_n),
      envir = calibration_cache
    )
  }
  get(key, envir = calibration_cache, inherits = FALSE)
}

objectives <- c("superiority", "noninferiority", "equivalence")
analyses <- c("Unadjusted", "Adjusted")

# Save every completed analysis row immediately. A later error or interruption
# therefore does not discard completed timing measurements, and rerunning the
# same command resumes from the checkpoint.
checkpoint_path <- file.path(OUTPUT_DIR, "timing_repeat_checkpoint.csv")
RESUME_CHECKPOINT <- identical(
  Sys.getenv("PWRT_RESUME_TIMING", unset = "1"), "1"
)

row_key <- function(dat) {
  paste(
    dat$outcome, dat$analysis, dat$factor_count,
    dat$target_per_group, dat$objective, dat$repeat_id,
    sep = "|"
  )
}

write_checkpoint <- function(dat, path) {
  tmp <- tempfile(
    pattern = "timing_checkpoint_",
    tmpdir = dirname(path),
    fileext = ".csv"
  )
  write.csv(dat, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) {
    if (!file.copy(tmp, path, overwrite = TRUE)) {
      unlink(tmp)
      stop("Could not write timing checkpoint: ", path, call. = FALSE)
    }
    unlink(tmp)
  }
  invisible(NULL)
}

checkpoint_results <- NULL
if (RESUME_CHECKPOINT && file.exists(checkpoint_path)) {
  checkpoint_results <- read.csv(
    checkpoint_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_checkpoint_columns <- c(
    "outcome", "analysis", "factor_count", "target_per_group",
    "objective", "repeat_id", "trials_per_repeat",
    "rerandomizations_per_rt", "allocation_calibration_paths"
  )
  if (!all(required_checkpoint_columns %in% names(checkpoint_results))) {
    stop(
      "Existing checkpoint has an incompatible structure: ", checkpoint_path,
      ". Delete it or run with PWRT_RESUME_TIMING=0.",
      call. = FALSE
    )
  }
  protocol_matches <-
    all(checkpoint_results$trials_per_repeat == TRIALS_PER_REPEAT) &&
    all(checkpoint_results$rerandomizations_per_rt == N_RERAND) &&
    all(checkpoint_results$allocation_calibration_paths == N_CALIBRATION)
  if (!protocol_matches) {
    stop(
      "Existing checkpoint was created under a different timing protocol: ",
      checkpoint_path,
      ". Delete it or run with PWRT_RESUME_TIMING=0.",
      call. = FALSE
    )
  }
  checkpoint_results <- checkpoint_results[
    !duplicated(row_key(checkpoint_results), fromLast = TRUE),
    , drop = FALSE
  ]
  message(
    "Loaded timing checkpoint with ", nrow(checkpoint_results),
    " completed analysis rows."
  )
}

completed_keys <- if (is.null(checkpoint_results)) {
  character()
} else {
  row_key(checkpoint_results)
}

for (design_index in seq_len(nrow(benchmark_designs))) {
  design <- benchmark_designs[design_index, , drop = FALSE]
  total_n <- 2L * design$target_per_group
  calibration <- get_cached_calibration(design$factor_count, total_n)

  for (objective in objectives) {
    for (repeat_id in seq_len(N_REPEATS)) {
      candidate_rows <- data.frame(
        outcome = rep(design$outcome, length(analyses)),
        analysis = analyses,
        factor_count = rep(design$factor_count, length(analyses)),
        target_per_group = rep(design$target_per_group, length(analyses)),
        objective = rep(objective, length(analyses)),
        repeat_id = rep(repeat_id, length(analyses)),
        stringsAsFactors = FALSE
      )
      pending_analyses <- analyses[!row_key(candidate_rows) %in% completed_keys]
      if (!length(pending_analyses)) {
        message(
          "Skipping completed timing block: outcome=", design$outcome,
          ", K=", design$factor_count,
          ", n/group=", design$target_per_group,
          ", objective=", objective,
          ", repeat=", repeat_id, "/", N_REPEATS
        )
        next
      }

      message(
        "Preparing trials: outcome=", design$outcome,
        ", K=", design$factor_count,
        ", n/group=", design$target_per_group,
        ", objective=", objective,
        ", repeat=", repeat_id, "/", N_REPEATS
      )
      trials <- prepare_trials(design, objective, repeat_id)

      for (analysis in pending_analyses) {
        message("  Timing ", analysis, " SIGA...")
        siga_time <- time_siga_repeat(
          trials, design$outcome, objective, analysis, calibration
        )
        message("  Timing ", analysis, " RT...")
        rt_time <- time_rt_repeat(
          trials, design, objective, analysis, repeat_id
        )

        new_row <- data.frame(
          outcome = design$outcome,
          analysis = analysis,
          factor_count = design$factor_count,
          target_per_group = design$target_per_group,
          objective = objective,
          repeat_id = repeat_id,
          trials_per_repeat = TRIALS_PER_REPEAT,
          siga_timed_calls = siga_time$calls,
          siga_total_seconds = siga_time$total_seconds,
          siga_seconds_per_trial = siga_time$seconds_per_trial,
          rt_timed_calls = rt_time$calls,
          rt_total_seconds = rt_time$total_seconds,
          rt_seconds_per_trial = rt_time$seconds_per_trial,
          rerandomizations_per_rt = N_RERAND,
          allocation_calibration_paths = N_CALIBRATION,
          stringsAsFactors = FALSE
        )
        checkpoint_results <- if (is.null(checkpoint_results)) {
          new_row
        } else {
          rbind(checkpoint_results, new_row)
        }
        checkpoint_results <- checkpoint_results[
          !duplicated(row_key(checkpoint_results), fromLast = TRUE),
          , drop = FALSE
        ]
        write_checkpoint(checkpoint_results, checkpoint_path)
        completed_keys <- row_key(checkpoint_results)
      }
    }
  }
}

requested_design_keys <- paste(
  benchmark_designs$outcome,
  benchmark_designs$factor_count,
  benchmark_designs$target_per_group,
  sep = "|"
)
checkpoint_design_keys <- paste(
  checkpoint_results$outcome,
  checkpoint_results$factor_count,
  checkpoint_results$target_per_group,
  sep = "|"
)
repeat_results <- checkpoint_results[
  checkpoint_design_keys %in% requested_design_keys &
    checkpoint_results$repeat_id %in% seq_len(N_REPEATS),
  , drop = FALSE
]
repeat_results <- repeat_results[order(
  repeat_results$outcome,
  repeat_results$factor_count,
  repeat_results$target_per_group,
  repeat_results$objective,
  repeat_results$analysis,
  repeat_results$repeat_id
), , drop = FALSE]

repeat_path <- file.path(OUTPUT_DIR, "clean_timing_repeat_means_rerun.csv")
write.csv(repeat_results, repeat_path, row.names = FALSE)

# -----------------------------------------------------------------------------
# Summarize using the exact column structure of the existing clean summary.
# -----------------------------------------------------------------------------

split_key <- interaction(
  repeat_results$outcome,
  repeat_results$analysis,
  repeat_results$factor_count,
  repeat_results$target_per_group,
  repeat_results$objective,
  drop = TRUE, lex.order = TRUE
)
groups <- split(repeat_results, split_key)
summary_rows <- lapply(groups, function(dat) {
  siga_mean <- mean(dat$siga_seconds_per_trial)
  siga_sd <- sd(dat$siga_seconds_per_trial)
  rt_mean <- mean(dat$rt_seconds_per_trial)
  rt_sd <- sd(dat$rt_seconds_per_trial)
  siga_cv <- if (siga_mean > 0) siga_sd / siga_mean else NA_real_
  rt_cv <- if (rt_mean > 0) rt_sd / rt_mean else NA_real_
  data.frame(
    outcome = dat$outcome[1L],
    analysis = dat$analysis[1L],
    factor_count = dat$factor_count[1L],
    target_per_group = dat$target_per_group[1L],
    objective = dat$objective[1L],
    benchmark_repeats = nrow(dat),
    trials_per_repeat = unique(dat$trials_per_repeat)[1L],
    siga_seconds_per_trial = siga_mean,
    siga_repeat_sd_seconds = siga_sd,
    rt_seconds_per_trial = rt_mean,
    rt_repeat_sd_seconds = rt_sd,
    siga_projected_minutes_100k = siga_mean * 100000 / 60,
    rt_projected_minutes_100k = rt_mean * 100000 / 60,
    siga_repeat_cv = siga_cv,
    rt_repeat_cv = rt_cv,
    timing_quality = if (
      is.finite(siga_cv) && is.finite(rt_cv) &&
        siga_cv <= 0.10 && rt_cv <= 0.10
    ) "Stable" else "Review variability",
    stringsAsFactors = FALSE
  )
})
summary_results <- do.call(rbind, summary_rows)
summary_results <- summary_results[order(
  summary_results$outcome,
  summary_results$analysis,
  summary_results$factor_count,
  summary_results$target_per_group,
  summary_results$objective
), , drop = FALSE]
rownames(summary_results) <- NULL

summary_path <- file.path(OUTPUT_DIR, "clean_timing_benchmark_summary_rerun.csv")
write.csv(summary_results, summary_path, row.names = FALSE)

# Merge with the previous summary. Under the default "needed" scope, replace
# every continuous row and add/replace binary 100- and 200-per-group rows while
# preserving the prior binary 500- and 1000-per-group rows.
complete_results <- summary_results
if (file.exists(EXISTING_TIMING_CSV)) {
  previous <- read.csv(EXISTING_TIMING_CSV, stringsAsFactors = FALSE, check.names = FALSE)
  required_columns <- names(summary_results)
  if (!all(required_columns %in% names(previous))) {
    warning(
      "Existing timing CSV does not have the expected columns; ",
      "the complete merged file will contain rerun rows only."
    )
  } else {
    previous <- previous[, required_columns, drop = FALSE]
    key_function <- function(dat) paste(
      dat$outcome, dat$analysis, dat$factor_count,
      dat$target_per_group, dat$objective, sep = "|"
    )
    replacement_keys <- key_function(summary_results)
    previous <- previous[!key_function(previous) %in% replacement_keys, , drop = FALSE]
    complete_results <- rbind(previous, summary_results)
  }
}
complete_results <- complete_results[order(
  complete_results$outcome,
  complete_results$factor_count,
  complete_results$target_per_group,
  complete_results$objective,
  complete_results$analysis
), , drop = FALSE]
rownames(complete_results) <- NULL

complete_path <- file.path(
  OUTPUT_DIR, "clean_timing_benchmark_summary_split_complete.csv"
)
write.csv(complete_results, complete_path, row.names = FALSE)

# -----------------------------------------------------------------------------
# Separately time the one-time, allocation-only calibration for each unique
# design. The same calibration is reusable across continuous and binary
# outcomes, objectives, treatment effects, and adjusted/unadjusted analyses
# whenever the factor distribution, allocation rule, and total sample size are
# unchanged.
# -----------------------------------------------------------------------------

calibration_summary <- NULL
if (TIME_CALIBRATION) {
  calibration_designs <- unique(benchmark_designs[, c(
    "factor_count", "target_per_group"
  ), drop = FALSE])
  calibration_designs <- calibration_designs[order(
    calibration_designs$factor_count,
    calibration_designs$target_per_group
  ), , drop = FALSE]
  rownames(calibration_designs) <- NULL

  calibration_rows <- list()
  calibration_position <- 0L
  for (design_index in seq_len(nrow(calibration_designs))) {
    factor_count <- calibration_designs$factor_count[design_index]
    target_per_group <- calibration_designs$target_per_group[design_index]
    total_n <- 2L * target_per_group

    warm_up_calibration(factor_count, total_n)
    for (repeat_id in seq_len(CALIBRATION_REPEATS)) {
      message(
        "Timing fresh allocation calibration: K=", factor_count,
        ", n/group=", target_per_group,
        ", repeat=", repeat_id, "/", CALIBRATION_REPEATS
      )
      elapsed <- time_fresh_calibration(
        factor_count = factor_count,
        total_n = total_n,
        repeat_id = repeat_id
      )
      calibration_position <- calibration_position + 1L
      calibration_rows[[calibration_position]] <- data.frame(
        factor_count = factor_count,
        target_per_group = target_per_group,
        total_n = total_n,
        calibration_repeat = repeat_id,
        allocation_calibration_paths = N_CALIBRATION,
        calibration_seconds = elapsed,
        calibration_minutes = elapsed / 60,
        stringsAsFactors = FALSE
      )
    }
  }

  calibration_repeats <- do.call(rbind, calibration_rows)
  calibration_repeat_path <- file.path(
    OUTPUT_DIR, "allocation_calibration_repeat_times.csv"
  )
  write.csv(calibration_repeats, calibration_repeat_path, row.names = FALSE)

  calibration_groups <- split(
    calibration_repeats,
    interaction(
      calibration_repeats$factor_count,
      calibration_repeats$target_per_group,
      drop = TRUE, lex.order = TRUE
    )
  )
  calibration_summary <- do.call(rbind, lapply(calibration_groups, function(dat) {
    mean_seconds <- mean(dat$calibration_seconds)
    sd_seconds <- if (nrow(dat) > 1L) sd(dat$calibration_seconds) else NA_real_
    data.frame(
      factor_count = dat$factor_count[1L],
      target_per_group = dat$target_per_group[1L],
      total_n = dat$total_n[1L],
      calibration_repeats = nrow(dat),
      allocation_calibration_paths = dat$allocation_calibration_paths[1L],
      calibration_mean_seconds = mean_seconds,
      calibration_sd_seconds = sd_seconds,
      calibration_mean_minutes = mean_seconds / 60,
      calibration_sd_minutes = sd_seconds / 60,
      stringsAsFactors = FALSE
    )
  }))
  calibration_summary <- calibration_summary[order(
    calibration_summary$factor_count,
    calibration_summary$target_per_group
  ), , drop = FALSE]
  rownames(calibration_summary) <- NULL

  calibration_summary_path <- file.path(
    OUTPUT_DIR, "allocation_calibration_timing_summary.csv"
  )
  write.csv(calibration_summary, calibration_summary_path, row.names = FALSE)

  calibration_tex_path <- file.path(
    OUTPUT_DIR, "allocation_calibration_timing_table.tex"
  )
  con <- file(calibration_tex_path, open = "wt")
  writeLines(c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\small",
    "\\begin{tabular}{rrrrr}",
    "\\toprule",
    "Factors & $n$/group & $B_0$ & Mean calibration time (min) & SD (min)\\\\",
    "\\midrule"
  ), con)
  for (i in seq_len(nrow(calibration_summary))) {
    row <- calibration_summary[i, ]
    sd_text <- if (is.finite(row$calibration_sd_minutes)) {
      sprintf("%.2f", row$calibration_sd_minutes)
    } else {
      "--"
    }
    writeLines(sprintf(
      "%d & %d & %s & %.2f & %s\\\\",
      row$factor_count,
      row$target_per_group,
      format(row$allocation_calibration_paths, big.mark = ",", scientific = FALSE),
      row$calibration_mean_minutes,
      sd_text
    ), con)
  }
  writeLines(c(
    "\\bottomrule",
    "\\end{tabular}",
    paste0(
      "\\caption{One-time single-core computation time for the allocation-only SIGA calibration. ",
      "The calibration uses $B_0=", format(N_CALIBRATION, big.mark = ",", scientific = FALSE),
      "$ simulated allocation paths and is reusable across outcomes, objectives, null boundaries, ",
      "and treatment-effect assumptions when the minimization design and total sample size are unchanged.}"
    ),
    "\\label{tab:allocation-calibration-time}",
    "\\end{table}"
  ), con)
  close(con)

  timing_with_calibration <- merge(
    complete_results,
    calibration_summary[, c(
      "factor_count", "target_per_group",
      "calibration_mean_seconds", "calibration_sd_seconds",
      "calibration_mean_minutes", "calibration_sd_minutes"
    ), drop = FALSE],
    by = c("factor_count", "target_per_group"),
    all.x = TRUE,
    sort = FALSE
  )
  timing_with_calibration$siga_total_minutes_100k_including_one_time_calibration <-
    timing_with_calibration$siga_projected_minutes_100k +
    timing_with_calibration$calibration_mean_minutes
  timing_with_calibration <- timing_with_calibration[order(
    timing_with_calibration$outcome,
    timing_with_calibration$factor_count,
    timing_with_calibration$target_per_group,
    timing_with_calibration$objective,
    timing_with_calibration$analysis
  ), , drop = FALSE]
  rownames(timing_with_calibration) <- NULL
  timing_with_calibration_path <- file.path(
    OUTPUT_DIR, "clean_timing_benchmark_with_calibration_complete.csv"
  )
  write.csv(timing_with_calibration, timing_with_calibration_path, row.names = FALSE)

  message("Wrote calibration repeat times: ", calibration_repeat_path)
  message("Wrote calibration timing summary: ", calibration_summary_path)
  message("Wrote calibration timing LaTeX table: ", calibration_tex_path)
  message("Wrote timing summary including one-time calibration: ", timing_with_calibration_path)
}

message("Wrote repeat-level timing results: ", repeat_path)
message("Wrote rerun timing summary: ", summary_path)
message("Wrote merged complete timing summary: ", complete_path)
message("Rows in merged summary: ", nrow(complete_results))
if (TIMING_SCOPE == "needed" || TIMING_SCOPE == "all") {
  message("Expected rows after adding binary n/group 100 and 200: 48")
}
