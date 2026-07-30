#!/usr/bin/env Rscript

# Additional simulation for the revised SIGA theory.
#
# The program separates:
#   SIGA-S: the trace-matched sampling-calibrated variance;
#   SIGA-R: SIGA-S plus the pair-path variance correction; and
#   RT:     the fixed-score conditional randomization test.
#
# It includes exact-alignment scenarios, the manuscript-style common-log-odds
# binary model, and generalized-eigenvector stress scenarios selected from the
# allocation-only pair-path covariance. The stress scenarios are deliberately
# designed to make the sampling/randomization variance distinction visible.
#
# Base R only. See the accompanying README for smoke, pilot and manuscript runs.

options(stringsAsFactors = FALSE, warn = 1)

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

BASE_DIR <- script_directory()
source(file.path(BASE_DIR, "siga_pair_path_engine.R"), local = FALSE)

ensure_directory <- function(path, attempts = 5L, wait_seconds = 1) {
  path <- path.expand(path)
  for (attempt in seq_len(attempts)) {
    if (dir.exists(path)) return(normalizePath(path, winslash = "/"))
    if (dir.create(path, recursive = TRUE, showWarnings = FALSE) && dir.exists(path)) {
      return(normalizePath(path, winslash = "/"))
    }
    if (attempt < attempts) Sys.sleep(wait_seconds)
  }
  stop("Could not create or access directory: ", path, call. = FALSE)
}

env_integer <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (is.na(value)) stop(name, " must be an integer.", call. = FALSE)
  value
}

env_numeric <- function(name, default) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  if (!is.finite(value)) stop(name, " must be finite.", call. = FALSE)
  value
}

env_integer_vector <- function(name) {
  text <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(text)) return(integer())
  value <- suppressWarnings(as.integer(trimws(strsplit(text, ",", fixed = TRUE)[[1L]])))
  if (anyNA(value)) stop(name, " must contain comma-separated integers.", call. = FALSE)
  unique(value)
}

append_csv <- function(data, path, attempts = 5L, wait_seconds = 1) {
  if (!nrow(data)) return(invisible(NULL))
  parent <- ensure_directory(dirname(path))
  last_error <- NULL

  for (attempt in seq_len(attempts)) {
    tmp <- NULL
    ok <- tryCatch({
      exists <- file.exists(path)
      tmp <- tempfile(
        pattern = paste0(basename(path), ".batch_"),
        tmpdir = parent,
        fileext = ".tmp"
      )
      write.table(
        data,
        file = tmp,
        sep = ",",
        row.names = FALSE,
        col.names = !exists,
        append = FALSE,
        quote = TRUE,
        qmethod = "double"
      )
      if (exists) {
        if (!file.append(path, tmp)) stop("file.append returned FALSE")
        unlink(tmp)
      } else if (!file.rename(tmp, path)) {
        if (!file.copy(tmp, path, overwrite = FALSE)) {
          stop("Could not move the first checkpoint batch into place")
        }
        unlink(tmp)
      }
      TRUE
    }, error = function(e) {
      last_error <<- conditionMessage(e)
      if (!is.null(tmp) && file.exists(tmp)) unlink(tmp)
      FALSE
    })
    if (isTRUE(ok)) return(invisible(NULL))
    if (attempt < attempts) Sys.sleep(wait_seconds)
  }
  stop(
    "Failed to append checkpoint: ", path,
    if (!is.null(last_error)) paste0("; ", last_error) else "",
    call. = FALSE
  )
}

safe_write_csv <- function(data, path) {
  parent <- ensure_directory(dirname(path))
  tmp <- tempfile(pattern = basename(path), tmpdir = parent, fileext = ".tmp")
  write.csv(data, tmp, row.names = FALSE, na = "")
  if (!file.rename(tmp, path)) {
    if (!file.copy(tmp, path, overwrite = TRUE)) {
      unlink(tmp)
      stop("Could not write ", path, call. = FALSE)
    }
    unlink(tmp)
  }
  invisible(path)
}

format_elapsed <- function(seconds) {
  if (!is.finite(seconds)) return("NA")
  if (seconds < 60) return(sprintf("%.1f seconds", seconds))
  if (seconds < 3600) return(sprintf("%.1f minutes", seconds / 60))
  sprintf("%.2f hours", seconds / 3600)
}

PROFILE <- tolower(trimws(Sys.getenv("PWRT_PROFILE", unset = "pilot")))
if (!PROFILE %in% c("smoke", "pilot", "manuscript")) {
  stop("PWRT_PROFILE must be smoke, pilot, or manuscript.", call. = FALSE)
}

profile_defaults <- switch(
  PROFILE,
  smoke = c(outer = 20L, rerand = 99L, calibration = 500L, batch = 10L),
  pilot = c(outer = 5000L, rerand = 999L, calibration = 20000L, batch = 20L),
  manuscript = c(outer = 100000L, rerand = 4999L, calibration = 100000L, batch = 10L)
)

RUN_MODE <- tolower(trimws(Sys.getenv("PWRT_MODE", unset = "run")))
if (!RUN_MODE %in% c("calibrate", "run", "aggregate")) {
  stop("PWRT_MODE must be calibrate, run, or aggregate.", call. = FALSE)
}

SCENARIO_SET <- tolower(trimws(Sys.getenv("PWRT_SCENARIO_SET", unset = "core")))
if (!SCENARIO_SET %in% c("core", "sensitivity", "all")) {
  stop("PWRT_SCENARIO_SET must be core, sensitivity, or all.", call. = FALSE)
}

N_OUTER <- env_integer("PWRT_N_OUTER", profile_defaults["outer"])
N_RERANDOMIZATIONS <- env_integer("PWRT_N_RERAND", profile_defaults["rerand"])
N_CALIBRATION <- env_integer("PWRT_N_CALIBRATION", profile_defaults["calibration"])
CALIBRATION_BATCH <- env_integer("PWRT_CALIBRATION_BATCH", 1000L)
OUTER_BATCH <- env_integer("PWRT_OUTER_BATCH", profile_defaults["batch"])
N_SHARDS <- env_integer("PWRT_N_SHARDS", 1L)
SHARD_ID <- env_integer("PWRT_SHARD_ID", 1L)
BASE_SEED <- env_integer("PWRT_SEED", 20260729L)
SCENARIO_FILTER <- env_integer_vector("PWRT_SCENARIO_IDS")
ALLOW_PARTIAL <- identical(Sys.getenv("PWRT_ALLOW_PARTIAL", unset = "0"), "1")

physical_cores <- parallel::detectCores(logical = FALSE)
if (is.na(physical_cores)) physical_cores <- parallel::detectCores(logical = TRUE)
if (is.na(physical_cores)) physical_cores <- 1L
DEFAULT_CORES <- if (.Platform$OS.type == "windows") {
  1L
} else {
  max(1L, min(8L, physical_cores - 1L))
}
N_CORES <- env_integer("PWRT_CORES", DEFAULT_CORES)
if (.Platform$OS.type == "windows") N_CORES <- 1L
USE_PARALLEL <- .Platform$OS.type == "unix" && N_CORES > 1L

PROJECT_DIR <- ensure_directory(Sys.getenv(
  "PWRT_PROJECT_DIR",
  unset = file.path(dirname(BASE_DIR), "pair_path_theory_output")
))
OUTPUT_ROOT <- ensure_directory(Sys.getenv(
  "PWRT_OUTPUT_DIR",
  unset = file.path(PROJECT_DIR, "simulation_output")
))
RUN_TAG <- paste0(
  "profile_", PROFILE,
  "_M", N_OUTER,
  "_B", N_RERANDOMIZATIONS,
  "_B0", N_CALIBRATION,
  "_S", N_SHARDS
)
OUTPUT_DIR <- ensure_directory(file.path(OUTPUT_ROOT, RUN_TAG))
CALIBRATION_DIR <- ensure_directory(file.path(PROJECT_DIR, "calibration_cache"))
SCENARIO_DIR <- ensure_directory(file.path(OUTPUT_DIR, "scenario_shards"))

stopifnot(
  N_OUTER >= 1L,
  N_RERANDOMIZATIONS >= 1L,
  N_CALIBRATION >= 2L,
  CALIBRATION_BATCH >= 1L,
  OUTER_BATCH >= 1L,
  N_SHARDS >= 1L,
  SHARD_ID >= 1L,
  SHARD_ID <= N_SHARDS,
  N_CORES >= 1L
)

if (PROFILE == "manuscript") {
  expected <- c(outer = 100000L, rerand = 4999L, calibration = 100000L)
  observed <- c(outer = N_OUTER, rerand = N_RERANDOMIZATIONS, calibration = N_CALIBRATION)
  if (any(expected != observed)) {
    warning(
      "Manuscript profile was selected but settings differ from the recommended values: ",
      paste(names(observed), observed, sep = "=", collapse = ", ")
    )
  }
}

seed_value <- function(scenario_id, replicate_id, stream_id) {
  modulus <- 2147483646
  value <- as.double(BASE_SEED) +
    1000003 * as.double(stream_id) +
    104729 * as.double(scenario_id) +
    1009 * as.double(replicate_id)
  as.integer(value %% modulus + 1)
}

calibration_seed <- function(design_id) {
  modulus <- 2147483646
  value <- as.double(BASE_SEED) + 7000003 + 99991 * as.double(design_id)
  as.integer(value %% modulus + 1)
}

make_design_table <- function() {
  data.frame(
    design_id = 1:4,
    design_label = c(
      "K2_n200_p080_independent",
      "K5_n400_p080_independent",
      "K5_n400_p080_correlated",
      "K2_n200_p095_independent"
    ),
    factor_count = c(2L, 5L, 5L, 2L),
    total_n = c(200L, 400L, 400L, 200L),
    pbc = c(0.80, 0.80, 0.80, 0.95),
    profile_type = c("independent", "independent", "latent_correlated", "independent"),
    latent_strength = c(0, 0, 1.0, 0),
    set = c("core", "core", "sensitivity", "sensitivity"),
    stringsAsFactors = FALSE
  )
}

FACTOR_PREVALENCE_MASTER <- c(0.50, 0.40, 0.30, 0.20, 0.10)

prepare_design <- function(row) {
  K <- as.integer(row$factor_count)
  patterns <- all_binary_patterns(K)
  factor_prob <- FACTOR_PREVALENCE_MASTER[seq_len(K)]
  profile_prob <- if (row$profile_type == "independent") {
    profile_prob_independent(factor_prob, patterns)
  } else {
    profile_prob_latent_correlated(
      factor_prob = factor_prob,
      latent_strength = row$latent_strength,
      mixing_prob = 0.5,
      patterns = patterns
    )
  }
  list(
    design_id = as.integer(row$design_id),
    design_label = as.character(row$design_label),
    factor_count = K,
    total_n = as.integer(row$total_n),
    pbc = as.numeric(row$pbc),
    weights = rep(1, K + 1L),
    profile_type = as.character(row$profile_type),
    latent_strength = as.numeric(row$latent_strength),
    patterns = patterns,
    factor_prob = factor_prob,
    profile_prob = profile_prob,
    marginal_prob = profile_marginals(patterns, profile_prob),
    correlation = profile_correlations(patterns, profile_prob),
    set = as.character(row$set)
  )
}

make_scenario_blueprints <- function() {
  rows <- list(
    c(1, "continuous", "homogeneous", "zero", "core"),
    c(1, "continuous", "pure_symmetric", "max_ratio", "core"),
    c(1, "continuous", "realistic", "max_ratio", "core"),
    c(1, "binary", "homogeneous", "zero", "core"),
    c(1, "binary", "symmetric_stress", "max_ratio", "core"),
    c(1, "binary", "common_log_odds", "common_log_odds", "core"),
    c(2, "continuous", "homogeneous", "zero", "core"),
    c(2, "continuous", "pure_symmetric", "min_ratio", "core"),
    c(2, "continuous", "pure_symmetric", "max_ratio", "core"),
    c(2, "continuous", "realistic", "max_ratio", "core"),
    c(2, "binary", "homogeneous", "zero", "core"),
    c(2, "binary", "symmetric_stress", "min_ratio", "core"),
    c(2, "binary", "symmetric_stress", "max_ratio", "core"),
    c(2, "binary", "common_log_odds", "common_log_odds", "core"),
    c(3, "continuous", "pure_symmetric", "min_ratio", "sensitivity"),
    c(3, "continuous", "pure_symmetric", "max_ratio", "sensitivity"),
    c(3, "binary", "symmetric_stress", "min_ratio", "sensitivity"),
    c(3, "binary", "symmetric_stress", "max_ratio", "sensitivity"),
    c(4, "continuous", "pure_symmetric", "max_ratio", "sensitivity"),
    c(4, "binary", "symmetric_stress", "max_ratio", "sensitivity")
  )
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- c("design_id", "outcome_type", "model_type", "direction", "set")
  out$design_id <- as.integer(out$design_id)
  out$scenario_id <- seq_len(nrow(out))
  out$boundary <- ifelse(out$outcome_type == "continuous", -0.20, -0.10)
  out$alpha <- 0.025
  out$alternative <- "greater"
  out$scenario_label <- paste0(
    "D", out$design_id, "_", out$outcome_type, "_", out$model_type, "_", out$direction
  )
  out
}

DESIGN_TABLE <- make_design_table()
DESIGNS <- setNames(
  lapply(seq_len(nrow(DESIGN_TABLE)), function(i) prepare_design(DESIGN_TABLE[i, ])),
  DESIGN_TABLE$design_id
)

SCENARIOS <- make_scenario_blueprints()
if (SCENARIO_SET != "all") {
  SCENARIOS <- SCENARIOS[SCENARIOS$set == SCENARIO_SET, , drop = FALSE]
}
if (length(SCENARIO_FILTER)) {
  SCENARIOS <- SCENARIOS[SCENARIOS$scenario_id %in% SCENARIO_FILTER, , drop = FALSE]
}
if (!nrow(SCENARIOS)) stop("No scenarios remain after filtering.", call. = FALSE)
rownames(SCENARIOS) <- NULL

calibration_cache_file <- function(design) {
  file.path(
    CALIBRATION_DIR,
    paste0(
      "pair_calibration_", design$design_label,
      "_B0", N_CALIBRATION,
      "_seed", calibration_seed(design$design_id),
      ".rds"
    )
  )
}

get_design_calibration <- function(design) {
  cache <- calibration_cache_file(design)
  if (file.exists(cache)) return(readRDS(cache))

  lock <- paste0(cache, ".lock")
  have_lock <- dir.create(lock, showWarnings = FALSE)
  if (!have_lock) {
    for (attempt in seq_len(720L)) {
      if (file.exists(cache)) return(readRDS(cache))
      Sys.sleep(5)
      have_lock <- dir.create(lock, showWarnings = FALSE)
      if (have_lock) break
    }
  }
  if (!have_lock) stop("Could not acquire calibration lock: ", lock, call. = FALSE)
  on.exit(unlink(lock, recursive = TRUE, force = TRUE), add = TRUE)
  if (file.exists(cache)) return(readRDS(cache))

  message(
    "Creating three-copy allocation calibration for ", design$design_label,
    ": n=", design$total_n, ", B0=", N_CALIBRATION
  )
  calibration <- calibrate_pair_path_design_R(
    B0 = N_CALIBRATION,
    n = design$total_n,
    patterns = design$patterns,
    profile_prob = design$profile_prob,
    pbc = design$pbc,
    weights_ = design$weights,
    seed = calibration_seed(design$design_id),
    batch_size = CALIBRATION_BATCH,
    progress = TRUE
  )
  calibration$directions <- generalized_pair_directions(
    calibration$psi,
    calibration$profile_prob
  )

  tmp <- tempfile(pattern = "pair_calibration_", tmpdir = CALIBRATION_DIR, fileext = ".rds")
  saveRDS(calibration, tmp)
  if (!file.rename(tmp, cache)) {
    if (!file.copy(tmp, cache, overwrite = TRUE)) {
      unlink(tmp)
      stop("Could not save calibration: ", cache, call. = FALSE)
    }
    unlink(tmp)
  }
  calibration
}

USED_DESIGN_IDS <- sort(unique(SCENARIOS$design_id))
CALIBRATIONS <- setNames(
  lapply(USED_DESIGN_IDS, function(id) get_design_calibration(DESIGNS[[as.character(id)]])),
  USED_DESIGN_IDS
)

direction_for_scenario <- function(scenario, calibration) {
  switch(
    scenario$direction,
    zero = rep(0, calibration$J),
    min_ratio = calibration$directions$min_direction,
    max_ratio = calibration$directions$max_direction,
    common_log_odds = rep(NA_real_, calibration$J),
    stop("Unknown direction: ", scenario$direction, call. = FALSE)
  )
}

build_scenario_model <- function(scenario) {
  design <- DESIGNS[[as.character(scenario$design_id)]]
  calibration <- CALIBRATIONS[[as.character(scenario$design_id)]]
  direction <- direction_for_scenario(scenario, calibration)
  boundary <- scenario$boundary

  model <- if (scenario$outcome_type == "continuous") {
    switch(
      scenario$model_type,
      homogeneous = make_continuous_homogeneous_model(
        design$patterns, design$profile_prob, boundary, outcome_sd = 1.0
      ),
      pure_symmetric = make_continuous_pure_model(
        design$patterns, design$profile_prob, boundary, direction,
        target_max_abs_d = 1.0, shared_noise_sd = 0.05
      ),
      realistic = make_continuous_realistic_model(
        design$patterns, design$profile_prob, boundary, direction,
        target_max_abs_d = 0.50, outcome_sd = 1.0, individual_effect_sd = 0.25
      ),
      stop("Unknown continuous model type: ", scenario$model_type, call. = FALSE)
    )
  } else {
    switch(
      scenario$model_type,
      homogeneous = make_binary_homogeneous_model(
        design$patterns, design$profile_prob, boundary, target_control_risk = 0.60
      ),
      symmetric_stress = make_binary_symmetric_model(
        design$patterns, design$profile_prob, boundary, direction,
        target_max_abs_d = 0.80, max_abs_delta = 0.90
      ),
      common_log_odds = calibrate_common_log_odds_model(
        design$patterns, design$profile_prob,
        target_control_risk = 0.60,
        target_risk_difference = boundary
      ),
      stop("Unknown binary model type: ", scenario$model_type, call. = FALSE)
    )
  }

  d <- as.numeric(model$d)
  weighted_mean_d <- sum(design$profile_prob * d)
  if (abs(weighted_mean_d) > 1e-9) {
    stop(
      "Scenario ", scenario$scenario_id,
      " does not satisfy the marginal boundary: weighted mean d=", weighted_mean_d,
      call. = FALSE
    )
  }
  d_pi <- sum(design$profile_prob * d^2)
  d_psi <- drop(crossprod(d, calibration$psi %*% d))
  ratio <- if (d_pi > 0) d_psi / d_pi else 1

  model$d_pi <- d_pi
  model$d_psi <- d_psi
  model$pair_ratio <- ratio
  model$pair_scalar_gap <- d_psi - d_pi
  model
}

SCENARIO_MODELS <- setNames(
  lapply(seq_len(nrow(SCENARIOS)), function(i) build_scenario_model(SCENARIOS[i, ])),
  SCENARIOS$scenario_id
)

build_design_diagnostics <- function() {
  rows <- lapply(USED_DESIGN_IDS, function(id) {
    design <- DESIGNS[[as.character(id)]]
    calibration <- CALIBRATIONS[[as.character(id)]]
    corr <- design$correlation
    offdiag <- corr[row(corr) != col(corr)]
    data.frame(
      design_id = design$design_id,
      design_label = design$design_label,
      factor_count = design$factor_count,
      total_n = design$total_n,
      pbc = design$pbc,
      profile_type = design$profile_type,
      latent_strength = design$latent_strength,
      calibration_paths = calibration$B0,
      calibration_seconds = calibration$elapsed_seconds,
      min_generalized_ratio = calibration$directions$min_ratio,
      max_generalized_ratio = calibration$directions$max_ratio,
      max_abs_mean_u = max(abs(calibration$mean_u)),
      max_abs_mean_g01 = max(abs(calibration$mean_g01)),
      max_abs_mean_g02 = max(abs(calibration$mean_g02)),
      max_abs_cross_D01 = matrix_max_abs(calibration$cross_d01),
      max_abs_cross_D02 = matrix_max_abs(calibration$cross_d02),
      max_abs_cross_D12 = matrix_max_abs(calibration$cross_d12),
      max_abs_cross_G0102 = matrix_max_abs(calibration$cross_g0102),
      frobenius_cross_G0102 = matrix_frobenius(calibration$cross_g0102),
      max_abs_factor_correlation = if (length(offdiag)) max(abs(offdiag)) else 0,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

build_scenario_definition_table <- function() {
  rows <- lapply(seq_len(nrow(SCENARIOS)), function(i) {
    scenario <- SCENARIOS[i, ]
    design <- DESIGNS[[as.character(scenario$design_id)]]
    model <- SCENARIO_MODELS[[as.character(scenario$scenario_id)]]
    data.frame(
      scenario_id = scenario$scenario_id,
      scenario_label = scenario$scenario_label,
      set = scenario$set,
      design_id = scenario$design_id,
      design_label = design$design_label,
      factor_count = design$factor_count,
      total_n = design$total_n,
      pbc = design$pbc,
      profile_type = design$profile_type,
      outcome_type = scenario$outcome_type,
      model_type = scenario$model_type,
      direction = scenario$direction,
      boundary = scenario$boundary,
      alpha = scenario$alpha,
      achieved_effect = model$achieved_effect,
      min_d = min(model$d),
      max_d = max(model$d),
      weighted_mean_d = sum(design$profile_prob * model$d),
      d_Pi_d = model$d_pi,
      d_Psi_d = model$d_psi,
      pair_ratio = model$pair_ratio,
      pair_scalar_gap = model$pair_scalar_gap,
      achieved_control_risk = model$achieved_control_risk %||% NA_real_,
      achieved_treatment_risk = model$achieved_treatment_risk %||% NA_real_,
      common_log_odds_theta = model$theta %||% NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

DESIGN_DIAGNOSTICS <- build_design_diagnostics()
SCENARIO_DEFINITIONS <- build_scenario_definition_table()
safe_write_csv(DESIGN_DIAGNOSTICS, file.path(OUTPUT_DIR, "design_pair_path_diagnostics.csv"))
safe_write_csv(SCENARIO_DEFINITIONS, file.path(OUTPUT_DIR, "scenario_definitions.csv"))

if (RUN_MODE == "calibrate") {
  message("Calibration and scenario-definition files were written to: ", OUTPUT_DIR)
  quit(save = "no", status = 0L)
}

simulate_trial <- function(scenario, replicate_id, design, model) {
  set.seed(seed_value(scenario$scenario_id, replicate_id, 1L))
  profile <- generate_profile_sequence(
    n = design$total_n,
    patterns = design$patterns,
    profile_prob = design$profile_prob
  )
  z <- ps_assign_R(
    X = profile$X,
    pbc = design$pbc,
    weights_ = design$weights,
    seed = seed_value(scenario$scenario_id, replicate_id, 2L)
  )
  A <- as.integer((z + 1L) / 2L)
  set.seed(seed_value(scenario$scenario_id, replicate_id, 3L))
  y <- generate_outcome_from_model(profile$id, A, z, model)
  list(
    profile_id = profile$id,
    X = profile$X,
    z = as.integer(z),
    A = A,
    y = y
  )
}

analyze_one_score <- function(score,
                              trial,
                              scenario,
                              calibration,
                              d) {
  statistic <- 0.5 * sum(trial$z * score)
  sampling <- siga_sampling_variance(score, trial$X, calibration)
  randomization <- siga_randomization_variance(sampling, d, calibration)
  p_s <- gaussian_pvalue(
    statistic,
    sampling$variance,
    alternative = scenario$alternative
  )
  p_r <- gaussian_pvalue(
    statistic,
    randomization$variance,
    alternative = scenario$alternative
  )
  list(
    statistic = statistic,
    sampling_variance = sampling$variance,
    randomization_variance = randomization$variance,
    randomization_variance_raw = randomization$variance_raw,
    pair_correction = randomization$correction,
    pair_scalar_gap_realized = randomization$scalar_gap,
    p_s = p_s$p,
    z_s = p_s$z,
    p_r = p_r$p,
    z_r = p_r$z,
    kappa = sampling$kappa,
    kappa_raw = sampling$kappa_raw,
    kappa_truncated = sampling$kappa_raw < 0,
    randomization_variance_truncated = randomization$truncated,
    low_component = sampling$low_component,
    within_component = sampling$within_component,
    counts = sampling$decomposition$counts
  )
}

run_replicate <- function(scenario, replicate_id, calibration, design, model) {
  generation_start <- proc.time()[3L]
  trial <- simulate_trial(scenario, replicate_id, design, model)
  generation_seconds <- proc.time()[3L] - generation_start

  score_start <- proc.time()[3L]
  scores <- make_boundary_scores(
    y = trial$y,
    A = trial$A,
    X = trial$X,
    boundary = scenario$boundary
  )
  score_seconds <- proc.time()[3L] - score_start

  siga_start <- proc.time()[3L]
  unadjusted <- analyze_one_score(
    score = scores$unadjusted,
    trial = trial,
    scenario = scenario,
    calibration = calibration,
    d = model$d
  )
  adjusted <- analyze_one_score(
    score = scores$adjusted,
    trial = trial,
    scenario = scenario,
    calibration = calibration,
    d = model$d
  )
  siga_seconds <- proc.time()[3L] - siga_start

  rt_start <- proc.time()[3L]
  rt <- rt_pvalues_R(
    X = trial$X,
    scores = cbind(unadjusted = scores$unadjusted, adjusted = scores$adjusted),
    z_obs = trial$z,
    B = N_RERANDOMIZATIONS,
    pbc = design$pbc,
    weights_ = design$weights,
    seed = seed_value(scenario$scenario_id, replicate_id, 4L),
    return_randomization_moments = TRUE
  )
  rt_seconds <- proc.time()[3L] - rt_start

  p_rt_u <- rt$greater[1L]
  p_rt_a <- rt$greater[2L]
  rt_var_u <- rt$randomization_variance[1L]
  rt_var_a <- rt$randomization_variance[2L]
  rt_mean_u <- rt$randomization_mean[1L]
  rt_mean_a <- rt$randomization_mean[2L]

  data.frame(
    scenario_id = scenario$scenario_id,
    scenario_label = scenario$scenario_label,
    replicate = replicate_id,
    set = scenario$set,
    design_id = scenario$design_id,
    design_label = design$design_label,
    factor_count = design$factor_count,
    total_n = design$total_n,
    pbc = design$pbc,
    profile_type = design$profile_type,
    outcome_type = scenario$outcome_type,
    model_type = scenario$model_type,
    direction = scenario$direction,
    boundary = scenario$boundary,
    alpha = scenario$alpha,
    achieved_effect = model$achieved_effect,
    pair_ratio_model = model$pair_ratio,
    pair_scalar_gap_model = model$pair_scalar_gap,
    configured_outer_trials = N_OUTER,
    rerandomizations_per_trial = N_RERANDOMIZATIONS,
    allocation_calibration_paths = N_CALIBRATION,
    configured_shards = N_SHARDS,
    statistic_unadjusted = unadjusted$statistic,
    statistic_adjusted = adjusted$statistic,
    p_siga_s_unadjusted = unadjusted$p_s,
    p_siga_r_unadjusted = unadjusted$p_r,
    p_rt_unadjusted = p_rt_u,
    p_siga_s_adjusted = adjusted$p_s,
    p_siga_r_adjusted = adjusted$p_r,
    p_rt_adjusted = p_rt_a,
    reject_siga_s_unadjusted = unadjusted$p_s <= scenario$alpha,
    reject_siga_r_unadjusted = unadjusted$p_r <= scenario$alpha,
    reject_rt_unadjusted = p_rt_u <= scenario$alpha,
    reject_siga_s_adjusted = adjusted$p_s <= scenario$alpha,
    reject_siga_r_adjusted = adjusted$p_r <= scenario$alpha,
    reject_rt_adjusted = p_rt_a <= scenario$alpha,
    variance_s_unadjusted = unadjusted$sampling_variance,
    variance_r_unadjusted = unadjusted$randomization_variance,
    variance_rt_unadjusted = rt_var_u,
    variance_s_adjusted = adjusted$sampling_variance,
    variance_r_adjusted = adjusted$randomization_variance,
    variance_rt_adjusted = rt_var_a,
    pair_correction_unadjusted = unadjusted$pair_correction,
    pair_correction_adjusted = adjusted$pair_correction,
    pair_scalar_gap_realized = unadjusted$pair_scalar_gap_realized,
    randomization_mean_unadjusted = rt_mean_u,
    randomization_mean_adjusted = rt_mean_a,
    kappa_raw_unadjusted = unadjusted$kappa_raw,
    kappa_raw_adjusted = adjusted$kappa_raw,
    kappa_truncated_unadjusted = unadjusted$kappa_truncated,
    kappa_truncated_adjusted = adjusted$kappa_truncated,
    variance_r_truncated_unadjusted = unadjusted$randomization_variance_truncated,
    variance_r_truncated_adjusted = adjusted$randomization_variance_truncated,
    low_component_unadjusted = unadjusted$low_component,
    within_component_unadjusted = unadjusted$within_component,
    low_component_adjusted = adjusted$low_component,
    within_component_adjusted = adjusted$within_component,
    observed_treated = sum(trial$A),
    data_generation_seconds = generation_seconds,
    score_construction_seconds = score_seconds,
    siga_analysis_seconds = siga_seconds,
    rt_analysis_seconds = rt_seconds,
    stringsAsFactors = FALSE
  )
}

scenario_paths <- function(scenario) {
  directory <- ensure_directory(file.path(
    SCENARIO_DIR,
    sprintf("scenario_%02d_%s", scenario$scenario_id, scenario$scenario_label)
  ))
  list(
    directory = directory,
    shard = file.path(directory, sprintf("shard_%03d_of_%03d.csv", SHARD_ID, N_SHARDS))
  )
}

read_completed_replicates <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) return(integer())
  dat <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(dat) || !"replicate" %in% names(dat)) return(integer())
  unique(as.integer(dat$replicate))
}

run_scenario_shard <- function(scenario) {
  design <- DESIGNS[[as.character(scenario$design_id)]]
  calibration <- CALIBRATIONS[[as.character(scenario$design_id)]]
  model <- SCENARIO_MODELS[[as.character(scenario$scenario_id)]]
  paths <- scenario_paths(scenario)

  assigned <- seq.int(from = SHARD_ID, to = N_OUTER, by = N_SHARDS)
  completed <- read_completed_replicates(paths$shard)
  pending <- setdiff(assigned, completed)

  message(
    "Scenario ", scenario$scenario_id, " [", scenario$scenario_label, "]",
    ": shard ", SHARD_ID, "/", N_SHARDS,
    "; assigned=", length(assigned),
    "; completed=", length(completed),
    "; pending=", length(pending)
  )
  if (!length(pending)) return(invisible(NULL))

  start <- proc.time()[3L]
  for (from in seq.int(1L, length(pending), by = OUTER_BATCH)) {
    ids <- pending[from:min(length(pending), from + OUTER_BATCH - 1L)]
    worker <- function(id) {
      run_replicate(scenario, id, calibration, design, model)
    }
    rows <- if (USE_PARALLEL && length(ids) > 1L) {
      parallel::mclapply(ids, worker, mc.cores = min(N_CORES, length(ids)))
    } else {
      lapply(ids, worker)
    }
    batch <- do.call(rbind, rows)
    append_csv(batch, paths$shard)

    completed_now <- length(completed) + min(length(pending), from + length(ids) - 1L)
    message(
      "  completed ", completed_now, "/", length(assigned),
      " assigned trials; elapsed ", format_elapsed(proc.time()[3L] - start)
    )
  }
  invisible(NULL)
}

read_scenario_details <- function(scenario) {
  pattern <- sprintf("scenario_%02d_%s", scenario$scenario_id, scenario$scenario_label)
  directory <- file.path(SCENARIO_DIR, pattern)
  if (!dir.exists(directory)) stop("Scenario directory is missing: ", directory, call. = FALSE)
  shard_pattern <- sprintf("^shard_[0-9]+_of_%03d\\.csv$", N_SHARDS)
  files <- list.files(directory, pattern = shard_pattern, full.names = TRUE)
  if (!length(files)) stop("No shard files found for ", scenario$scenario_label, call. = FALSE)

  pieces <- lapply(files, function(path) {
    dat <- read.csv(path, stringsAsFactors = FALSE)
    dat$source_file <- basename(path)
    dat
  })
  dat <- do.call(rbind, pieces)
  required_metadata <- c(
    "scenario_id", "configured_outer_trials", "rerandomizations_per_trial",
    "allocation_calibration_paths", "configured_shards"
  )
  missing_metadata <- setdiff(required_metadata, names(dat))
  if (length(missing_metadata)) {
    stop("Missing metadata columns: ", paste(missing_metadata, collapse = ", "), call. = FALSE)
  }
  expected_metadata <- list(
    scenario_id = scenario$scenario_id,
    configured_outer_trials = N_OUTER,
    rerandomizations_per_trial = N_RERANDOMIZATIONS,
    allocation_calibration_paths = N_CALIBRATION,
    configured_shards = N_SHARDS
  )
  for (field in names(expected_metadata)) {
    observed <- unique(as.numeric(dat[[field]]))
    if (length(observed) != 1L || !isTRUE(all.equal(observed, as.numeric(expected_metadata[[field]])))) {
      stop(
        "Incompatible metadata in ", scenario$scenario_label, ": ", field,
        "=", paste(observed, collapse = ","),
        "; expected ", expected_metadata[[field]],
        call. = FALSE
      )
    }
  }
  dat$replicate <- as.integer(dat$replicate)
  dat <- dat[order(dat$replicate), , drop = FALSE]

  duplicate_id <- duplicated(dat$replicate)
  if (any(duplicate_id)) {
    duplicates <- unique(dat$replicate[duplicate_id])
    for (id in duplicates) {
      block <- dat[dat$replicate == id, setdiff(names(dat), "source_file"), drop = FALSE]
      if (nrow(unique(block)) > 1L) {
        stop("Conflicting duplicate replicate ", id, " in ", scenario$scenario_label, call. = FALSE)
      }
    }
    dat <- dat[!duplicated(dat$replicate), , drop = FALSE]
  }

  if (!ALLOW_PARTIAL && nrow(dat) != N_OUTER) {
    stop(
      scenario$scenario_label, " contains ", nrow(dat),
      " unique replicates; expected ", N_OUTER, ".",
      call. = FALSE
    )
  }
  dat
}

rejection_summary <- function(x) {
  x <- as.logical(x)
  x <- x[!is.na(x)]
  if (!length(x)) {
    return(c(probability = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  ci <- wilson_interval(sum(x), length(x))
  c(probability = ci["estimate"], lower = ci["lower"], upper = ci["upper"])
}

safe_ratio <- function(numerator, denominator) {
  out <- numerator / denominator
  out[!is.finite(out)] <- NA_real_
  out
}

summarize_analysis <- function(dat, analysis, alpha) {
  suffix <- if (analysis == "unadjusted") "unadjusted" else "adjusted"
  r_s <- rejection_summary(dat[[paste0("reject_siga_s_", suffix)]])
  r_r <- rejection_summary(dat[[paste0("reject_siga_r_", suffix)]])
  r_rt <- rejection_summary(dat[[paste0("reject_rt_", suffix)]])

  diff_s_rt <- paired_difference_interval(
    dat[[paste0("reject_siga_s_", suffix)]],
    dat[[paste0("reject_rt_", suffix)]]
  )
  diff_r_rt <- paired_difference_interval(
    dat[[paste0("reject_siga_r_", suffix)]],
    dat[[paste0("reject_rt_", suffix)]]
  )

  v_s <- as.numeric(dat[[paste0("variance_s_", suffix)]])
  v_r <- as.numeric(dat[[paste0("variance_r_", suffix)]])
  v_rt <- as.numeric(dat[[paste0("variance_rt_", suffix)]])
  p_s <- as.numeric(dat[[paste0("p_siga_s_", suffix)]])
  p_r <- as.numeric(dat[[paste0("p_siga_r_", suffix)]])
  p_rt <- as.numeric(dat[[paste0("p_rt_", suffix)]])

  ratio_r_s <- safe_ratio(v_r, v_s)
  ratio_s_rt <- safe_ratio(v_s, v_rt)
  ratio_r_rt <- safe_ratio(v_r, v_rt)
  mean_ratio <- mean(ratio_r_s, na.rm = TRUE)
  predicted_rt_size <- if (is.finite(mean_ratio) && mean_ratio > 0) {
    1 - pnorm(qnorm(1 - alpha) * sqrt(mean_ratio))
  } else {
    NA_real_
  }

  data.frame(
    analysis = analysis,
    siga_s = unname(r_s["probability"]),
    siga_s_lower_95 = unname(r_s["lower"]),
    siga_s_upper_95 = unname(r_s["upper"]),
    siga_r = unname(r_r["probability"]),
    siga_r_lower_95 = unname(r_r["lower"]),
    siga_r_upper_95 = unname(r_r["upper"]),
    rt = unname(r_rt["probability"]),
    rt_lower_95 = unname(r_rt["lower"]),
    rt_upper_95 = unname(r_rt["upper"]),
    difference_siga_s_minus_rt = unname(diff_s_rt["estimate"]),
    difference_siga_s_minus_rt_se = unname(diff_s_rt["se"]),
    difference_siga_s_minus_rt_lower_95 = unname(diff_s_rt["lower"]),
    difference_siga_s_minus_rt_upper_95 = unname(diff_s_rt["upper"]),
    difference_siga_r_minus_rt = unname(diff_r_rt["estimate"]),
    difference_siga_r_minus_rt_se = unname(diff_r_rt["se"]),
    difference_siga_r_minus_rt_lower_95 = unname(diff_r_rt["lower"]),
    difference_siga_r_minus_rt_upper_95 = unname(diff_r_rt["upper"]),
    mean_variance_s = mean(v_s, na.rm = TRUE),
    mean_variance_r = mean(v_r, na.rm = TRUE),
    mean_variance_rt = mean(v_rt, na.rm = TRUE),
    mean_variance_ratio_r_over_s = mean(ratio_r_s, na.rm = TRUE),
    mean_variance_ratio_s_over_rt = mean(ratio_s_rt, na.rm = TRUE),
    mean_variance_ratio_r_over_rt = mean(ratio_r_rt, na.rm = TRUE),
    median_variance_ratio_r_over_s = median(ratio_r_s, na.rm = TRUE),
    predicted_rt_size_from_mean_ratio = predicted_rt_size,
    mean_abs_p_difference_siga_s_rt = mean(abs(p_s - p_rt), na.rm = TRUE),
    mean_abs_p_difference_siga_r_rt = mean(abs(p_r - p_rt), na.rm = TRUE),
    p95_abs_p_difference_siga_s_rt = unname(quantile(abs(p_s - p_rt), 0.95, na.rm = TRUE)),
    p95_abs_p_difference_siga_r_rt = unname(quantile(abs(p_r - p_rt), 0.95, na.rm = TRUE)),
    kappa_raw_mean = mean(dat[[paste0("kappa_raw_", suffix)]], na.rm = TRUE),
    kappa_truncation_rate = mean(dat[[paste0("kappa_truncated_", suffix)]], na.rm = TRUE),
    randomization_variance_truncation_rate = mean(
      dat[[paste0("variance_r_truncated_", suffix)]], na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}

summarize_scenario <- function(dat, scenario) {
  design <- DESIGNS[[as.character(scenario$design_id)]]
  model <- SCENARIO_MODELS[[as.character(scenario$scenario_id)]]
  by_analysis <- rbind(
    summarize_analysis(dat, "unadjusted", scenario$alpha),
    summarize_analysis(dat, "adjusted", scenario$alpha)
  )
  fixed <- data.frame(
    scenario_id = scenario$scenario_id,
    scenario_label = scenario$scenario_label,
    set = scenario$set,
    design_id = scenario$design_id,
    design_label = design$design_label,
    factor_count = design$factor_count,
    total_n = design$total_n,
    pbc = design$pbc,
    profile_type = design$profile_type,
    outcome_type = scenario$outcome_type,
    model_type = scenario$model_type,
    direction = scenario$direction,
    boundary = scenario$boundary,
    alpha = scenario$alpha,
    achieved_effect = model$achieved_effect,
    pair_ratio_model = model$pair_ratio,
    pair_scalar_gap_model = model$pair_scalar_gap,
    n_outer = nrow(dat),
    rerandomizations_per_trial = N_RERANDOMIZATIONS,
    allocation_calibration_paths = N_CALIBRATION,
    mean_realized_pair_scalar_gap = mean(dat$pair_scalar_gap_realized, na.rm = TRUE),
    mean_observed_treated = mean(dat$observed_treated, na.rm = TRUE),
    data_generation_minutes = sum(dat$data_generation_seconds, na.rm = TRUE) / 60,
    score_construction_minutes = sum(dat$score_construction_seconds, na.rm = TRUE) / 60,
    siga_analysis_minutes = sum(dat$siga_analysis_seconds, na.rm = TRUE) / 60,
    rt_analysis_minutes = sum(dat$rt_analysis_seconds, na.rm = TRUE) / 60,
    stringsAsFactors = FALSE
  )
  cbind(fixed[rep(1L, nrow(by_analysis)), , drop = FALSE], by_analysis)
}

write_latex_summary <- function(summary, path) {
  fmt <- function(x, digits = 2) {
    ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "--")
  }
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Additional pair-path theory validation. Rejection probabilities are percentages.}",
    "\\label{tab:pair-path-validation}",
    "\\small",
    "\\begin{tabular}{rllllrrrr}",
    "\\toprule",
    "ID & Outcome & Model & Direction & Analysis & $V_R/V_S$ & SIGA-S & SIGA-R & RT \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(summary))) {
    row <- summary[i, ]
    lines <- c(lines, paste0(
      row$scenario_id, " & ",
      row$outcome_type, " & ",
      gsub("_", "-", row$model_type, fixed = TRUE), " & ",
      gsub("_", "-", row$direction, fixed = TRUE), " & ",
      row$analysis, " & ",
      fmt(row$mean_variance_ratio_r_over_s, 3), " & ",
      fmt(100 * row$siga_s, 2), " & ",
      fmt(100 * row$siga_r, 2), " & ",
      fmt(100 * row$rt, 2), " \\\\"
    ))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, path)
}

aggregate_results <- function() {
  all_summary <- list()
  all_details <- list()
  for (i in seq_len(nrow(SCENARIOS))) {
    scenario <- SCENARIOS[i, ]
    message("Aggregating ", scenario$scenario_label)
    dat <- read_scenario_details(scenario)
    all_details[[i]] <- dat
    all_summary[[i]] <- summarize_scenario(dat, scenario)
  }
  summary <- do.call(rbind, all_summary)
  details <- do.call(rbind, all_details)
  summary <- summary[order(summary$scenario_id, summary$analysis), , drop = FALSE]
  details <- details[order(details$scenario_id, details$replicate), , drop = FALSE]

  safe_write_csv(summary, file.path(OUTPUT_DIR, "pair_path_theory_validation_summary.csv"))
  safe_write_csv(details, file.path(OUTPUT_DIR, "pair_path_theory_validation_details.csv"))
  write_latex_summary(summary, file.path(OUTPUT_DIR, "pair_path_theory_validation_table.tex"))

  message("Summary written to: ", file.path(OUTPUT_DIR, "pair_path_theory_validation_summary.csv"))
  invisible(summary)
}

message("SIGA pair-path theory validation")
message("Profile: ", PROFILE)
message("Mode: ", RUN_MODE)
message("Scenario set: ", SCENARIO_SET)
message("Outer trials per scenario: ", N_OUTER)
message("Rerandomizations per trial: ", N_RERANDOMIZATIONS)
message("Allocation-only calibration paths: ", N_CALIBRATION)
message("Selected scenarios: ", paste(SCENARIOS$scenario_id, collapse = ","))
message("Output directory: ", OUTPUT_DIR)

if (RUN_MODE == "run") {
  for (i in seq_len(nrow(SCENARIOS))) run_scenario_shard(SCENARIOS[i, ])
} else if (RUN_MODE == "aggregate") {
  aggregate_results()
}
