#!/usr/bin/env Rscript

# R3 simulation: weak-null heterogeneity in the main-study designs, the
# variance-corrected reference test (CRT), and Monte Carlo standard errors.
#
# Methods compared on common outer trials:
#   RT          fixed-score reference randomisation test (plus-one Monte Carlo p-value)
#   CRT-model   RT with critical value rescaled by lambda = sqrt(V_R/V_S), model-based d(b)
#   CRT-data    same with the data-based d_hat(b) (stratum mean differences)
#   CRT-data-db same with the noise-debiased quadratic form
#   SIGA-S      Gaussian tail with the sampling variance V_S
#   SIGA-R      Gaussian tail with the conditional reference variance V_R (model-based d)
#
# Scenario sets (R3_SCENARIO_SET):
#   null         heterogeneous effects at the four null values (superiority 0,
#                non-inferiority margin, both equivalence limits), first-factor direction
#   power        heterogeneous effects at the manuscript alternatives
#   interaction  as "null" but with the highest-order-interaction direction (stress)
#   manuscript   the original homogeneous / common-log-odds models (all seven roles)
#   all          everything
#
# Base R only.  Environment variables are prefixed R3_.  See README.md.

options(stringsAsFactors = FALSE, warn = 1)

script_directory <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
BASE_DIR <- script_directory()
source(file.path(BASE_DIR, "siga_pair_path_engine.R"), local = FALSE)
source(file.path(BASE_DIR, "r3_engine_additions.R"), local = FALSE)

# -----------------------------------------------------------------------------
# Utilities (file handling identical in spirit to the theory-extension driver)
# -----------------------------------------------------------------------------
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
      tmp <- tempfile(pattern = paste0(basename(path), ".batch_"), tmpdir = parent, fileext = ".tmp")
      write.table(data, file = tmp, sep = ",", row.names = FALSE, col.names = !exists,
                  append = FALSE, quote = TRUE, qmethod = "double")
      if (exists) {
        if (!file.append(path, tmp)) stop("file.append returned FALSE")
        unlink(tmp)
      } else if (!file.rename(tmp, path)) {
        if (!file.copy(tmp, path, overwrite = FALSE)) stop("Could not move the first batch into place")
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
  stop("Failed to append checkpoint: ", path, if (!is.null(last_error)) paste0("; ", last_error) else "", call. = FALSE)
}
safe_write_csv <- function(data, path) {
  parent <- ensure_directory(dirname(path))
  tmp <- tempfile(pattern = basename(path), tmpdir = parent, fileext = ".tmp")
  write.csv(data, tmp, row.names = FALSE, na = "")
  if (!file.rename(tmp, path)) {
    if (!file.copy(tmp, path, overwrite = TRUE)) { unlink(tmp); stop("Could not write ", path, call. = FALSE) }
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

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
PROFILE <- tolower(trimws(Sys.getenv("R3_PROFILE", unset = "pilot")))
if (!PROFILE %in% c("smoke", "pilot", "manuscript")) stop("R3_PROFILE must be smoke, pilot, or manuscript.", call. = FALSE)
profile_defaults <- switch(
  PROFILE,
  smoke = c(outer = 20L, rerand = 99L, calibration = 500L, batch = 10L),
  pilot = c(outer = 2000L, rerand = 999L, calibration = 20000L, batch = 20L),
  manuscript = c(outer = 100000L, rerand = 4999L, calibration = 100000L, batch = 10L)
)
RUN_MODE <- tolower(trimws(Sys.getenv("R3_MODE", unset = "run")))
if (!RUN_MODE %in% c("calibrate", "run", "aggregate")) stop("R3_MODE must be calibrate, run, or aggregate.", call. = FALSE)
SCENARIO_SET <- tolower(trimws(Sys.getenv("R3_SCENARIO_SET", unset = "null")))
if (!SCENARIO_SET %in% c("null", "power", "interaction", "manuscript", "all")) {
  stop("R3_SCENARIO_SET must be null, power, interaction, manuscript, or all.", call. = FALSE)
}
N_OUTER <- env_integer("R3_N_OUTER", profile_defaults["outer"])
N_RERANDOMIZATIONS <- env_integer("R3_N_RERAND", profile_defaults["rerand"])
N_CALIBRATION <- env_integer("R3_N_CALIBRATION", profile_defaults["calibration"])
CALIBRATION_BATCH <- env_integer("R3_CALIBRATION_BATCH", 1000L)
OUTER_BATCH <- env_integer("R3_OUTER_BATCH", 0L)   # 0 = choose from the number of cores (below)
N_SHARDS <- env_integer("R3_N_SHARDS", 1L)
SHARD_ID <- env_integer("R3_SHARD_ID", 1L)
BASE_SEED <- env_integer("R3_SEED", 20261002L)
SCENARIO_FILTER <- env_integer_vector("R3_SCENARIO_IDS")
ALLOW_PARTIAL <- identical(Sys.getenv("R3_ALLOW_PARTIAL", unset = "0"), "1")
CONT_MAX_D <- env_numeric("R3_CONT_MAX_D", 0.50)        # max_s |d_s| for continuous heterogeneity
CONT_ETA_SD <- env_numeric("R3_CONT_ETA_SD", 0.25)      # individual effect SD
BIN_MAX_D <- env_numeric("R3_BIN_MAX_D", 0.15)          # max_s |d_s| for binary heterogeneity (risk-difference scale)
BINARY_CONTINUITY <- identical(Sys.getenv("R3_BINARY_CONTINUITY", unset = "0"), "1")
PBC_OVERRIDE <- env_numeric("R3_PBC", 0.80)             # biased-coin probability for all designs
CONT_SD <- env_numeric("R3_CONT_SD", 1.0)               # continuous outcome standard deviation
if (!(PBC_OVERRIDE > 0.5 && PBC_OVERRIDE <= 1)) stop("R3_PBC must lie in (0.5, 1].", call. = FALSE)

physical_cores <- parallel::detectCores(logical = FALSE)
if (is.na(physical_cores)) physical_cores <- parallel::detectCores(logical = TRUE)
if (is.na(physical_cores)) physical_cores <- 1L
DEFAULT_CORES <- if (.Platform$OS.type == "windows") 1L else max(1L, min(8L, physical_cores - 1L))
N_CORES <- env_integer("R3_CORES", DEFAULT_CORES)
if (.Platform$OS.type == "windows") N_CORES <- 1L
USE_PARALLEL <- .Platform$OS.type == "unix" && N_CORES > 1L
# The batch is the unit of parallel work and of checkpointing: mclapply runs at most one
# trial per core within a batch, so the batch must be a multiple of the core count.
if (OUTER_BATCH <= 0L) OUTER_BATCH <- max(as.integer(profile_defaults["batch"]), 10L * N_CORES)

PROJECT_DIR <- ensure_directory(Sys.getenv("R3_PROJECT_DIR", unset = file.path(dirname(BASE_DIR), "r3_output")))
OUTPUT_ROOT <- ensure_directory(Sys.getenv("R3_OUTPUT_DIR", unset = file.path(PROJECT_DIR, "simulation_output")))
RUN_TAG <- paste0("profile_", PROFILE, "_M", N_OUTER, "_B", N_RERANDOMIZATIONS, "_B0", N_CALIBRATION, "_S", N_SHARDS,
                  "_pbc", PBC_OVERRIDE, "_sd", CONT_SD, "_dC", CONT_MAX_D, "_dB", BIN_MAX_D)
OUTPUT_DIR <- ensure_directory(file.path(OUTPUT_ROOT, RUN_TAG))
CALIBRATION_DIR <- ensure_directory(file.path(PROJECT_DIR, "calibration_cache"))
SCENARIO_DIR <- ensure_directory(file.path(OUTPUT_DIR, "scenario_shards"))

stopifnot(N_OUTER >= 1L, N_RERANDOMIZATIONS >= 1L, N_CALIBRATION >= 2L, CALIBRATION_BATCH >= 1L,
          OUTER_BATCH >= 1L, N_SHARDS >= 1L, SHARD_ID >= 1L, SHARD_ID <= N_SHARDS, N_CORES >= 1L)
if (USE_PARALLEL && OUTER_BATCH < N_CORES) warning("R3_OUTER_BATCH is smaller than the number of cores; only ", OUTER_BATCH, " cores will be used.")

seed_value <- function(scenario_id, replicate_id, stream_id) {
  modulus <- 2147483646
  value <- as.double(BASE_SEED) + 1000003 * as.double(stream_id) +
    104729 * as.double(scenario_id) + 1009 * as.double(replicate_id)
  as.integer(value %% modulus + 1)
}
calibration_seed <- function(design_id) {
  modulus <- 2147483646
  as.integer((as.double(BASE_SEED) + 7000003 + 99991 * as.double(design_id)) %% modulus + 1)
}

# -----------------------------------------------------------------------------
# Designs (the four main-study designs)
# -----------------------------------------------------------------------------
FACTOR_PREVALENCE_MASTER <- c(0.50, 0.40, 0.30, 0.20, 0.10)
DESIGN_TABLE <- data.frame(
  design_id = 1:4,
  design_label = c("K2_n200", "K2_n1000", "K5_n400", "K5_n2000"),
  factor_count = c(2L, 2L, 5L, 5L),
  total_n = c(200L, 1000L, 400L, 2000L),
  pbc = PBC_OVERRIDE,
  binary_equivalence_limit = c(0.21, 0.10, 0.15, 0.10),
  cont_sup_power = c(0.430, 0.190, 0.300, 0.135),
  cont_ni_power = c(0.230, -0.010, 0.100, -0.065),
  cont_eq_power = c(0.000, 0.280, 0.180, 0.330),
  bin_sup_power = c(0.20, 0.10, 0.14, 0.07),
  bin_ni_power = c(0.10, 0.00, 0.05, -0.03),
  bin_eq_power = c(0.00, 0.00, 0.00, 0.04),
  stringsAsFactors = FALSE
)
CONT_NI_MARGIN <- 0.20
CONT_EQ_LIMIT <- 0.45
BIN_NI_MARGIN <- 0.10

prepare_design <- function(row) {
  K <- as.integer(row$factor_count)
  patterns <- all_binary_patterns(K)
  factor_prob <- FACTOR_PREVALENCE_MASTER[seq_len(K)]
  profile_prob <- profile_prob_independent(factor_prob, patterns)
  list(
    design_id = as.integer(row$design_id), design_label = row$design_label,
    factor_count = K, total_n = as.integer(row$total_n), pbc = as.numeric(row$pbc),
    weights = rep(1, K + 1L), patterns = patterns, factor_prob = factor_prob,
    profile_prob = profile_prob, row = row
  )
}
DESIGNS <- setNames(lapply(seq_len(nrow(DESIGN_TABLE)), function(i) prepare_design(DESIGN_TABLE[i, ])), DESIGN_TABLE$design_id)

# -----------------------------------------------------------------------------
# Scenario grid: 4 designs x 2 outcomes x 7 roles x 3 model variants = 168 ids
# -----------------------------------------------------------------------------
ROLES <- data.frame(
  objective = c("superiority", "superiority", "noninferiority", "noninferiority", "equivalence", "equivalence", "equivalence"),
  role = c("type1", "power", "type1", "power", "type1_lower", "type1_upper", "power"),
  stringsAsFactors = FALSE
)
VARIANTS <- c("reference", "hetero_first_factor", "hetero_interaction")

scenario_true_effect <- function(design_row, outcome_type, objective, role) {
  if (outcome_type == "continuous") {
    switch(paste(objective, role),
      "superiority type1" = 0, "superiority power" = design_row$cont_sup_power,
      "noninferiority type1" = -CONT_NI_MARGIN, "noninferiority power" = design_row$cont_ni_power,
      "equivalence type1_lower" = -CONT_EQ_LIMIT, "equivalence type1_upper" = CONT_EQ_LIMIT,
      "equivalence power" = design_row$cont_eq_power)
  } else {
    lim <- design_row$binary_equivalence_limit
    switch(paste(objective, role),
      "superiority type1" = 0, "superiority power" = design_row$bin_sup_power,
      "noninferiority type1" = -BIN_NI_MARGIN, "noninferiority power" = design_row$bin_ni_power,
      "equivalence type1_lower" = -lim, "equivalence type1_upper" = lim,
      "equivalence power" = design_row$bin_eq_power)
  }
}
scenario_boundaries <- function(design_row, outcome_type, objective) {
  if (objective == "superiority") return(c(zero = 0))
  if (objective == "noninferiority") return(c(lower = if (outcome_type == "continuous") -CONT_NI_MARGIN else -BIN_NI_MARGIN))
  lim <- if (outcome_type == "continuous") CONT_EQ_LIMIT else design_row$binary_equivalence_limit
  c(lower = -lim, upper = lim)
}
scenario_alpha <- function(objective) switch(objective, superiority = 0.05, noninferiority = 0.025, equivalence = 0.05)

make_scenario_grid <- function() {
  rows <- list(); id <- 0L
  for (d in seq_len(nrow(DESIGN_TABLE))) for (oc in c("continuous", "binary")) for (r in seq_len(nrow(ROLES))) for (v in VARIANTS) {
    id <- id + 1L
    drow <- DESIGN_TABLE[d, ]
    rows[[id]] <- data.frame(
      scenario_id = id, design_id = drow$design_id, design_label = drow$design_label,
      factor_count = drow$factor_count, total_n = drow$total_n, pbc = drow$pbc,
      outcome_type = oc, objective = ROLES$objective[r], role = ROLES$role[r], variant = v,
      true_effect = scenario_true_effect(drow, oc, ROLES$objective[r], ROLES$role[r]),
      alpha = scenario_alpha(ROLES$objective[r]),
      is_null = ROLES$role[r] != "power",
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  out$direction <- ifelse(out$variant == "hetero_first_factor", "first_factor",
                   ifelse(out$variant == "hetero_interaction", "interaction", "zero"))
  out$model_type <- ifelse(out$variant == "reference",
                           ifelse(out$outcome_type == "continuous", "homogeneous", "common_log_odds"),
                           ifelse(out$outcome_type == "continuous", "continuous_hetero", "binary_direct_hetero"))
  out$set <- ifelse(out$variant == "reference", "manuscript",
             ifelse(out$variant == "hetero_first_factor" & out$is_null, "null",
             ifelse(out$variant == "hetero_first_factor" & !out$is_null, "power",
             ifelse(out$variant == "hetero_interaction" & out$is_null, "interaction", "interaction_power"))))
  out$scenario_label <- paste0(out$design_label, "_", substr(out$outcome_type, 1L, 4L), "_",
                               out$objective, "_", out$role, "_", out$variant)
  out
}
SCENARIO_GRID <- make_scenario_grid()
SCENARIOS <- if (SCENARIO_SET == "all") SCENARIO_GRID else SCENARIO_GRID[SCENARIO_GRID$set == SCENARIO_SET, , drop = FALSE]
if (length(SCENARIO_FILTER)) SCENARIOS <- SCENARIO_GRID[SCENARIO_GRID$scenario_id %in% SCENARIO_FILTER, , drop = FALSE]
if (!nrow(SCENARIOS)) stop("No scenarios remain after filtering.", call. = FALSE)
rownames(SCENARIOS) <- NULL

# -----------------------------------------------------------------------------
# Three-copy allocation calibration (cached per design)
# -----------------------------------------------------------------------------
calibration_cache_file <- function(design) {
  file.path(CALIBRATION_DIR, paste0("pair_calibration_", design$design_label, "_p", design$pbc,
                                    "_B0", N_CALIBRATION, "_seed", calibration_seed(design$design_id), ".rds"))
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
  message("Creating three-copy allocation calibration for ", design$design_label, ": n=", design$total_n, ", B0=", N_CALIBRATION)
  calibration <- calibrate_pair_path_design_R(
    B0 = N_CALIBRATION, n = design$total_n, patterns = design$patterns, profile_prob = design$profile_prob,
    pbc = design$pbc, weights_ = design$weights, seed = calibration_seed(design$design_id),
    batch_size = CALIBRATION_BATCH, progress = TRUE
  )
  calibration$directions <- generalized_pair_directions(calibration$psi, calibration$profile_prob)
  tmp <- tempfile(pattern = "pair_calibration_", tmpdir = CALIBRATION_DIR, fileext = ".rds")
  saveRDS(calibration, tmp)
  if (!file.rename(tmp, cache)) {
    if (!file.copy(tmp, cache, overwrite = TRUE)) { unlink(tmp); stop("Could not save calibration: ", cache, call. = FALSE) }
    unlink(tmp)
  }
  calibration
}
USED_DESIGN_IDS <- sort(unique(SCENARIOS$design_id))
CALIBRATIONS <- setNames(lapply(USED_DESIGN_IDS, function(id) get_design_calibration(DESIGNS[[as.character(id)]])), USED_DESIGN_IDS)

# -----------------------------------------------------------------------------
# Outcome models
# -----------------------------------------------------------------------------
build_scenario_model <- function(scenario) {
  design <- DESIGNS[[as.character(scenario$design_id)]]
  calibration <- CALIBRATIONS[[as.character(scenario$design_id)]]
  v <- make_direction(scenario$direction, design$patterns, design$profile_prob, calibration)
  if (scenario$outcome_type == "continuous") {
    d0 <- if (scenario$variant == "reference") rep(0, nrow(design$patterns)) else scale_centred_direction(v, design$profile_prob, CONT_MAX_D)
    model <- make_continuous_model_r3(design$patterns, design$profile_prob, scenario$true_effect, d0,
                                      outcome_sd = CONT_SD, individual_effect_sd = CONT_ETA_SD)
  } else if (scenario$variant == "reference") {
    model <- make_binary_common_log_odds_model_r3(design$patterns, design$profile_prob, scenario$true_effect)
  } else {
    d0 <- scale_centred_direction(v, design$profile_prob, BIN_MAX_D)
    model <- make_binary_direct_model_r3(design$patterns, design$profile_prob, scenario$true_effect, d0)
  }
  # diagnostics at each tested null value
  bounds <- scenario_boundaries(design$row, scenario$outcome_type, scenario$objective)
  diag <- lapply(bounds, function(b) {
    d <- model_based_deviation(model, b)
    d_pi <- sum(design$profile_prob * d^2)
    d_psi <- drop(crossprod(d, calibration$psi %*% d))
    c(d_pi = d_pi, d_psi = d_psi, ratio = if (d_pi > 0) d_psi / d_pi else 1)
  })
  model$boundaries <- bounds
  model$boundary_diagnostics <- diag
  model
}
SCENARIO_MODELS <- setNames(lapply(seq_len(nrow(SCENARIOS)), function(i) build_scenario_model(SCENARIOS[i, ])), SCENARIOS$scenario_id)

scenario_definition_table <- function() {
  rows <- lapply(seq_len(nrow(SCENARIOS)), function(i) {
    sc <- SCENARIOS[i, ]; m <- SCENARIO_MODELS[[as.character(sc$scenario_id)]]
    bd <- m$boundary_diagnostics
    data.frame(
      sc[, c("scenario_id", "scenario_label", "set", "design_label", "factor_count", "total_n", "pbc",
             "outcome_type", "objective", "role", "variant", "model_type", "direction", "true_effect", "alpha")],
      boundaries = paste(sprintf("%s=%.3f", names(m$boundaries), m$boundaries), collapse = ";"),
      achieved_effect = m$achieved_effect,
      max_abs_d0 = max(abs(m$d0)),
      heterogeneity_shrink = m$heterogeneity_shrink %||% 1,
      individual_effect_sd = m$individual_effect_sd %||% NA_real_,
      pair_ratio_at_first_boundary = unname(bd[[1L]]["ratio"]),
      pair_ratio_at_last_boundary = unname(bd[[length(bd)]]["ratio"]),
      control_risk = m$achieved_control_risk %||% NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
design_diagnostics_table <- function() {
  do.call(rbind, lapply(USED_DESIGN_IDS, function(id) {
    d <- DESIGNS[[as.character(id)]]; cal <- CALIBRATIONS[[as.character(id)]]
    data.frame(design_id = id, design_label = d$design_label, factor_count = d$factor_count, total_n = d$total_n,
               pbc = d$pbc, calibration_paths = cal$B0, calibration_seconds = cal$elapsed_seconds,
               min_generalized_ratio = cal$directions$min_ratio, max_generalized_ratio = cal$directions$max_ratio,
               stringsAsFactors = FALSE)
  }))
}
safe_write_csv(scenario_definition_table(), file.path(OUTPUT_DIR, "r3_scenario_definitions.csv"))
safe_write_csv(design_diagnostics_table(), file.path(OUTPUT_DIR, "r3_design_diagnostics.csv"))
if (RUN_MODE == "calibrate") {
  message("Calibration and scenario definitions written to: ", OUTPUT_DIR)
  quit(save = "no", status = 0L)
}

# -----------------------------------------------------------------------------
# One outer trial
# -----------------------------------------------------------------------------
METHODS <- c("rt", "crt_model", "crt_data", "crt_data_db", "siga_s", "siga_r")
ANALYSES <- c("unadjusted", "adjusted")
LAMBDA_COLS <- c("rt", "crt_model", "crt_data", "crt_data_db")

simulate_trial <- function(scenario, replicate_id, design, model) {
  set.seed(seed_value(scenario$scenario_id, replicate_id, 1L))
  profile <- generate_profile_sequence(design$total_n, design$patterns, design$profile_prob)
  z <- ps_assign_R(profile$X, pbc = design$pbc, weights_ = design$weights,
                   seed = seed_value(scenario$scenario_id, replicate_id, 2L))
  A <- as.integer((z + 1L) / 2L)
  set.seed(seed_value(scenario$scenario_id, replicate_id, 3L))
  y <- generate_outcome_from_model(profile$id, A, z, model)
  list(profile_id = profile$id, X = profile$X, z = as.integer(z), A = A, y = y)
}

gaussian_p <- function(statistic, variance, alternative) gaussian_pvalue(statistic, variance, alternative)$p

run_replicate <- function(scenario, replicate_id, calibration, design, model) {
  t0 <- proc.time()[3L]
  trial <- simulate_trial(scenario, replicate_id, design, model)
  J <- calibration$J
  bounds <- model$boundaries
  nb <- length(bounds)
  alternatives <- if (scenario$objective == "superiority") "two.sided" else if (scenario$objective == "noninferiority") "greater" else c("greater", "less")

  # per boundary and analysis: residual vectors, variances, lambdas, Gaussian p-values
  score_mat <- matrix(0, nrow = design$total_n, ncol = 2L * nb)
  lambda <- matrix(1, nrow = 2L * nb, ncol = length(LAMBDA_COLS), dimnames = list(NULL, LAMBDA_COLS))
  col_names <- character(2L * nb)
  per <- list()
  for (bi in seq_len(nb)) {
    b <- bounds[bi]; bname <- names(bounds)[bi]
    scores <- make_boundary_scores(trial$y, trial$A, trial$X, b)
    d_model <- model_based_deviation(model, b)
    dh <- estimate_stratum_effect_deviations(trial$y, trial$A, trial$profile_id, J, b)
    for (ai in seq_along(ANALYSES)) {
      a <- ANALYSES[ai]
      col <- 2L * (bi - 1L) + ai
      col_names[col] <- paste0(bname, "_", a)
      sc <- scores[[a]]
      score_mat[, col] <- sc
      statistic <- 0.5 * sum(trial$z * sc)
      samp <- siga_sampling_variance(sc, trial$X, calibration)
      vr_model <- corrected_variance(samp, d_model, calibration = calibration)
      vr_data <- corrected_variance(samp, dh$d, calibration = calibration)
      vr_data_db <- corrected_variance(samp, dh$d, var_d = dh$var_d, calibration = calibration, debias = TRUE)
      lambda[col, ] <- c(1, vr_model$lambda, vr_data$lambda, vr_data_db$lambda)
      per[[col_names[col]]] <- list(
        boundary = b, analysis = a, alternative = alternatives[bi], statistic = statistic,
        v_s = samp$variance, v_r_model = vr_model$variance, v_r_model_raw = vr_model$variance_raw,
        v_r_data = vr_data$variance, v_r_data_db = vr_data_db$variance,
        lambda_model = vr_model$lambda, lambda_data = vr_data$lambda, lambda_data_db = vr_data_db$lambda,
        floored_model = vr_model$floored, floored_data = vr_data$floored, floored_data_db = vr_data_db$floored,
        kappa_raw = samp$kappa_raw,
        p_siga_s = gaussian_p(statistic, samp$variance, alternatives[bi]),
        p_siga_r = gaussian_p(statistic, vr_model$variance, alternatives[bi]),
        dhat_max_abs_error = max(abs(dh$d - d_model)),
        dhat_rms_error = sqrt(sum(design$profile_prob * (dh$d - d_model)^2)),
        strata_both_arms = sum(dh$both_arms)
      )
      # optional lattice continuity correction (unadjusted binary superiority at b = 0 only)
      if (BINARY_CONTINUITY && scenario$outcome_type == "binary" && scenario$objective == "superiority" && a == "unadjusted") {
        yp <- as.integer(sum(trial$y))
        per[[col_names[col]]]$p_siga_s <- lattice_mixture_pvalue(statistic, samp$variance, design$total_n, yp, calibration$treated_count_prob, "two.sided")
        per[[col_names[col]]]$p_siga_r <- lattice_mixture_pvalue(statistic, vr_model$variance, design$total_n, yp, calibration$treated_count_prob, "two.sided")
      }
    }
  }
  colnames(score_mat) <- col_names
  rownames(lambda) <- col_names
  t_siga <- proc.time()[3L]

  rt <- rt_corrected_pvalues(trial$X, score_mat, trial$z, lambda, B = N_RERANDOMIZATIONS, pbc = design$pbc,
                               weights_ = design$weights, seed = seed_value(scenario$scenario_id, replicate_id, 4L))
  t_rt <- proc.time()[3L]

  # combine boundaries into the objective-level p-value per analysis and method
  pick <- function(matrix_list, col, lam) {
    alt <- per[[col]]$alternative
    switch(alt, two.sided = rt$two_sided[col, lam], greater = rt$greater[col, lam], less = rt$less[col, lam])
  }
  out <- list(
    scenario_id = scenario$scenario_id, scenario_label = scenario$scenario_label, replicate = replicate_id,
    set = scenario$set, design_id = scenario$design_id, design_label = design$design_label,
    factor_count = design$factor_count, total_n = design$total_n, pbc = design$pbc,
    outcome_type = scenario$outcome_type, objective = scenario$objective, role = scenario$role,
    variant = scenario$variant, model_type = scenario$model_type, direction = scenario$direction,
    true_effect = scenario$true_effect, alpha = scenario$alpha, is_null = scenario$is_null,
    configured_outer_trials = N_OUTER, rerandomizations_per_trial = N_RERANDOMIZATIONS,
    allocation_calibration_paths = N_CALIBRATION, configured_shards = N_SHARDS,
    binary_continuity = BINARY_CONTINUITY
  )
  for (a in ANALYSES) {
    cols <- paste0(names(bounds), "_", a)
    p_rt <- sapply(LAMBDA_COLS, function(lam) {
      ps <- sapply(cols, function(col) pick(rt, col, lam))
      if (scenario$objective == "equivalence") max(ps) else ps[[1L]]
    })
    p_s <- sapply(cols, function(col) per[[col]]$p_siga_s); p_s <- if (scenario$objective == "equivalence") max(p_s) else p_s[[1L]]
    p_r <- sapply(cols, function(col) per[[col]]$p_siga_r); p_r <- if (scenario$objective == "equivalence") max(p_r) else p_r[[1L]]
    p_all <- c(rt = unname(p_rt["rt"]), crt_model = unname(p_rt["crt_model"]), crt_data = unname(p_rt["crt_data"]),
               crt_data_db = unname(p_rt["crt_data_db"]), siga_s = p_s, siga_r = p_r)
    for (m in METHODS) {
      out[[paste0("p_", m, "_", a)]] <- p_all[[m]]
      out[[paste0("reject_", m, "_", a)]] <- p_all[[m]] <= scenario$alpha
    }
    # per-boundary diagnostics
    for (col in cols) {
      pr <- per[[col]]; bname <- sub(paste0("_", a, "$"), "", col)
      pre <- paste0(bname, "_", a, "_")
      out[[paste0(pre, "statistic")]] <- pr$statistic
      out[[paste0(pre, "v_s")]] <- pr$v_s
      out[[paste0(pre, "v_r_model")]] <- pr$v_r_model
      out[[paste0(pre, "v_r_data")]] <- pr$v_r_data
      out[[paste0(pre, "v_r_data_db")]] <- pr$v_r_data_db
      out[[paste0(pre, "v_rt_empirical")]] <- unname(rt$randomization_variance[col])
      out[[paste0(pre, "lambda_model")]] <- pr$lambda_model
      out[[paste0(pre, "lambda_data")]] <- pr$lambda_data
      out[[paste0(pre, "lambda_data_db")]] <- pr$lambda_data_db
      out[[paste0(pre, "floored_any")]] <- pr$floored_model || pr$floored_data || pr$floored_data_db
      out[[paste0(pre, "kappa_raw")]] <- pr$kappa_raw
      out[[paste0(pre, "dhat_max_abs_error")]] <- pr$dhat_max_abs_error
      out[[paste0(pre, "dhat_rms_error")]] <- pr$dhat_rms_error
      out[[paste0(pre, "strata_both_arms")]] <- pr$strata_both_arms
    }
  }
  out$observed_treated <- sum(trial$A)
  out$siga_seconds <- t_siga - t0
  out$rt_seconds <- t_rt - t_siga
  as.data.frame(out, stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------
# Sharded execution
# -----------------------------------------------------------------------------
scenario_paths <- function(scenario) {
  directory <- ensure_directory(file.path(SCENARIO_DIR, sprintf("scenario_%03d_%s", scenario$scenario_id, scenario$scenario_label)))
  list(directory = directory, shard = file.path(directory, sprintf("shard_%03d_of_%03d.csv", SHARD_ID, N_SHARDS)))
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
  message("Scenario ", scenario$scenario_id, " [", scenario$scenario_label, "]: shard ", SHARD_ID, "/", N_SHARDS,
          "; assigned=", length(assigned), "; completed=", length(completed), "; pending=", length(pending))
  if (!length(pending)) return(invisible(NULL))
  start <- proc.time()[3L]
  for (from in seq.int(1L, length(pending), by = OUTER_BATCH)) {
    ids <- pending[from:min(length(pending), from + OUTER_BATCH - 1L)]
    worker <- function(id) run_replicate(scenario, id, calibration, design, model)
    rows <- if (USE_PARALLEL && length(ids) > 1L) parallel::mclapply(ids, worker, mc.cores = min(N_CORES, length(ids))) else lapply(ids, worker)
    failed <- vapply(rows, function(r) inherits(r, "try-error") || !is.data.frame(r), logical(1L))
    if (any(failed)) stop("A worker failed in scenario ", scenario$scenario_id, ": ", paste(unlist(rows[failed]), collapse = "; "), call. = FALSE)
    append_csv(do.call(rbind, rows), paths$shard)
    message("  completed ", length(completed) + min(length(pending), from + length(ids) - 1L), "/", length(assigned),
            " assigned trials; elapsed ", format_elapsed(proc.time()[3L] - start))
  }
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# Aggregation with Monte Carlo standard errors
# -----------------------------------------------------------------------------
read_scenario_details <- function(scenario) {
  directory <- file.path(SCENARIO_DIR, sprintf("scenario_%03d_%s", scenario$scenario_id, scenario$scenario_label))
  if (!dir.exists(directory)) stop("Scenario directory is missing: ", directory, call. = FALSE)
  files <- list.files(directory, pattern = sprintf("^shard_[0-9]+_of_%03d\\.csv$", N_SHARDS), full.names = TRUE)
  if (!length(files)) stop("No shard files found for ", scenario$scenario_label, call. = FALSE)
  dat <- do.call(rbind, lapply(files, function(p) read.csv(p, stringsAsFactors = FALSE)))
  dat$replicate <- as.integer(dat$replicate)
  dat <- dat[order(dat$replicate), , drop = FALSE]
  dat <- dat[!duplicated(dat$replicate), , drop = FALSE]
  if (!ALLOW_PARTIAL && nrow(dat) != N_OUTER) stop(scenario$scenario_label, " contains ", nrow(dat), " unique replicates; expected ", N_OUTER, call. = FALSE)
  dat
}

summarize_analysis <- function(dat, scenario, a) {
  M <- nrow(dat)
  rej <- lapply(METHODS, function(m) as.logical(dat[[paste0("reject_", m, "_", a)]]))
  names(rej) <- METHODS
  out <- data.frame(analysis = a, n_outer = M, stringsAsFactors = FALSE)
  for (m in METHODS) {
    p <- mean(rej[[m]], na.rm = TRUE)
    out[[m]] <- p
    out[[paste0(m, "_se")]] <- binomial_se(p, M)
  }
  for (m in setdiff(METHODS, "rt")) {
    out[[paste0("diff_", m, "_minus_rt")]] <- mean(rej[[m]] - rej[["rt"]], na.rm = TRUE)
    out[[paste0("diff_", m, "_minus_rt_se")]] <- paired_difference_se(rej[[m]], rej[["rt"]])
  }
  out[["diff_siga_s_minus_crt_data"]] <- mean(rej[["siga_s"]] - rej[["crt_data"]], na.rm = TRUE)
  out[["diff_siga_s_minus_crt_data_se"]] <- paired_difference_se(rej[["siga_s"]], rej[["crt_data"]])
  if (isTRUE(scenario$is_null)) {
    for (m in METHODS) out[[paste0(m, "_minus_nominal")]] <- out[[m]] - scenario$alpha
  }
  # variance diagnostics at the first boundary (the tested null value for type-1 rows)
  bounds <- SCENARIO_MODELS[[as.character(scenario$scenario_id)]]$boundaries
  bname <- if (isTRUE(scenario$is_null)) names(bounds)[which.min(abs(bounds - scenario$true_effect))] else names(bounds)[1L]
  out$diagnostic_boundary <- bname
  pre <- paste0(bname, "_", a, "_")
  g <- function(x) as.numeric(dat[[paste0(pre, x)]])
  ratio <- function(x, y) { r <- x / y; mean(r[is.finite(r)]) }
  out$mean_lambda_model <- mean(g("lambda_model"), na.rm = TRUE)
  out$mean_lambda_data <- mean(g("lambda_data"), na.rm = TRUE)
  out$mean_lambda_data_db <- mean(g("lambda_data_db"), na.rm = TRUE)
  out$mean_ratio_vr_model_over_vs <- ratio(g("v_r_model"), g("v_s"))
  out$mean_ratio_vrt_empirical_over_vs <- ratio(g("v_rt_empirical"), g("v_s"))
  out$mean_ratio_vrt_empirical_over_vr_model <- ratio(g("v_rt_empirical"), g("v_r_model"))
  out$mean_dhat_max_abs_error <- mean(g("dhat_max_abs_error"), na.rm = TRUE)
  out$mean_dhat_rms_error <- mean(g("dhat_rms_error"), na.rm = TRUE)
  out$floor_rate <- mean(as.logical(dat[[paste0(pre, "floored_any")]]), na.rm = TRUE)
  out$mean_strata_both_arms <- mean(g("strata_both_arms"), na.rm = TRUE)
  out$siga_minutes <- sum(dat$siga_seconds, na.rm = TRUE) / 60
  out$rt_minutes <- sum(dat$rt_seconds, na.rm = TRUE) / 60
  out
}
summarize_scenario <- function(dat, scenario) {
  model <- SCENARIO_MODELS[[as.character(scenario$scenario_id)]]
  fixed <- scenario[, c("scenario_id", "scenario_label", "set", "design_label", "factor_count", "total_n", "pbc",
                        "outcome_type", "objective", "role", "variant", "model_type", "direction", "true_effect", "alpha", "is_null")]
  fixed$max_abs_d0 <- max(abs(model$d0))
  bsel <- if (isTRUE(scenario$is_null)) which.min(abs(model$boundaries - scenario$true_effect)) else 1L
  fixed$pair_ratio_model <- unname(model$boundary_diagnostics[[bsel]]["ratio"])
  by <- rbind(summarize_analysis(dat, scenario, "unadjusted"), summarize_analysis(dat, scenario, "adjusted"))
  cbind(fixed[rep(1L, nrow(by)), , drop = FALSE], by)
}

fmt_pct <- function(p, se, digits = 2) ifelse(is.finite(p), sprintf(paste0("%.", digits, "f (%.", digits, "f)"), 100 * p, 100 * se), "--")
latex_escape <- function(x) gsub("_", "\\\\_", x)

write_rejection_table <- function(summary, path, caption, label) {
  lines <- c("\\begin{table}[!htbp]", "\\centering", paste0("\\caption{", caption, "}"), paste0("\\label{", label, "}"),
             "\\small", "\\begin{tabular}{@{}lllllllll@{}}", "\\toprule",
             "Design & Outcome & Objective (role) & Analysis & RT & CRT (model $d$) & CRT (data $\\widehat d$) & SIGA-S & SIGA-R\\\\", "\\midrule")
  for (i in seq_len(nrow(summary))) {
    r <- summary[i, ]
    lines <- c(lines, paste0(
      latex_escape(r$design_label), " & ", r$outcome_type, " & ", r$objective, " (", latex_escape(r$role), ") & ", r$analysis, " & ",
      fmt_pct(r$rt, r$rt_se), " & ", fmt_pct(r$crt_model, r$crt_model_se), " & ", fmt_pct(r$crt_data, r$crt_data_se), " & ",
      fmt_pct(r$siga_s, r$siga_s_se), " & ", fmt_pct(r$siga_r, r$siga_r_se), "\\\\"))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}",
             "\\begin{flushleft}\\footnotesize Rejection probabilities in percent with Monte Carlo standard errors in parentheses. CRT denotes the reference test with the critical value rescaled by $\\widehat\\lambda_n$; the debiased data-based variant is reported in the CSV output.\\end{flushleft}",
             "\\end{table}")
  writeLines(lines, path)
}

write_table1_summary <- function(summary, path) {
  comparisons <- list(
    c("diff_siga_s_minus_rt", "SIGA-S versus RT"),
    c("diff_siga_r_minus_rt", "SIGA-R versus RT"),
    c("diff_crt_data_minus_rt", "CRT (data) versus RT"),
    c("diff_siga_s_minus_crt_data", "SIGA-S versus CRT (data)")
  )
  block <- function(rows, title) {
    out <- character()
    for (cmp in comparisons) {
      v <- rows[[cmp[1L]]]; se <- rows[[paste0(cmp[1L], "_se")]]
      if (!length(v) || all(!is.finite(v))) next
      k <- which.max(abs(v))
      out <- c(out, paste0(cmp[2L], " & ", title, " & ", sprintf("%.2f", 100 * abs(v[k])), " & ", sprintf("%.2f", 100 * se[k]),
                           " & ", latex_escape(rows$scenario_label[k]), " (", rows$analysis[k], ")\\\\"))
    }
    out
  }
  null_rows <- summary[summary$is_null, , drop = FALSE]
  power_rows <- summary[!summary$is_null, , drop = FALSE]
  nominal <- character()
  for (m in c("rt", "crt_model", "crt_data", "crt_data_db", "siga_s", "siga_r")) {
    v <- null_rows[[paste0(m, "_minus_nominal")]]; se <- null_rows[[paste0(m, "_se")]]
    if (!length(v)) next
    k <- which.max(abs(v))
    nominal <- c(nominal, paste0(toupper(gsub("_", "-", m)), " versus nominal level & null values & ", sprintf("%+.2f", 100 * v[k]), " & ",
                                 sprintf("%.2f", 100 * se[k]), " & ", latex_escape(null_rows$scenario_label[k]), " (", null_rows$analysis[k], ")\\\\"))
  }
  lines <- c("\\begin{table}[!htbp]", "\\centering",
             "\\caption{Largest absolute differences in rejection probability (percentage points) across scenarios and both residual vectors, with the Monte Carlo standard error of the difference at the maximising scenario.}",
             "\\label{tab:r3-summary}", "\\small", "\\begin{tabular}{@{}lllll@{}}", "\\toprule",
             "Comparison & Rows & Max $|$difference$|$ & SE & Scenario\\\\", "\\midrule",
             block(null_rows, "null values"), block(power_rows, "power"), "\\midrule", nominal,
             "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, path)
}

rbind_fill <- function(frames) {
  cols <- unique(unlist(lapply(frames, names)))
  do.call(rbind, lapply(frames, function(d) { for (cc in setdiff(cols, names(d))) d[[cc]] <- NA; d[, cols, drop = FALSE] }))
}

aggregate_results <- function() {
  summaries <- list(); details <- list()
  for (i in seq_len(nrow(SCENARIOS))) {
    sc <- SCENARIOS[i, ]
    message("Aggregating ", sc$scenario_label)
    dat <- read_scenario_details(sc)
    details[[i]] <- dat
    summaries[[i]] <- summarize_scenario(dat, sc)
  }
  summary <- rbind_fill(summaries)
  summary <- summary[order(summary$scenario_id, summary$analysis), , drop = FALSE]
  tag <- paste0("r3_", SCENARIO_SET)
  safe_write_csv(summary, file.path(OUTPUT_DIR, paste0(tag, "_summary.csv")))
  if (!identical(Sys.getenv("R3_SKIP_DETAILS", unset = "0"), "1")) {
    safe_write_csv(rbind_fill(details), file.path(OUTPUT_DIR, paste0(tag, "_details.csv")))
  }
  nulls <- summary[summary$is_null, , drop = FALSE]
  powers <- summary[!summary$is_null, , drop = FALSE]
  if (nrow(nulls)) write_rejection_table(nulls, file.path(OUTPUT_DIR, paste0(tag, "_null_table.tex")),
                                         "Rejection probabilities at null values under treatment-effect heterogeneity.", "tab:r3-null")
  if (nrow(powers)) write_rejection_table(powers, file.path(OUTPUT_DIR, paste0(tag, "_power_table.tex")),
                                          "Power at the fixed alternatives under treatment-effect heterogeneity.", "tab:r3-power")
  write_table1_summary(summary, file.path(OUTPUT_DIR, paste0(tag, "_table1_summary.tex")))
  message("Summary written to: ", file.path(OUTPUT_DIR, paste0(tag, "_summary.csv")))
  invisible(summary)
}

KERNEL_SETTING <- tolower(Sys.getenv("R3_KERNEL", unset = "auto"))
KERNEL_AVAILABLE <- if (KERNEL_SETTING %in% c("auto", "c")) load_rt_kernel(BASE_DIR, quiet = TRUE) else FALSE
if (KERNEL_SETTING == "c" && !KERNEL_AVAILABLE) stop("R3_KERNEL=C requested but rt_kernel.c could not be compiled; install a C compiler or set R3_KERNEL=R.", call. = FALSE)
message("R3 weak-null / corrected-reference-test simulation")
message("Regeneration kernel: ", if (KERNEL_AVAILABLE) "compiled C (rt_kernel.c)" else "pure R")
message("Profile: ", PROFILE, "; mode: ", RUN_MODE, "; scenario set: ", SCENARIO_SET)
message("Outer trials: ", N_OUTER, "; rerandomisations: ", N_RERANDOMIZATIONS, "; calibration paths: ", N_CALIBRATION)
message("Heterogeneity: continuous max|d| = ", CONT_MAX_D, " (eta sd ", CONT_ETA_SD, "); binary max|d| = ", BIN_MAX_D)
message("Biased coin p_bc: ", PBC_OVERRIDE, "; continuous outcome SD: ", CONT_SD)
message("Binary continuity correction: ", BINARY_CONTINUITY, "; cores: ", N_CORES, "; outer batch: ", OUTER_BATCH, "; shard ", SHARD_ID, "/", N_SHARDS)
message("Scenarios: ", paste(SCENARIOS$scenario_id, collapse = ","))
message("Output directory: ", OUTPUT_DIR)

if (RUN_MODE == "run") {
  for (i in seq_len(nrow(SCENARIOS))) run_scenario_shard(SCENARIOS[i, ])
} else if (RUN_MODE == "aggregate") {
  aggregate_results()
}
