#!/usr/bin/env Rscript

# =============================================================================
# Comprehensive SIGA-R validation for the manuscript benchmark
# =============================================================================
#
# Scientific purpose
# ------------------
# Re-run the complete continuous- and binary-outcome benchmark on common outer
# trials and compare:
#   1. SIGA-S, targeting the repeated-sampling distribution;
#   2. SIGA-R, targeting the conditional fixed-score randomization distribution;
#   3. the reference fixed-score randomization test (RT).
#
# The scenario grid exactly matches the manuscript benchmark:
#   * continuous and binary outcomes;
#   * 2 factors with 100 or 500 participants/group;
#   * 5 factors with 200 or 1000 participants/group;
#   * superiority, non-inferiority, and equivalence;
#   * type-I-error boundary scenarios and power scenarios;
#   * unadjusted and treatment-blind adjusted scores.
#
# Full manuscript settings (defaults)
# -----------------------------------
#   outer trials                    100,000 per scenario
#   regenerated RT paths             4,999 per outer trial
#   three-copy allocation calibration 100,000 per design/sample size
#   biased-coin probability             0.80
#
# Required companion file
# -----------------------
#   siga_pair_path_engine.R
# Place it in the same directory as this script, or set PWRT_PAIR_ENGINE.
#
# Modes
# -----
#   PWRT_MODE=calibrate  create the four reusable three-copy calibrations
#   PWRT_MODE=run        run/resume one outer-simulation shard
#
# Example smoke test
# ------------------
#   PWRT_N_OUTER=20 PWRT_N_RERAND=99 PWRT_N_CALIBRATION=500 \
#   PWRT_N_SHARDS=1 PWRT_SHARD_ID=1 PWRT_MODE=run \
#   Rscript simulation_siga_r_full_benchmark.R
#
# Full sharded run, shard 3 of 24
# --------------------------------
#   PWRT_MODE=run PWRT_N_SHARDS=24 PWRT_SHARD_ID=3 \
#   Rscript simulation_siga_r_full_benchmark.R
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

# Prevent implicit BLAS/OpenMP oversubscription when many R processes are used.
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

BASE_DIR <- script_directory()
PAIR_ENGINE <- path.expand(Sys.getenv(
  "PWRT_PAIR_ENGINE",
  unset = file.path(BASE_DIR, "siga_pair_path_engine.R")
))
if (!file.exists(PAIR_ENGINE)) {
  stop(
    "Required file not found: ", PAIR_ENGINE,
    "\nPlace siga_pair_path_engine.R beside this script or set PWRT_PAIR_ENGINE.",
    call. = FALSE
  )
}
source(PAIR_ENGINE, local = FALSE)

# -----------------------------------------------------------------------------
# Environment and file helpers
# -----------------------------------------------------------------------------

env_integer <- function(name, default) {
  text <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(text)) return(as.integer(default))
  value <- suppressWarnings(as.integer(text))
  if (is.na(value)) stop(name, " must be an integer.", call. = FALSE)
  value
}

env_numeric <- function(name, default) {
  text <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(text)) return(as.numeric(default))
  value <- suppressWarnings(as.numeric(text))
  if (!is.finite(value)) stop(name, " must be finite numeric.", call. = FALSE)
  value
}

env_integer_vector <- function(name) {
  text <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(text)) return(integer())
  value <- suppressWarnings(as.integer(trimws(strsplit(text, ",", fixed = TRUE)[[1L]])))
  if (anyNA(value)) stop(name, " must contain comma-separated integers.", call. = FALSE)
  unique(value)
}

ensure_directory <- function(path, attempts = 5L, wait_seconds = 1) {
  path <- path.expand(path)
  for (attempt in seq_len(attempts)) {
    if (dir.exists(path)) return(normalizePath(path, winslash = "/", mustWork = TRUE))
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (dir.exists(path)) return(normalizePath(path, winslash = "/", mustWork = TRUE))
    if (attempt < attempts) Sys.sleep(wait_seconds)
  }
  stop("Could not create directory: ", path, call. = FALSE)
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
          stop("Could not move first batch into place")
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
    "Could not append checkpoint: ", path,
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
      stop("Could not write: ", path, call. = FALSE)
    }
    unlink(tmp)
  }
  invisible(path)
}

# -----------------------------------------------------------------------------
# Fixed manuscript settings, with smoke-test overrides
# -----------------------------------------------------------------------------

RUN_MODE <- tolower(trimws(Sys.getenv("PWRT_MODE", unset = "run")))
if (!RUN_MODE %in% c("calibrate", "run")) {
  stop("PWRT_MODE must be 'calibrate' or 'run'.", call. = FALSE)
}

N_OUTER <- env_integer("PWRT_N_OUTER", 100000L)
N_RERANDOMIZATIONS <- env_integer("PWRT_N_RERAND", 4999L)
N_CALIBRATION <- env_integer("PWRT_N_CALIBRATION", 100000L)
CALIBRATION_BATCH <- env_integer("PWRT_CALIBRATION_BATCH", 1000L)
OUTER_BATCH <- env_integer("PWRT_OUTER_BATCH", 10L)
N_SHARDS <- env_integer("PWRT_N_SHARDS", 1L)
SHARD_ID <- env_integer("PWRT_SHARD_ID", 1L)
BASE_SEED <- env_integer("PWRT_SEED", 20260801L)
P_BIASED_COIN <- env_numeric("PWRT_P_BIASED_COIN", 0.80)
EPSILON_EXPONENT <- env_numeric("PWRT_EPSILON_EXPONENT", 1.0)
SCENARIO_FILTER <- env_integer_vector("PWRT_SCENARIO_IDS")

stopifnot(
  N_OUTER >= 1L,
  N_RERANDOMIZATIONS >= 1L,
  N_CALIBRATION >= 2L,
  CALIBRATION_BATCH >= 1L,
  OUTER_BATCH >= 1L,
  N_SHARDS >= 1L,
  SHARD_ID >= 1L,
  SHARD_ID <= N_SHARDS,
  P_BIASED_COIN > 0.5,
  P_BIASED_COIN < 1,
  EPSILON_EXPONENT > 0
)

if (N_OUTER != 100000L || N_RERANDOMIZATIONS != 4999L || N_CALIBRATION != 100000L) {
  warning(
    "Non-manuscript simulation settings are active: N_OUTER=", N_OUTER,
    ", N_RERAND=", N_RERANDOMIZATIONS,
    ", N_CALIBRATION=", N_CALIBRATION, "."
  )
}

OUTPUT_ROOT <- ensure_directory(path.expand(Sys.getenv(
  "PWRT_OUTPUT_DIR",
  unset = file.path(BASE_DIR, "siga_r_full_benchmark_output")
)))
RUN_TAG <- paste0(
  "M", N_OUTER,
  "_B", N_RERANDOMIZATIONS,
  "_Bpsi", N_CALIBRATION,
  "_pbc", gsub("[.]", "p", formatC(P_BIASED_COIN, format = "f", digits = 3)),
  "_eps", gsub("[.]", "p", formatC(EPSILON_EXPONENT, format = "f", digits = 2)),
  "_seed", BASE_SEED
)
RUN_DIR <- ensure_directory(file.path(OUTPUT_ROOT, RUN_TAG))
CALIBRATION_DIR <- ensure_directory(file.path(OUTPUT_ROOT, "pair_calibration_cache"))
SCENARIO_DIR <- ensure_directory(file.path(RUN_DIR, "scenario_shards"))
CONFIG_DIR <- ensure_directory(file.path(RUN_DIR, "configuration"))

# -----------------------------------------------------------------------------
# Reproducible random-number streams
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# Manuscript scenario grid
# -----------------------------------------------------------------------------

FACTOR_PREVALENCE_MASTER <- c(0.50, 0.40, 0.30, 0.20, 0.10)
CONTINUOUS_BETA_MASTER <- c(0.50, 0.40, 0.30, 0.20, 0.10)
CONTINUOUS_INTERACTION <- 0.25
CONTINUOUS_OUTCOME_SD <- 1.00
BINARY_BETA_MASTER <- c(0.35, -0.25, 0.20, -0.15, 0.10)
BINARY_INTERACTION <- 0.20
TARGET_CONTROL_RISK <- 0.60

make_design_table <- function() {
  data.frame(
    design_id = 1:4,
    factor_count = c(2L, 2L, 5L, 5L),
    target_per_group = c(100L, 500L, 200L, 1000L),
    total_n = c(200L, 1000L, 400L, 2000L),
    stringsAsFactors = FALSE
  )
}

make_scenarios <- function() {
  designs <- make_design_table()

  continuous_effects <- data.frame(
    design_id = 1:4,
    superiority = c(0.430, 0.190, 0.300, 0.135),
    noninferiority = c(0.230, -0.010, 0.100, -0.065),
    equivalence = c(0.000, 0.280, 0.180, 0.330),
    stringsAsFactors = FALSE
  )

  binary_effects <- data.frame(
    design_id = 1:4,
    superiority = c(0.200, 0.100, 0.140, 0.070),
    noninferiority = c(0.100, 0.000, 0.050, -0.030),
    equivalence = c(0.000, 0.000, 0.000, 0.040),
    equivalence_lower = c(-0.21, -0.10, -0.15, -0.10),
    equivalence_upper = c(0.21, 0.10, 0.15, 0.10),
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
    alpha = c(0.05, 0.05, 0.025, 0.025, 0.05, 0.05, 0.05),
    stringsAsFactors = FALSE
  )

  rows <- list()
  pos <- 0L
  for (outcome in c("continuous", "binary")) {
    for (d in seq_len(nrow(designs))) {
      design <- designs[d, ]
      for (r in seq_len(nrow(roles))) {
        role <- roles[r, ]
        if (outcome == "continuous") {
          eff <- continuous_effects[continuous_effects$design_id == design$design_id, ]
          ni_margin <- 0.20
          eq_lower <- -0.45
          eq_upper <- 0.45
          true_effect <- switch(
            role$scenario_role,
            type1 = if (role$objective == "superiority") 0 else -ni_margin,
            power = eff[[role$objective]],
            type1_lower = eq_lower,
            type1_upper = eq_upper
          )
        } else {
          eff <- binary_effects[binary_effects$design_id == design$design_id, ]
          ni_margin <- 0.10
          eq_lower <- eff$equivalence_lower
          eq_upper <- eff$equivalence_upper
          true_effect <- switch(
            role$scenario_role,
            type1 = if (role$objective == "superiority") 0 else -ni_margin,
            power = eff[[role$objective]],
            type1_lower = eq_lower,
            type1_upper = eq_upper
          )
        }

        pos <- pos + 1L
        rows[[pos]] <- data.frame(
          outcome = outcome,
          design_id = design$design_id,
          factor_count = design$factor_count,
          target_per_group = design$target_per_group,
          total_n = design$total_n,
          objective = role$objective,
          scenario_role = role$scenario_role,
          true_effect = as.numeric(true_effect),
          alpha = role$alpha,
          noninferiority_margin = ni_margin,
          equivalence_lower = eq_lower,
          equivalence_upper = eq_upper,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  out <- do.call(rbind, rows)
  out$scenario_id <- seq_len(nrow(out))
  out$scenario_label <- paste0(
    substr(out$outcome, 1L, 1L),
    "_K", out$factor_count,
    "_npg", out$target_per_group,
    "_", out$objective,
    "_", out$scenario_role
  )
  if (nrow(out) != 56L) stop("Internal error: the full grid must contain 56 scenarios.")
  if (length(SCENARIO_FILTER)) {
    invalid <- setdiff(SCENARIO_FILTER, out$scenario_id)
    if (length(invalid)) stop("Unknown scenario IDs: ", paste(invalid, collapse = ","))
    out <- out[out$scenario_id %in% SCENARIO_FILTER, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

SCENARIOS <- make_scenarios()
DESIGNS <- make_design_table()

scenario_specification <- function(scenario) {
  if (scenario$objective == "superiority") {
    return(list(boundaries = 0, alternatives = "two.sided"))
  }
  if (scenario$objective == "noninferiority") {
    return(list(
      boundaries = -abs(scenario$noninferiority_margin),
      alternatives = "greater"
    ))
  }
  list(
    boundaries = c(scenario$equivalence_lower, scenario$equivalence_upper),
    alternatives = c("greater", "less")
  )
}

safe_write_csv(SCENARIOS, file.path(CONFIG_DIR, "full_scenario_grid.csv"))
safe_write_csv(DESIGNS, file.path(CONFIG_DIR, "design_grid.csv"))

# -----------------------------------------------------------------------------
# Pair-path allocation calibration
# -----------------------------------------------------------------------------

prepare_design <- function(design_row) {
  K <- as.integer(design_row$factor_count)
  patterns <- all_binary_patterns(K)
  factor_prob <- FACTOR_PREVALENCE_MASTER[seq_len(K)]
  profile_prob <- profile_prob_independent(factor_prob, patterns)
  list(
    design_id = as.integer(design_row$design_id),
    factor_count = K,
    target_per_group = as.integer(design_row$target_per_group),
    total_n = as.integer(design_row$total_n),
    patterns = patterns,
    factor_prob = factor_prob,
    profile_prob = profile_prob,
    weights = rep(1, K + 1L)
  )
}

calibration_file <- function(design) {
  ptag <- paste(formatC(design$factor_prob, format = "f", digits = 3), collapse = "-")
  file.path(
    CALIBRATION_DIR,
    sprintf(
      "pair_cal_K%d_n%d_p%s_pbc%.3f_B%d_seed%d.rds",
      design$factor_count,
      design$total_n,
      ptag,
      P_BIASED_COIN,
      N_CALIBRATION,
      calibration_seed(design$design_id)
    )
  )
}

CALIBRATION_MEMORY <- new.env(parent = emptyenv())

get_pair_calibration <- function(design) {
  key <- paste0("D", design$design_id)
  if (exists(key, envir = CALIBRATION_MEMORY, inherits = FALSE)) {
    return(get(key, envir = CALIBRATION_MEMORY, inherits = FALSE))
  }

  path <- calibration_file(design)
  if (file.exists(path)) {
    calibration <- readRDS(path)
    assign(key, calibration, envir = CALIBRATION_MEMORY)
    return(calibration)
  }

  lock <- paste0(path, ".lock")
  have_lock <- dir.create(lock, showWarnings = FALSE)
  if (!have_lock) {
    for (attempt in seq_len(1440L)) {
      if (file.exists(path)) {
        calibration <- readRDS(path)
        assign(key, calibration, envir = CALIBRATION_MEMORY)
        return(calibration)
      }
      Sys.sleep(5)
      have_lock <- dir.create(lock, showWarnings = FALSE)
      if (have_lock) break
    }
  }
  if (!have_lock) stop("Could not acquire calibration lock: ", lock, call. = FALSE)
  on.exit(unlink(lock, recursive = TRUE, force = TRUE), add = TRUE)

  if (file.exists(path)) {
    calibration <- readRDS(path)
  } else {
    message(
      "Creating three-copy calibration for design ", design$design_id,
      ": K=", design$factor_count,
      ", n=", design$total_n,
      ", B=", N_CALIBRATION
    )
    calibration <- calibrate_pair_path_design_R(
      B0 = N_CALIBRATION,
      n = design$total_n,
      patterns = design$patterns,
      profile_prob = design$profile_prob,
      pbc = P_BIASED_COIN,
      weights_ = design$weights,
      seed = calibration_seed(design$design_id),
      batch_size = CALIBRATION_BATCH,
      progress = TRUE
    )
    tmp <- tempfile(pattern = "pair_cal_", tmpdir = CALIBRATION_DIR, fileext = ".rds")
    saveRDS(calibration, tmp)
    if (!file.rename(tmp, path)) {
      if (!file.copy(tmp, path, overwrite = TRUE)) {
        unlink(tmp)
        stop("Could not save calibration: ", path, call. = FALSE)
      }
      unlink(tmp)
    }
  }

  assign(key, calibration, envir = CALIBRATION_MEMORY)
  calibration
}

if (RUN_MODE == "calibrate") {
  for (i in seq_len(nrow(DESIGNS))) {
    design <- prepare_design(DESIGNS[i, ])
    invisible(get_pair_calibration(design))
  }
  message("All four design calibrations are complete.")
  quit(save = "no", status = 0L)
}

# -----------------------------------------------------------------------------
# Outcome models and trial generation
# -----------------------------------------------------------------------------

prepare_scenario <- function(scenario) {
  design <- prepare_design(DESIGNS[DESIGNS$design_id == scenario$design_id, ])
  if (scenario$outcome == "continuous") {
    model <- list(
      type = "continuous_benchmark",
      mu0 = pattern_baseline_mean(
        patterns = design$patterns,
        profile_prob = design$profile_prob,
        beta_master = CONTINUOUS_BETA_MASTER,
        interaction = CONTINUOUS_INTERACTION
      ),
      true_effect = scenario$true_effect,
      outcome_sd = CONTINUOUS_OUTCOME_SD
    )
  } else {
    model <- calibrate_common_log_odds_model(
      patterns = design$patterns,
      profile_prob = design$profile_prob,
      target_control_risk = TARGET_CONTROL_RISK,
      target_risk_difference = scenario$true_effect,
      beta_master = BINARY_BETA_MASTER,
      interaction = BINARY_INTERACTION
    )
  }
  list(scenario = scenario, design = design, model = model)
}

simulate_trial <- function(prepared, replicate_id) {
  scenario <- prepared$scenario
  design <- prepared$design
  model <- prepared$model

  set.seed(seed_value(scenario$scenario_id, replicate_id, 1L))
  profile <- generate_profile_sequence(
    n = design$total_n,
    patterns = design$patterns,
    profile_prob = design$profile_prob
  )
  z <- ps_assign_R(
    X = profile$X,
    pbc = P_BIASED_COIN,
    weights_ = design$weights,
    seed = seed_value(scenario$scenario_id, replicate_id, 2L)
  )
  A <- as.integer((z + 1L) / 2L)

  set.seed(seed_value(scenario$scenario_id, replicate_id, 3L))
  if (scenario$outcome == "continuous") {
    y0 <- model$mu0[profile$id] + rnorm(design$total_n, sd = model$outcome_sd)
    y <- y0 + A * model$true_effect
  } else {
    p <- ifelse(A == 1L, model$p1[profile$id], model$p0[profile$id])
    y <- rbinom(design$total_n, size = 1L, prob = p)
  }

  list(
    X = profile$X,
    profile_id = profile$id,
    z = as.integer(z),
    A = A,
    y = y
  )
}

boundary_deviation <- function(prepared, boundary) {
  scenario <- prepared$scenario
  model <- prepared$model
  J <- nrow(prepared$design$patterns)
  if (scenario$outcome == "continuous") {
    return(rep(scenario$true_effect - boundary, J))
  }
  as.numeric(model$p1 - model$p0 - boundary)
}

# -----------------------------------------------------------------------------
# P-values and the manuscript safeguard
# -----------------------------------------------------------------------------

gaussian_pvalue_local <- function(statistic, variance, alternative) {
  z <- statistic / sqrt(max(variance, .Machine$double.eps))
  switch(
    alternative,
    two.sided = 2 * pnorm(-abs(z)),
    greater = pnorm(z, lower.tail = FALSE),
    less = pnorm(z),
    stop("Unknown alternative: ", alternative, call. = FALSE)
  )
}

lattice_normal_pvalue <- function(statistic,
                                  total_success,
                                  variance,
                                  treated_count_prob,
                                  alternative,
                                  tolerance = 1e-12) {
  n <- length(treated_count_prob) - 1L
  S <- as.integer(total_success)
  if (S < 0L || S > n) stop("total_success must be between 0 and n.")
  if (S == 0L || S == n) return(1)

  wk <- as.numeric(treated_count_prob)
  wk <- wk / sum(wk)
  sdv <- sqrt(max(variance, .Machine$double.eps))
  ans <- 0
  for (k in 0:n) {
    if (wk[k + 1L] <= 0) next
    lo <- max(0L, k - (n - S))
    hi <- min(k, S)
    rr <- lo:hi
    support <- rr - k * S / n
    mass <- pnorm((support + 0.5) / sdv) - pnorm((support - 0.5) / sdv)
    normalizer <- sum(mass)
    if (!is.finite(normalizer) || normalizer <= 0) next
    keep <- switch(
      alternative,
      two.sided = abs(support) + tolerance >= abs(statistic),
      greater = support + tolerance >= statistic,
      less = support - tolerance <= statistic
    )
    ans <- ans + wk[k + 1L] * sum(mass[keep]) / normalizer
  }
  min(max(ans, 0), 1)
}

randomization_variance_with_safeguard <- function(sampling_result, d, calibration) {
  d <- as.numeric(d)
  if (length(d) != calibration$J) stop("d has the wrong length.", call. = FALSE)
  counts <- sampling_result$decomposition$counts
  pi_hat <- counts / sum(counts)
  correction_matrix <- calibration$psi - diag(pi_hat, nrow = length(pi_hat))
  scalar_gap <- drop(crossprod(d, correction_matrix %*% d))
  correction <- calibration$n * scalar_gap / 16
  raw <- sampling_result$variance + correction
  epsilon_n <- calibration$n ^ (-EPSILON_EXPONENT)
  lower_bound <- epsilon_n * sampling_result$variance
  variance <- if (is.finite(raw) && raw >= lower_bound) raw else lower_bound
  list(
    variance = variance,
    variance_raw = raw,
    correction = correction,
    scalar_gap = scalar_gap,
    pi_hat = pi_hat,
    epsilon_n = epsilon_n,
    lower_bound = lower_bound,
    safeguarded = !is.finite(raw) || raw < lower_bound
  )
}

combine_objective <- function(component_p, objective) {
  if (objective == "equivalence") max(component_p) else component_p[1L]
}

rt_component_pvalue <- function(rt, column, alternative) {
  switch(
    alternative,
    two.sided = rt$two_sided[column],
    greater = rt$greater[column],
    less = rt$less[column]
  )
}

# -----------------------------------------------------------------------------
# Analysis of one outer trial
# -----------------------------------------------------------------------------

analyze_trial <- function(prepared, trial, calibration, replicate_id) {
  scenario <- prepared$scenario
  spec <- scenario_specification(scenario)
  nb <- length(spec$boundaries)

  scores_u <- matrix(NA_real_, nrow = scenario$total_n, ncol = nb)
  scores_a <- matrix(NA_real_, nrow = scenario$total_n, ncol = nb)
  for (j in seq_len(nb)) {
    ss <- make_boundary_scores(
      y = trial$y,
      A = trial$A,
      X = trial$X,
      boundary = spec$boundaries[j]
    )
    scores_u[, j] <- ss$unadjusted
    scores_a[, j] <- ss$adjusted
  }

  score_matrix <- cbind(scores_u, scores_a)
  rt <- rt_pvalues_R(
    X = trial$X,
    scores = score_matrix,
    z_obs = trial$z,
    B = N_RERANDOMIZATIONS,
    pbc = P_BIASED_COIN,
    weights_ = prepared$design$weights,
    seed = seed_value(scenario$scenario_id, replicate_id, 4L),
    return_randomization_moments = TRUE
  )

  p_s_u <- p_r_u <- p_rt_u <- numeric(nb)
  p_s_a <- p_r_a <- p_rt_a <- numeric(nb)
  vs_u <- vr_u <- vrt_u <- ratio_rs_u <- ratio_rrt_u <- numeric(nb)
  vs_a <- vr_a <- vrt_a <- ratio_rs_a <- ratio_rrt_a <- numeric(nb)
  safeguard_u <- safeguard_a <- logical(nb)
  d_mean <- d_max_abs <- numeric(nb)

  for (j in seq_len(nb)) {
    boundary <- spec$boundaries[j]
    alternative <- spec$alternatives[j]
    d <- boundary_deviation(prepared, boundary)
    d_mean[j] <- sum(prepared$design$profile_prob * d)
    d_max_abs[j] <- max(abs(d))

    stat_u <- 0.5 * sum(trial$z * scores_u[, j])
    s_u <- siga_sampling_variance(scores_u[, j], trial$X, calibration)
    r_u <- randomization_variance_with_safeguard(s_u, d, calibration)

    use_lattice <- scenario$outcome == "binary" &&
      scenario$objective == "superiority" &&
      abs(boundary) < 1e-14
    if (use_lattice) {
      p_s_u[j] <- lattice_normal_pvalue(
        statistic = stat_u,
        total_success = sum(trial$y),
        variance = s_u$variance,
        treated_count_prob = calibration$treated_count_prob,
        alternative = alternative
      )
      p_r_u[j] <- lattice_normal_pvalue(
        statistic = stat_u,
        total_success = sum(trial$y),
        variance = r_u$variance,
        treated_count_prob = calibration$treated_count_prob,
        alternative = alternative
      )
    } else {
      p_s_u[j] <- gaussian_pvalue_local(stat_u, s_u$variance, alternative)
      p_r_u[j] <- gaussian_pvalue_local(stat_u, r_u$variance, alternative)
    }
    p_rt_u[j] <- rt_component_pvalue(rt, j, alternative)
    vs_u[j] <- s_u$variance
    vr_u[j] <- r_u$variance
    vrt_u[j] <- rt$randomization_variance[j]
    ratio_rs_u[j] <- vr_u[j] / vs_u[j]
    ratio_rrt_u[j] <- vr_u[j] / vrt_u[j]
    safeguard_u[j] <- r_u$safeguarded

    stat_a <- 0.5 * sum(trial$z * scores_a[, j])
    s_a <- siga_sampling_variance(scores_a[, j], trial$X, calibration)
    r_a <- randomization_variance_with_safeguard(s_a, d, calibration)
    p_s_a[j] <- gaussian_pvalue_local(stat_a, s_a$variance, alternative)
    p_r_a[j] <- gaussian_pvalue_local(stat_a, r_a$variance, alternative)
    p_rt_a[j] <- rt_component_pvalue(rt, nb + j, alternative)
    vs_a[j] <- s_a$variance
    vr_a[j] <- r_a$variance
    vrt_a[j] <- rt$randomization_variance[nb + j]
    ratio_rs_a[j] <- vr_a[j] / vs_a[j]
    ratio_rrt_a[j] <- vr_a[j] / vrt_a[j]
    safeguard_a[j] <- r_a$safeguarded
  }

  pad2 <- function(x) {
    x <- as.vector(x)
    if (length(x) < 2L) x <- c(x, NA)
    x[1:2]
  }

  list(
    p_s_unadjusted = combine_objective(p_s_u, scenario$objective),
    p_r_unadjusted = combine_objective(p_r_u, scenario$objective),
    p_rt_unadjusted = combine_objective(p_rt_u, scenario$objective),
    p_s_adjusted = combine_objective(p_s_a, scenario$objective),
    p_r_adjusted = combine_objective(p_r_a, scenario$objective),
    p_rt_adjusted = combine_objective(p_rt_a, scenario$objective),
    component = list(
      p_s_u = pad2(p_s_u), p_r_u = pad2(p_r_u), p_rt_u = pad2(p_rt_u),
      p_s_a = pad2(p_s_a), p_r_a = pad2(p_r_a), p_rt_a = pad2(p_rt_a),
      vs_u = pad2(vs_u), vr_u = pad2(vr_u), vrt_u = pad2(vrt_u),
      vs_a = pad2(vs_a), vr_a = pad2(vr_a), vrt_a = pad2(vrt_a),
      ratio_rs_u = pad2(ratio_rs_u), ratio_rrt_u = pad2(ratio_rrt_u),
      ratio_rs_a = pad2(ratio_rs_a), ratio_rrt_a = pad2(ratio_rrt_a),
      safeguard_u = pad2(safeguard_u), safeguard_a = pad2(safeguard_a),
      d_mean = pad2(d_mean), d_max_abs = pad2(d_max_abs)
    )
  )
}

run_replicate <- function(prepared, replicate_id, calibration) {
  scenario <- prepared$scenario
  generation_start <- proc.time()[3L]
  trial <- simulate_trial(prepared, replicate_id)
  generation_seconds <- proc.time()[3L] - generation_start

  analysis_start <- proc.time()[3L]
  ans <- analyze_trial(prepared, trial, calibration, replicate_id)
  analysis_seconds <- proc.time()[3L] - analysis_start
  cpt <- ans$component

  data.frame(
    scenario_id = scenario$scenario_id,
    scenario_label = scenario$scenario_label,
    replicate = replicate_id,
    outcome = scenario$outcome,
    design_id = scenario$design_id,
    factor_count = scenario$factor_count,
    target_per_group = scenario$target_per_group,
    total_n = scenario$total_n,
    objective = scenario$objective,
    scenario_role = scenario$scenario_role,
    true_effect = scenario$true_effect,
    alpha = scenario$alpha,
    p_s_unadjusted = ans$p_s_unadjusted,
    p_r_unadjusted = ans$p_r_unadjusted,
    p_rt_unadjusted = ans$p_rt_unadjusted,
    p_s_adjusted = ans$p_s_adjusted,
    p_r_adjusted = ans$p_r_adjusted,
    p_rt_adjusted = ans$p_rt_adjusted,
    reject_s_unadjusted = ans$p_s_unadjusted <= scenario$alpha,
    reject_r_unadjusted = ans$p_r_unadjusted <= scenario$alpha,
    reject_rt_unadjusted = ans$p_rt_unadjusted <= scenario$alpha,
    reject_s_adjusted = ans$p_s_adjusted <= scenario$alpha,
    reject_r_adjusted = ans$p_r_adjusted <= scenario$alpha,
    reject_rt_adjusted = ans$p_rt_adjusted <= scenario$alpha,
    p_s_unadjusted_1 = cpt$p_s_u[1L],
    p_s_unadjusted_2 = cpt$p_s_u[2L],
    p_r_unadjusted_1 = cpt$p_r_u[1L],
    p_r_unadjusted_2 = cpt$p_r_u[2L],
    p_rt_unadjusted_1 = cpt$p_rt_u[1L],
    p_rt_unadjusted_2 = cpt$p_rt_u[2L],
    p_s_adjusted_1 = cpt$p_s_a[1L],
    p_s_adjusted_2 = cpt$p_s_a[2L],
    p_r_adjusted_1 = cpt$p_r_a[1L],
    p_r_adjusted_2 = cpt$p_r_a[2L],
    p_rt_adjusted_1 = cpt$p_rt_a[1L],
    p_rt_adjusted_2 = cpt$p_rt_a[2L],
    variance_s_unadjusted_1 = cpt$vs_u[1L],
    variance_s_unadjusted_2 = cpt$vs_u[2L],
    variance_r_unadjusted_1 = cpt$vr_u[1L],
    variance_r_unadjusted_2 = cpt$vr_u[2L],
    variance_rt_unadjusted_1 = cpt$vrt_u[1L],
    variance_rt_unadjusted_2 = cpt$vrt_u[2L],
    variance_s_adjusted_1 = cpt$vs_a[1L],
    variance_s_adjusted_2 = cpt$vs_a[2L],
    variance_r_adjusted_1 = cpt$vr_a[1L],
    variance_r_adjusted_2 = cpt$vr_a[2L],
    variance_rt_adjusted_1 = cpt$vrt_a[1L],
    variance_rt_adjusted_2 = cpt$vrt_a[2L],
    ratio_r_over_s_unadjusted_1 = cpt$ratio_rs_u[1L],
    ratio_r_over_s_unadjusted_2 = cpt$ratio_rs_u[2L],
    ratio_r_over_rt_unadjusted_1 = cpt$ratio_rrt_u[1L],
    ratio_r_over_rt_unadjusted_2 = cpt$ratio_rrt_u[2L],
    ratio_r_over_s_adjusted_1 = cpt$ratio_rs_a[1L],
    ratio_r_over_s_adjusted_2 = cpt$ratio_rs_a[2L],
    ratio_r_over_rt_adjusted_1 = cpt$ratio_rrt_a[1L],
    ratio_r_over_rt_adjusted_2 = cpt$ratio_rrt_a[2L],
    safeguard_unadjusted_1 = cpt$safeguard_u[1L],
    safeguard_unadjusted_2 = cpt$safeguard_u[2L],
    safeguard_adjusted_1 = cpt$safeguard_a[1L],
    safeguard_adjusted_2 = cpt$safeguard_a[2L],
    d_weighted_mean_1 = cpt$d_mean[1L],
    d_weighted_mean_2 = cpt$d_mean[2L],
    d_max_abs_1 = cpt$d_max_abs[1L],
    d_max_abs_2 = cpt$d_max_abs[2L],
    data_generation_seconds = generation_seconds,
    complete_analysis_seconds = analysis_seconds,
    configured_outer_trials = N_OUTER,
    rerandomizations_per_trial = N_RERANDOMIZATIONS,
    pair_calibration_replicates = N_CALIBRATION,
    epsilon_exponent = EPSILON_EXPONENT,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Sharded checkpointed execution
# -----------------------------------------------------------------------------

scenario_paths <- function(scenario) {
  directory <- ensure_directory(file.path(
    SCENARIO_DIR,
    sprintf("scenario_%02d_%s", scenario$scenario_id, scenario$scenario_label)
  ))
  list(
    directory = directory,
    detail = file.path(
      directory,
      sprintf("shard_%04d_of_%04d.csv", SHARD_ID, N_SHARDS)
    ),
    timing = file.path(
      directory,
      sprintf("shard_%04d_of_%04d_timing.csv", SHARD_ID, N_SHARDS)
    )
  )
}

read_completed_replicates <- function(path) {
  if (!file.exists(path)) return(integer())
  x <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(x) || !"replicate" %in% names(x)) {
    stop("Could not read checkpoint: ", path, call. = FALSE)
  }
  unique(as.integer(x$replicate))
}

write_timing <- function(path, data) {
  safe_write_csv(data, path)
}

run_scenario_shard <- function(scenario) {
  prepared <- prepare_scenario(scenario)
  calibration <- get_pair_calibration(prepared$design)
  paths <- scenario_paths(scenario)
  assigned <- seq.int(from = SHARD_ID, to = N_OUTER, by = N_SHARDS)
  completed <- read_completed_replicates(paths$detail)
  pending <- setdiff(assigned, completed)

  message(
    "Scenario ", scenario$scenario_id, "/56 [", scenario$scenario_label, "]",
    "; shard ", SHARD_ID, "/", N_SHARDS,
    "; assigned=", length(assigned),
    "; completed=", length(intersect(assigned, completed)),
    "; pending=", length(pending)
  )
  if (!length(pending)) return(invisible(NULL))

  start <- proc.time()[3L]
  for (from in seq.int(1L, length(pending), by = OUTER_BATCH)) {
    ids <- pending[from:min(from + OUTER_BATCH - 1L, length(pending))]
    batch <- do.call(rbind, lapply(ids, function(id) {
      run_replicate(prepared, id, calibration)
    }))
    append_csv(batch, paths$detail)
    done <- min(from + length(ids) - 1L, length(pending))
    if (done %% max(OUTER_BATCH, 100L) == 0L || done == length(pending)) {
      message(
        "  scenario ", scenario$scenario_id,
        ", shard ", SHARD_ID,
        ": ", done, "/", length(pending),
        " pending replicates completed"
      )
    }
  }

  timing <- data.frame(
    scenario_id = scenario$scenario_id,
    scenario_label = scenario$scenario_label,
    shard_id = SHARD_ID,
    n_shards = N_SHARDS,
    assigned_replicates = length(assigned),
    newly_completed_replicates = length(pending),
    elapsed_seconds = proc.time()[3L] - start,
    completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    stringsAsFactors = FALSE
  )
  write_timing(paths$timing, timing)
  invisible(NULL)
}

message("Comprehensive SIGA-R benchmark")
message("  output base               : ", OUTPUT_ROOT)
message("  run directory             : ", RUN_DIR)
message("  outer trials/scenario     : ", N_OUTER)
message("  RT paths/outer trial      : ", N_RERANDOMIZATIONS)
message("  three-copy calibrations   : ", N_CALIBRATION)
message("  biased-coin probability   : ", P_BIASED_COIN)
message("  epsilon_n                 : n^(-", EPSILON_EXPONENT, ")")
message("  shard                     : ", SHARD_ID, "/", N_SHARDS)
message("  selected scenarios        : ", paste(SCENARIOS$scenario_id, collapse = ","))

for (i in seq_len(nrow(SCENARIOS))) {
  run_scenario_shard(SCENARIOS[i, ])
}

message("Shard completed successfully.")
