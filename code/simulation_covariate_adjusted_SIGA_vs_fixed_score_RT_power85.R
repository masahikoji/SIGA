#!/usr/bin/env Rscript

# =============================================================================
# Full operating-characteristic comparison under biased-coin minimization
# METHOD: covariate-adjusted SIGA approximation versus the same fixed-score randomization test
# =============================================================================
#
# Scientific objective
# --------------------
# Compare, on the same independently generated trials:
#   1. the proposed SIGA allocation-calibrated approximation; and
#   2. an ordinary conditional randomization test that regenerates allocation
#      paths from the prespecified minimization algorithm.
#
# IMPORTANT: the comparator is the fixed-score conditional randomization test
# targeted by the theory and manuscript. The observed and rerandomized paths use
# exactly T = sum_i (A_i - 1/2) r_i, with the boundary-specific score r held
# fixed. The comparator is not a mean-difference test with a varying treatment
# count and is not a refitted ANCOVA t test.
#
# The proposed method uses the fixed score T = sum_i (A_i - 1/2) r_i and an
# allocation-only calibration of the joint-stratum imbalance covariance. The
# calibration is generated once per (n, factor count, prevalence, p_bc) design,
# cached, and reused across all testing objectives and outer trials.
#
# The program remains pure base R. No Rcpp or external package is required.
#
# Examples
# --------
# Smoke test:
#   PWRT_N_OUTER=20 PWRT_N_RERAND=99 PWRT_N_CALIBRATION=1000 Rscript simulation_covariate_adjusted_SIGA_vs_fixed_score_RT_power85.R
#
# Full requested settings are the defaults:
#   PWRT_N_OUTER=100000, PWRT_N_RERAND=4999, PWRT_N_CALIBRATION=100000.
#
# Existing sharding, checkpointing, resumption, objective filtering, and
# aggregation interfaces are retained. New output directories are used so that
# old pathwise results cannot be mixed with the corrected analysis.
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

ADJUSTED <- TRUE
METHOD_ID <- "covariate_adjusted_siga_vs_fixed_score_rt_power85"
METHOD_TITLE <- "covariate-adjusted SIGA approximation versus fixed-score conditional randomization test with power targeted at 80--90 percent"
PROPOSED_LABEL <- "SIGA"
REFERENCE_LABEL <- "Fixed-score randomization test"

# =============================================================================
# Environment and paths
# =============================================================================

script_directory <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    path <- sub("^--file=", "", file_arg[1L])
    return(dirname(normalizePath(path, winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

env_integer <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.integer(default))
  out <- suppressWarnings(as.integer(value))
  if (is.na(out)) stop(name, " must be an integer.", call. = FALSE)
  out
}

env_numeric <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.numeric(default))
  out <- suppressWarnings(as.numeric(value))
  if (!is.finite(out)) stop(name, " must be numeric.", call. = FALSE)
  out
}

env_character_vector <- function(name, default = character()) {
  value <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) return(default)
  out <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  out[nzchar(out)]
}

BASE_DIR <- script_directory()

# The project directory is fixed explicitly so that moving or launching this
# script from another folder cannot redirect a long simulation to an unrelated
# project. PWRT_PROJECT_DIR and PWRT_OUTPUT_DIR may still override these defaults.
PROJECT_DIR <- path.expand(Sys.getenv(
  "PWRT_PROJECT_DIR",
  unset = "~/Dropbox/Tex/03_Post_PhD/029_exact_minimization/02_program"
))
OUTPUT_DIR <- path.expand(Sys.getenv(
  "PWRT_OUTPUT_DIR",
  unset = file.path(
    PROJECT_DIR,
    paste0(METHOD_ID, "_output")
  )
))

ensure_directory <- function(path, attempts = 5L, wait_seconds = 2) {
  path <- path.expand(path)
  for (attempt in seq_len(attempts)) {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }
    if (dir.exists(path)) {
      return(normalizePath(path, winslash = "/", mustWork = TRUE))
    }
    if (attempt < attempts) Sys.sleep(wait_seconds)
  }
  stop("Could not create or access directory: ", path, call. = FALSE)
}

PROJECT_DIR <- ensure_directory(PROJECT_DIR)
OUTPUT_DIR <- ensure_directory(OUTPUT_DIR)
CHUNK_DIR <- ensure_directory(file.path(OUTPUT_DIR, "scenario_shards"))
CALIBRATION_DIR <- ensure_directory(file.path(
  PROJECT_DIR, "siga_allocation_calibration_cache"
))

# =============================================================================
# Requested full-run settings
# =============================================================================

RUN_MODE <- tolower(trimws(Sys.getenv("PWRT_MODE", unset = "run")))
if (!RUN_MODE %in% c("run", "aggregate")) {
  stop("PWRT_MODE must be 'run' or 'aggregate'.", call. = FALSE)
}

N_OUTER <- env_integer("PWRT_N_OUTER", 100000L)
N_RERANDOMIZATIONS <- env_integer("PWRT_N_RERAND", 4999L)
N_SHARDS <- env_integer("PWRT_N_SHARDS", 1L)
SHARD_ID <- env_integer("PWRT_SHARD_ID", 1L)
OUTER_BATCH_SIZE <- env_integer("PWRT_OUTER_BATCH", 10L)
P_BIASED_COIN <- env_numeric("PWRT_P_BIASED_COIN", 0.80)
BASE_SEED <- env_integer("PWRT_SEED", 20260725L)
N_CALIBRATION <- env_integer("PWRT_N_CALIBRATION", 100000L)
CALIBRATION_BATCH_SIZE <- env_integer("PWRT_CALIBRATION_BATCH", 5000L)

physical_cores <- parallel::detectCores(logical = FALSE)
if (is.na(physical_cores)) physical_cores <- parallel::detectCores(logical = TRUE)
if (is.na(physical_cores)) physical_cores <- 1L
DEFAULT_CORES <- if (.Platform$OS.type == "windows") 1L else max(1L, min(8L, physical_cores - 1L))
N_CORES <- env_integer("PWRT_CORES", DEFAULT_CORES)
if (.Platform$OS.type == "windows") N_CORES <- 1L
USE_PARALLEL <- N_CORES > 1L && .Platform$OS.type == "unix"

OBJECTIVE_FILTER <- env_character_vector(
  "PWRT_OBJECTIVES",
  c("superiority", "noninferiority", "equivalence")
)
SCENARIO_ID_FILTER_RAW <- env_character_vector("PWRT_SCENARIO_IDS", character())
SCENARIO_ID_FILTER <- if (length(SCENARIO_ID_FILTER_RAW)) {
  as.integer(SCENARIO_ID_FILTER_RAW)
} else {
  integer()
}
if (anyNA(SCENARIO_ID_FILTER)) {
  stop("PWRT_SCENARIO_IDS must contain comma-separated integers.", call. = FALSE)
}

# Testing objectives and outcome model.
ALPHA_SUPERIORITY <- 0.05
ALPHA_NONINFERIORITY <- 0.025
ALPHA_EQUIVALENCE <- 0.05
NONINFERIORITY_MARGIN <- 0.20
EQUIVALENCE_LIMITS <- c(-0.45, 0.45)
POWER_EFFECT_SUPERIORITY <- 0.20
POWER_EFFECT_NONINFERIORITY <- 0.00
POWER_EFFECT_EQUIVALENCE <- 0.00
OUTCOME_SD <- 1.00
INTERACTION_COEFFICIENT <- 0.25

FACTOR_PREVALENCE_MASTER <- c(0.50, 0.40, 0.30, 0.20, 0.10)
FACTOR_COEFFICIENT_MASTER <- c(0.50, 0.40, 0.30, 0.20, 0.10)

# =============================================================================
# Validation
# =============================================================================

validate_settings <- function() {
  stopifnot(
    N_OUTER >= 1L,
    N_RERANDOMIZATIONS >= 1L,
    N_CALIBRATION >= 2L,
    CALIBRATION_BATCH_SIZE >= 1L,
    N_SHARDS >= 1L,
    SHARD_ID >= 1L,
    SHARD_ID <= N_SHARDS,
    OUTER_BATCH_SIZE >= 1L,
    P_BIASED_COIN > 0.5,
    P_BIASED_COIN < 1,
    N_CORES >= 1L,
    all(OBJECTIVE_FILTER %in% c("superiority", "noninferiority", "equivalence"))
  )
  if (N_RERANDOMIZATIONS != 4999L) {
    warning(
      "The manuscript design specifies 4,999 rerandomizations; current value is ",
      N_RERANDOMIZATIONS, "."
    )
  }
  invisible(TRUE)
}
validate_settings()

# =============================================================================
# Reproducible seeds and statistical utilities
# =============================================================================

seed_value <- function(scenario_id, replicate_id, stream_id) {
  modulus <- 2147483646
  value <- as.double(BASE_SEED) +
    1000003 * as.double(stream_id) +
    104729 * as.double(scenario_id) +
    1009 * as.double(replicate_id)
  as.integer(value %% modulus + 1)
}

wilson_interval <- function(x, n, conf.level = 0.95) {
  if (!is.finite(n) || n <= 0L) {
    return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  z <- qnorm(1 - (1 - conf.level) / 2)
  p <- x / n
  den <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(
    estimate = p,
    lower = max(0, center - half),
    upper = min(1, center + half)
  )
}

append_csv <- function(data, path, attempts = 5L, wait_seconds = 2) {
  if (!nrow(data)) return(invisible(NULL))

  parent <- dirname(path)
  last_error <- NULL
  for (attempt in seq_len(attempts)) {
    tmp <- NULL
    ok <- tryCatch({
      ensure_directory(parent, attempts = 1L, wait_seconds = 0)
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
        if (!file.append(path, tmp)) {
          stop("file.append returned FALSE")
        }
        unlink(tmp)
      } else {
        if (!file.rename(tmp, path)) {
          if (!file.copy(tmp, path, overwrite = FALSE)) {
            stop("Could not move the first batch into place")
          }
          unlink(tmp)
        }
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
    "Failed to append checkpoint after ", attempts, " attempts: ", path,
    if (!is.null(last_error)) paste0("; last error: ", last_error) else "",
    call. = FALSE
  )
}

format_seconds <- function(x) {
  if (!is.finite(x)) return("NA")
  sprintf("%.2f min", x / 60)
}

# =============================================================================
# Scenario grid
# =============================================================================

make_scenarios <- function() {
  # Design-specific alternatives were selected from a separate pilot to
  # target approximately 85% power (acceptable band 80%--90%) under the
  # fixed-score reference randomization test. For K=2 and n/group=100,
  # equivalence power is maximized near true effect 0 and was already in
  # the target band, so the centered alternative is retained.
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
    true_effect_code = c(
      "superiority_type1", "superiority_power",
      "noninferiority_type1", "noninferiority_power",
      "equivalence_type1_lower", "equivalence_type1_upper", "equivalence_power"
    ),
    alpha = c(
      ALPHA_SUPERIORITY, ALPHA_SUPERIORITY,
      ALPHA_NONINFERIORITY, ALPHA_NONINFERIORITY,
      ALPHA_EQUIVALENCE, ALPHA_EQUIVALENCE, ALPHA_EQUIVALENCE
    ),
    stringsAsFactors = FALSE
  )

  rows <- vector("list", nrow(designs) * nrow(roles))
  pos <- 0L
  for (d in seq_len(nrow(designs))) {
    for (r in seq_len(nrow(roles))) {
      pos <- pos + 1L
      rows[[pos]] <- data.frame(
        factor_count = designs$factor_count[d],
        target_per_group = designs$target_per_group[d],
        total_n = 2L * designs$target_per_group[d],
        objective = roles$objective[r],
        scenario_role = roles$scenario_role[r],
        true_effect = switch(
          roles$true_effect_code[r],
          superiority_type1 = 0,
          superiority_power = designs$superiority_power_effect[d],
          noninferiority_type1 = -NONINFERIORITY_MARGIN,
          noninferiority_power = designs$noninferiority_power_effect[d],
          equivalence_type1_lower = EQUIVALENCE_LIMITS[1L],
          equivalence_type1_upper = EQUIVALENCE_LIMITS[2L],
          equivalence_power = designs$equivalence_power_effect[d]
        ),
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
  out <- out[out$objective %in% OBJECTIVE_FILTER, , drop = FALSE]
  if (length(SCENARIO_ID_FILTER)) {
    out <- out[out$scenario_id %in% SCENARIO_ID_FILTER, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

SCENARIOS <- make_scenarios()
if (!nrow(SCENARIOS)) stop("No scenarios remain after filtering.", call. = FALSE)

scenario_boundaries <- function(scenario) {
  if (scenario$objective == "superiority") return(c(zero = 0))
  if (scenario$objective == "noninferiority") {
    return(c(lower = -NONINFERIORITY_MARGIN))
  }
  c(lower = EQUIVALENCE_LIMITS[1L], upper = EQUIVALENCE_LIMITS[2L])
}

# =============================================================================
# Covariate and outcome generation
# =============================================================================

generate_factors <- function(n, factor_count) {
  prevalence <- FACTOR_PREVALENCE_MASTER[seq_len(factor_count)]
  X <- matrix(0L, nrow = n, ncol = factor_count)
  for (j in seq_len(factor_count)) {
    X[, j] <- rbinom(n, size = 1L, prob = prevalence[j])
  }
  colnames(X) <- paste0("X", seq_len(factor_count))
  X
}

baseline_mean <- function(X) {
  factor_count <- ncol(X)
  prevalence <- FACTOR_PREVALENCE_MASTER[seq_len(factor_count)]
  coefficients <- FACTOR_COEFFICIENT_MASTER[seq_len(factor_count)]
  centered <- sweep(X, 2L, prevalence, "-")
  mu <- drop(centered %*% coefficients)
  if (factor_count >= 2L) {
    mu <- mu + INTERACTION_COEFFICIENT *
      (X[, 1L] * X[, 2L] - prevalence[1L] * prevalence[2L])
  }
  mu
}

# =============================================================================
# General K-factor biased-coin Pocock--Simon minimization
# =============================================================================

allocation_probability <- function(x, d_overall, d_level1) {
  score_if_1 <- abs(d_overall + 1L)
  score_if_0 <- abs(d_overall - 1L)
  for (j in seq_along(x)) {
    level_imbalance <- if (x[j] == 1L) d_level1[j] else d_overall - d_level1[j]
    score_if_1 <- score_if_1 + abs(level_imbalance + 1L)
    score_if_0 <- score_if_0 + abs(level_imbalance - 1L)
  }
  if (score_if_1 < score_if_0) return(P_BIASED_COIN)
  if (score_if_1 > score_if_0) return(1 - P_BIASED_COIN)
  0.5
}

generate_observed_allocation <- function(X) {
  n <- nrow(X)
  factor_count <- ncol(X)
  A <- integer(n)
  q <- numeric(n)
  d_overall <- 0L
  d_level1 <- integer(factor_count)

  for (i in seq_len(n)) {
    q_i <- allocation_probability(X[i, ], d_overall, d_level1)
    A_i <- as.integer(runif(1L) < q_i)
    increment <- 2L * A_i - 1L
    A[i] <- A_i
    q[i] <- q_i
    d_overall <- d_overall + increment
    d_level1 <- d_level1 + increment * X[i, ]
  }
  list(A = A, q = q)
}

simulate_trial <- function(scenario) {
  X <- generate_factors(scenario$total_n, scenario$factor_count)
  allocation <- generate_observed_allocation(X)
  Y0 <- baseline_mean(X) + rnorm(scenario$total_n, sd = OUTCOME_SD)
  Y <- Y0 + allocation$A * scenario$true_effect
  list(X = X, A = allocation$A, q = allocation$q, Y = Y)
}

# =============================================================================
# Allocation-only SIGA calibration
# =============================================================================

joint_stratum_id <- function(X) {
  X <- as.matrix(X)
  powers <- 2L ^ (0:(ncol(X) - 1L))
  as.integer(1L + drop(X %*% powers))
}

nearest_psd <- function(M, tolerance = 1e-12) {
  M <- (M + t(M)) / 2
  ee <- eigen(M, symmetric = TRUE)
  cutoff <- tolerance * max(1, max(abs(ee$values)))
  values <- pmax(ee$values, 0)
  values[values < cutoff] <- 0
  out <- ee$vectors %*% (values * t(ee$vectors))
  (out + t(out)) / 2
}

calibration_seed <- function(total_n, factor_count) {
  modulus <- 2147483646
  value <- as.double(BASE_SEED) + 7000003 +
    1009 * as.double(total_n) + 104729 * as.double(factor_count)
  as.integer(value %% modulus + 1)
}

calibration_cache_file <- function(scenario) {
  probs <- FACTOR_PREVALENCE_MASTER[seq_len(scenario$factor_count)]
  prob_tag <- paste(formatC(probs, format = "f", digits = 3), collapse = "-")
  pbc_tag <- formatC(P_BIASED_COIN, format = "f", digits = 3)
  file.path(
    CALIBRATION_DIR,
    sprintf(
      "siga_cal_K%d_n%d_p%s_pbc%s_B%d_seed%d.rds",
      scenario$factor_count, scenario$total_n, prob_tag, pbc_tag,
      N_CALIBRATION, calibration_seed(scenario$total_n, scenario$factor_count)
    )
  )
}

run_allocation_calibration <- function(scenario) {
  n <- scenario$total_n
  K <- scenario$factor_count
  J <- 2L ^ K
  probs <- FACTOR_PREVALENCE_MASTER[seq_len(K)]
  B0 <- N_CALIBRATION
  batch_size <- min(CALIBRATION_BATCH_SIZE, B0)

  sum_u <- numeric(J)
  sum_uu <- matrix(0, nrow = J, ncol = J)
  treated_count <- numeric(n + 1L)
  powers <- 2L ^ (0:(K - 1L))

  set.seed(calibration_seed(n, K))
  start <- proc.time()[3L]
  completed <- 0L

  while (completed < B0) {
    B <- min(batch_size, B0 - completed)
    d_overall <- integer(B)
    d_level1 <- matrix(0L, nrow = B, ncol = K)
    stratum_counts <- matrix(0L, nrow = B, ncol = J)
    stratum_imbalance <- matrix(0L, nrow = B, ncol = J)
    row_index <- seq_len(B)

    for (i in seq_len(n)) {
      x <- matrix(0L, nrow = B, ncol = K)
      for (j in seq_len(K)) {
        x[, j] <- as.integer(runif(B) < probs[j])
      }

      score_if_1 <- abs(d_overall + 1L)
      score_if_0 <- abs(d_overall - 1L)
      for (j in seq_len(K)) {
        level_imbalance <- ifelse(
          x[, j] == 1L,
          d_level1[, j],
          d_overall - d_level1[, j]
        )
        score_if_1 <- score_if_1 + abs(level_imbalance + 1L)
        score_if_0 <- score_if_0 + abs(level_imbalance - 1L)
      }

      q <- rep.int(0.5, B)
      q[score_if_1 < score_if_0] <- P_BIASED_COIN
      q[score_if_1 > score_if_0] <- 1 - P_BIASED_COIN
      z <- ifelse(runif(B) < q, 1L, -1L)

      id <- as.integer(1L + drop(x %*% powers))
      idx <- cbind(row_index, id)
      stratum_counts[idx] <- stratum_counts[idx] + 1L
      stratum_imbalance[idx] <- stratum_imbalance[idx] + z

      d_overall <- d_overall + z
      d_level1 <- d_level1 + x * z
    }

    U <- matrix(0, nrow = B, ncol = J)
    nonempty <- stratum_counts > 0L
    U[nonempty] <- stratum_imbalance[nonempty] /
      sqrt(stratum_counts[nonempty])

    sum_u <- sum_u + colSums(U)
    sum_uu <- sum_uu + crossprod(U)
    n_treated <- as.integer((n + d_overall) / 2L)
    treated_count <- treated_count + tabulate(n_treated + 1L, nbins = n + 1L)

    completed <- completed + B
    message(
      "  calibration: ", format(completed, big.mark = ","), "/",
      format(B0, big.mark = ","), " allocation paths"
    )
  }

  mean_u <- sum_u / B0
  gamma <- (sum_uu - B0 * tcrossprod(mean_u)) / (B0 - 1L)
  gamma <- nearest_psd(gamma)

  structure(
    list(
      gamma = gamma,
      mean_u = mean_u,
      treated_count_prob = treated_count / B0,
      B0 = B0,
      n = n,
      K = K,
      J = J,
      probs = probs,
      pbc = P_BIASED_COIN,
      seed = calibration_seed(n, K),
      elapsed_seconds = proc.time()[3L] - start
    ),
    class = c("siga_calibration", "list")
  )
}

get_design_calibration <- function(scenario) {
  cache <- calibration_cache_file(scenario)
  if (file.exists(cache)) return(readRDS(cache))

  lock <- paste0(cache, ".lock")
  have_lock <- dir.create(lock, showWarnings = FALSE)
  if (!have_lock) {
    for (attempt in seq_len(360L)) {
      if (file.exists(cache)) return(readRDS(cache))
      Sys.sleep(5)
      have_lock <- dir.create(lock, showWarnings = FALSE)
      if (have_lock) break
    }
  }
  if (!have_lock) stop("Could not acquire calibration lock: ", lock)
  on.exit(unlink(lock, recursive = TRUE, force = TRUE), add = TRUE)

  if (file.exists(cache)) return(readRDS(cache))
  message(
    "Creating allocation-only calibration for K=", scenario$factor_count,
    ", n=", scenario$total_n, ", B0=", N_CALIBRATION
  )
  calibration <- run_allocation_calibration(scenario)
  tmp <- tempfile(pattern = "calibration_", tmpdir = CALIBRATION_DIR, fileext = ".rds")
  saveRDS(calibration, tmp)
  if (!file.rename(tmp, cache)) {
    if (!file.copy(tmp, cache, overwrite = TRUE)) {
      unlink(tmp)
      stop("Could not save calibration: ", cache)
    }
    unlink(tmp)
  }
  calibration
}

# =============================================================================
# Boundary outcomes and fixed scores for the proposed approximation
# =============================================================================

construct_boundary_data <- function(trial, boundaries) {
  n <- length(trial$Y)
  boundary_count <- length(boundaries)
  W <- matrix(NA_real_, nrow = n, ncol = boundary_count)
  for (j in seq_along(boundaries)) {
    W[, j] <- trial$Y - boundaries[j] * trial$A
  }

  if (!ADJUSTED) {
    scores <- sweep(W, 2L, colMeans(W), "-")
  } else {
    design <- cbind(Intercept = 1, trial$X)
    qr_design <- qr(design)
    scores <- matrix(NA_real_, nrow = n, ncol = boundary_count)
    for (j in seq_len(boundary_count)) {
      scores[, j] <- qr.resid(qr_design, W[, j])
    }
  }

  if (is.null(dim(scores))) scores <- matrix(scores, ncol = 1L)
  colnames(W) <- colnames(scores) <- names(boundaries)
  if (any(!is.finite(scores)) || any(!is.finite(W))) {
    stop("Non-finite boundary data generated.")
  }
  list(W = W, scores = scores)
}

score_decomposition <- function(score, X, J) {
  id <- joint_stratum_id(X)
  counts <- tabulate(id, nbins = J)
  sums <- numeric(J)
  tmp <- rowsum(score, group = id, reorder = FALSE)
  sums[as.integer(rownames(tmp))] <- drop(tmp)
  means <- numeric(J)
  nonempty <- counts > 0L
  means[nonempty] <- sums[nonempty] / counts[nonempty]
  residual <- score - means[id]
  list(counts = counts, means = means, residual = residual)
}

siga_variance_one <- function(score, X, calibration) {
  n <- length(score)
  dec <- score_decomposition(score, X, calibration$J)
  sqrt_counts <- sqrt(dec$counts)
  Omega <- tcrossprod(sqrt_counts) * calibration$gamma
  nonempty <- dec$counts > 0L
  Jplus <- sum(nonempty)
  trace_low <- sum(diag(Omega)[nonempty] / dec$counts[nonempty])
  denominator <- n - Jplus
  kappa <- if (denominator > 0L) (n - trace_low) / denominator else 0
  kappa <- max(kappa, 0)

  low <- drop(crossprod(dec$means, Omega %*% dec$means))
  within <- sum(dec$residual ^ 2)
  variance <- 0.25 * (low + kappa * within)
  if (!is.finite(variance) || variance <= .Machine$double.eps) {
    stop("Degenerate SIGA variance.")
  }
  variance
}

pvalue_from_z <- function(z, objective) {
  if (objective == "superiority") {
    return(list(overall = min(1, 2 * pnorm(-abs(z[1L]))),
                lower = NA_real_, upper = NA_real_))
  }
  if (objective == "noninferiority") {
    p <- pnorm(z[1L], lower.tail = FALSE)
    return(list(overall = p, lower = p, upper = NA_real_))
  }
  p_lower <- pnorm(z[1L], lower.tail = FALSE)
  p_upper <- pnorm(z[2L])
  list(overall = max(p_lower, p_upper), lower = p_lower, upper = p_upper)
}

siga_analysis <- function(trial, scenario, calibration) {
  boundaries <- scenario_boundaries(scenario)
  boundary_data <- construct_boundary_data(trial, boundaries)
  scores <- boundary_data$scores
  statistic <- colSums(scores * (trial$A - 0.5))
  variance <- vapply(
    seq_len(ncol(scores)),
    function(j) siga_variance_one(scores[, j], trial$X, calibration),
    numeric(1L)
  )
  z <- statistic / sqrt(variance)
  p <- pvalue_from_z(z, scenario$objective)
  list(
    p = p$overall,
    p_lower = p$lower,
    p_upper = p$upper,
    statistic = statistic,
    variance = variance,
    z = z
  )
}

# =============================================================================
# Fixed-score conditional randomization test
# =============================================================================

rerandomized_fixed_score_statistics <- function(X, scores, B) {
  n <- nrow(X)
  K <- ncol(X)
  boundary_count <- ncol(scores)
  d_overall <- integer(B)
  d_level1 <- matrix(0L, nrow = B, ncol = K)
  statistic <- matrix(0, nrow = B, ncol = boundary_count)

  for (i in seq_len(n)) {
    score_if_1 <- abs(d_overall + 1L)
    score_if_0 <- abs(d_overall - 1L)
    for (j in seq_len(K)) {
      level_imbalance <- if (X[i, j] == 1L) {
        d_level1[, j]
      } else {
        d_overall - d_level1[, j]
      }
      score_if_1 <- score_if_1 + abs(level_imbalance + 1L)
      score_if_0 <- score_if_0 + abs(level_imbalance - 1L)
    }

    q_star <- rep.int(0.5, B)
    q_star[score_if_1 < score_if_0] <- P_BIASED_COIN
    q_star[score_if_1 > score_if_0] <- 1 - P_BIASED_COIN
    A_star <- as.integer(runif(B) < q_star)
    z_star <- 2L * A_star - 1L
    centered_assignment <- A_star - 0.5

    for (j in seq_len(boundary_count)) {
      statistic[, j] <- statistic[, j] + centered_assignment * scores[i, j]
    }

    d_overall <- d_overall + z_star
    active_levels <- which(X[i, ] == 1L)
    if (length(active_levels)) {
      for (j in active_levels) d_level1[, j] <- d_level1[, j] + z_star
    }
  }

  statistic
}

monte_carlo_pvalue <- function(statistic_star, statistic_observed, objective) {
  B <- nrow(statistic_star)
  if (objective == "superiority") {
    p <- (1 + sum(abs(statistic_star[, 1L]) >=
                   abs(statistic_observed[1L]))) / (B + 1)
    return(list(overall = p, lower = NA_real_, upper = NA_real_))
  }
  if (objective == "noninferiority") {
    p <- (1 + sum(statistic_star[, 1L] >= statistic_observed[1L])) / (B + 1)
    return(list(overall = p, lower = p, upper = NA_real_))
  }
  p_lower <- (1 + sum(statistic_star[, 1L] >= statistic_observed[1L])) / (B + 1)
  p_upper <- (1 + sum(statistic_star[, 2L] <= statistic_observed[2L])) / (B + 1)
  list(overall = max(p_lower, p_upper), lower = p_lower, upper = p_upper)
}

randomization_test_analysis <- function(trial, scenario) {
  boundaries <- scenario_boundaries(scenario)
  boundary_data <- construct_boundary_data(trial, boundaries)
  scores <- boundary_data$scores
  statistic_observed <- colSums(scores * (trial$A - 0.5))
  statistic_star <- rerandomized_fixed_score_statistics(
    trial$X, scores, N_RERANDOMIZATIONS
  )
  p <- monte_carlo_pvalue(
    statistic_star, statistic_observed, scenario$objective
  )
  list(
    p = p$overall,
    p_lower = p$lower,
    p_upper = p$upper,
    statistic_observed = statistic_observed
  )
}

# =============================================================================
# Sharded, checkpointed execution
# =============================================================================

scenario_paths <- function(scenario) {
  scenario_dir <- ensure_directory(file.path(
    CHUNK_DIR,
    sprintf("scenario_%02d_%s", scenario$scenario_id, scenario$scenario_label)
  ))
  list(
    dir = scenario_dir,
    detail = file.path(
      scenario_dir,
      sprintf("shard_%04d_of_%04d.csv", SHARD_ID, N_SHARDS)
    ),
    timing = file.path(
      scenario_dir,
      sprintf("shard_%04d_of_%04d_timing.csv", SHARD_ID, N_SHARDS)
    )
  )
}

read_completed_replicates <- function(detail_file, attempts = 5L, wait_seconds = 2) {
  if (!file.exists(detail_file)) return(integer())
  for (attempt in seq_len(attempts)) {
    x <- tryCatch(
      read.csv(detail_file, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (!is.null(x) && "replicate" %in% names(x)) {
      ids <- unique(as.integer(x$replicate))
      return(ids[is.finite(ids)])
    }
    if (attempt < attempts) Sys.sleep(wait_seconds)
  }
  stop("Could not read checkpoint file: ", detail_file, call. = FALSE)
}

read_timing_checkpoint <- function(timing_file) {
  if (!file.exists(timing_file)) return(NULL)
  tryCatch(read.csv(timing_file, stringsAsFactors = FALSE), error = function(e) NULL)
}

write_timing_checkpoint <- function(path, timing) {
  ensure_directory(dirname(path))
  tmp <- tempfile(pattern = "timing_", tmpdir = dirname(path), fileext = ".tmp")
  write.csv(timing, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) {
    if (!file.copy(tmp, path, overwrite = TRUE)) {
      unlink(tmp)
      stop("Could not write timing checkpoint: ", path, call. = FALSE)
    }
    unlink(tmp)
  }
  invisible(NULL)
}

run_replicate <- function(scenario, replicate_id, calibration) {
  set.seed(seed_value(scenario$scenario_id, replicate_id, 1L))
  generation_start <- proc.time()[3L]
  trial <- simulate_trial(scenario)
  generation_seconds <- proc.time()[3L] - generation_start

  proposed_start <- proc.time()[3L]
  proposed <- siga_analysis(trial, scenario, calibration)
  proposed_seconds <- proc.time()[3L] - proposed_start

  set.seed(seed_value(scenario$scenario_id, replicate_id, 2L))
  randomization_start <- proc.time()[3L]
  randomization <- randomization_test_analysis(trial, scenario)
  randomization_seconds <- proc.time()[3L] - randomization_start

  data.frame(
    scenario_id = scenario$scenario_id,
    scenario_label = scenario$scenario_label,
    replicate = replicate_id,
    factor_count = scenario$factor_count,
    target_per_group = scenario$target_per_group,
    total_n = scenario$total_n,
    objective = scenario$objective,
    scenario_role = scenario$scenario_role,
    true_effect = scenario$true_effect,
    alpha = scenario$alpha,
    p_proposed = proposed$p,
    p_proposed_lower = proposed$p_lower,
    p_proposed_upper = proposed$p_upper,
    p_randomization = randomization$p,
    p_randomization_lower = randomization$p_lower,
    p_randomization_upper = randomization$p_upper,
    reject_proposed = proposed$p <= scenario$alpha,
    reject_randomization = randomization$p <= scenario$alpha,
    data_generation_seconds = generation_seconds,
    proposed_analysis_seconds = proposed_seconds,
    randomization_test_seconds = randomization_seconds,
    stringsAsFactors = FALSE
  )
}

run_scenario_shard <- function(scenario) {
  paths <- scenario_paths(scenario)
  assigned <- seq.int(from = SHARD_ID, to = N_OUTER, by = N_SHARDS)
  completed <- read_completed_replicates(paths$detail)
  pending <- setdiff(assigned, completed)

  message(
    "Scenario ", scenario$scenario_id, " [", scenario$scenario_label, "]",
    "; shard ", SHARD_ID, "/", N_SHARDS,
    "; assigned=", length(assigned),
    "; completed=", length(completed),
    "; pending=", length(pending)
  )

  previous_timing <- read_timing_checkpoint(paths$timing)
  previous_wall <- 0
  if (!is.null(previous_timing) &&
      "shard_wall_clock_seconds" %in% names(previous_timing) &&
      is.finite(previous_timing$shard_wall_clock_seconds[1L])) {
    previous_wall <- previous_timing$shard_wall_clock_seconds[1L]
  }

  # A completed scenario is left untouched. This avoids overwriting its original
  # timing information when a full script is rerun only to resume later scenarios.
  if (!length(pending)) {
    message("  already complete; no simulation was rerun.")
    if (!is.null(previous_timing)) return(invisible(previous_timing))
  }

  session_start <- proc.time()[3L]
  if (length(pending)) {
    calibration <- get_design_calibration(scenario)
    batches <- split(pending, ceiling(seq_along(pending) / OUTER_BATCH_SIZE))
    for (batch_index in seq_along(batches)) {
      ids <- batches[[batch_index]]
      batch_start <- proc.time()[3L]

      worker <- function(id) run_replicate(scenario, id, calibration)
      pieces <- if (USE_PARALLEL) {
        parallel::mclapply(
          ids,
          worker,
          mc.cores = N_CORES,
          mc.preschedule = FALSE,
          mc.set.seed = FALSE
        )
      } else {
        lapply(ids, worker)
      }
      batch <- do.call(rbind, pieces)
      batch <- batch[order(batch$replicate), , drop = FALSE]

      # Recreate the scenario directory immediately before every checkpoint.
      # This protects a long run from a transient cloud-storage directory loss.
      paths <- scenario_paths(scenario)
      append_csv(batch, paths$detail)

      completed_now <- length(unique(c(completed, unlist(batches[seq_len(batch_index)]))))
      cumulative_wall <- previous_wall + (proc.time()[3L] - session_start)
      timing_checkpoint <- data.frame(
        method_id = METHOD_ID,
        scenario_id = scenario$scenario_id,
        scenario_label = scenario$scenario_label,
        shard_id = SHARD_ID,
        n_shards = N_SHARDS,
        assigned_replicates = length(assigned),
        completed_replicates = completed_now,
        shard_wall_clock_seconds = cumulative_wall,
        shard_completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
        stringsAsFactors = FALSE
      )
      write_timing_checkpoint(paths$timing, timing_checkpoint)

      message(
        "  batch ", batch_index, "/", length(batches),
        ": replicates ", min(ids), "--", max(ids),
        "; elapsed=", format_seconds(proc.time()[3L] - batch_start)
      )
      rm(pieces, batch)
      gc(verbose = FALSE)
    }
  }

  paths <- scenario_paths(scenario)
  if (!file.exists(paths$detail)) {
    stop("Checkpoint file is missing after scenario execution: ", paths$detail,
         call. = FALSE)
  }
  detail <- read.csv(paths$detail, stringsAsFactors = FALSE)
  detail <- detail[!duplicated(detail$replicate), , drop = FALSE]
  detail <- detail[detail$replicate %in% assigned, , drop = FALSE]

  cumulative_wall <- previous_wall + (proc.time()[3L] - session_start)
  timing <- data.frame(
    method_id = METHOD_ID,
    scenario_id = scenario$scenario_id,
    scenario_label = scenario$scenario_label,
    shard_id = SHARD_ID,
    n_shards = N_SHARDS,
    assigned_replicates = length(assigned),
    completed_replicates = nrow(detail),
    shard_wall_clock_seconds = cumulative_wall,
    shard_completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    stringsAsFactors = FALSE
  )
  write_timing_checkpoint(paths$timing, timing)
  invisible(timing)
}

# =============================================================================
# Aggregation and manuscript-ready tables
# =============================================================================

summarize_scenario <- function(scenario) {
  scenario_dir_pattern <- sprintf("scenario_%02d_", scenario$scenario_id)
  dirs <- list.dirs(CHUNK_DIR, recursive = FALSE, full.names = TRUE)
  scenario_dirs <- dirs[startsWith(basename(dirs), scenario_dir_pattern)]
  if (!length(scenario_dirs)) return(NULL)

  detail_files <- unlist(lapply(
    scenario_dirs,
    function(d) list.files(d, pattern = "^shard_[0-9]+_of_[0-9]+\\.csv$", full.names = TRUE)
  ))
  if (!length(detail_files)) return(NULL)

  pieces <- lapply(detail_files, function(f) {
    tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (!length(pieces)) return(NULL)
  detail <- do.call(rbind, pieces)
  detail <- detail[detail$scenario_id == scenario$scenario_id, , drop = FALSE]
  detail <- detail[!duplicated(detail$replicate), , drop = FALSE]
  detail <- detail[order(detail$replicate), , drop = FALSE]

  n <- nrow(detail)
  proposed_count <- sum(detail$reject_proposed)
  rt_count <- sum(detail$reject_randomization)
  proposed_ci <- wilson_interval(proposed_count, n)
  rt_ci <- wilson_interval(rt_count, n)

  timing_files <- unlist(lapply(
    scenario_dirs,
    function(d) list.files(d, pattern = "_timing\\.csv$", full.names = TRUE)
  ))
  timing_pieces <- lapply(timing_files, function(f) {
    tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  })
  timing_pieces <- Filter(Negate(is.null), timing_pieces)
  timing <- if (length(timing_pieces)) do.call(rbind, timing_pieces) else NULL

  calibration <- get_design_calibration(scenario)

  data.frame(
    method_id = METHOD_ID,
    adjusted = ADJUSTED,
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
    requested_outer_trials = N_OUTER,
    completed_outer_trials = n,
    rerandomizations_per_trial = N_RERANDOMIZATIONS,
    calibration_paths = calibration$B0,
    calibration_minutes = calibration$elapsed_seconds / 60,
    proposed_rejection = unname(proposed_ci["estimate"]),
    proposed_lower_95 = unname(proposed_ci["lower"]),
    proposed_upper_95 = unname(proposed_ci["upper"]),
    randomization_rejection = unname(rt_ci["estimate"]),
    randomization_lower_95 = unname(rt_ci["lower"]),
    randomization_upper_95 = unname(rt_ci["upper"]),
    pvalue_mean_absolute_error = mean(abs(detail$p_proposed - detail$p_randomization)),
    pvalue_root_mean_squared_error = sqrt(mean((detail$p_proposed - detail$p_randomization)^2)),
    pvalue_correlation = if (n >= 2L) cor(detail$p_proposed, detail$p_randomization) else NA_real_,
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

tex_escape <- function(x) {
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x
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

write_tex_results <- function(summary, path) {
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
    paste0("\\label{tab:", METHOD_ID, "-operating-characteristics}"),
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
      N_CALIBRATION, " paths and its time is reported separately in minutes in the CSV output.}"
    ),
    paste0("\\label{tab:", METHOD_ID, "-timing}"),
    "\\end{table}"
  )
  writeLines(lines, con = path, useBytes = TRUE)
}

write_settings <- function(path) {
  total_updates <- sum(SCENARIOS$total_n) * N_OUTER * N_RERANDOMIZATIONS
  lines <- c(
    paste0("Method: ", METHOD_TITLE),
    paste0("Adjusted: ", ADJUSTED),
    paste0("Reference statistic: fixed treatment-score statistic T=sum(A-1/2)r"),
    paste0("Outer trials per scenario: ", N_OUTER),
    paste0("Rerandomizations per outer trial: ", N_RERANDOMIZATIONS),
    paste0("Allocation-only calibration paths per design: ", N_CALIBRATION),
    paste0("Calibration cache directory: ", CALIBRATION_DIR),
    paste0("Number of selected scenarios: ", nrow(SCENARIOS)),
    paste0("Biased-coin probability: ", P_BIASED_COIN),
    paste0("Power target band: 80%--90% (nominal target approximately 85%)"),
    paste0("Shards: ", N_SHARDS),
    paste0("Current shard: ", SHARD_ID),
    paste0("Cores: ", N_CORES),
    paste0("Approximate requested participant-by-rerandomization updates: ",
           format(total_updates, scientific = TRUE)),
    paste0("Output directory: ", OUTPUT_DIR),
    paste0("Generated at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  )
  writeLines(lines, con = path, useBytes = TRUE)
}

aggregate_results <- function() {
  summaries <- lapply(seq_len(nrow(SCENARIOS)), function(i) summarize_scenario(SCENARIOS[i, ]))
  summaries <- Filter(Negate(is.null), summaries)
  if (!length(summaries)) {
    stop("No shard result files were found for aggregation.", call. = FALSE)
  }
  summary <- do.call(rbind, summaries)
  summary <- summary[order(summary$scenario_id), , drop = FALSE]
  write.csv(
    summary,
    file.path(OUTPUT_DIR, paste0(METHOD_ID, "_summary.csv")),
    row.names = FALSE
  )
  write_tex_results(
    summary,
    file.path(OUTPUT_DIR, paste0(METHOD_ID, "_results.tex"))
  )
  write_settings(file.path(OUTPUT_DIR, paste0(METHOD_ID, "_settings.txt")))
  summary
}

# =============================================================================
# Main
# =============================================================================

message("============================================================")
message(METHOD_TITLE)
message("Mode: ", RUN_MODE)
message("Project directory: ", PROJECT_DIR)
message("Output directory: ", OUTPUT_DIR)
message("Outer trials per scenario: ", format(N_OUTER, big.mark = ","))
message("Rerandomizations per trial: ", format(N_RERANDOMIZATIONS, big.mark = ","))
message("Allocation-only calibration paths: ", format(N_CALIBRATION, big.mark = ","))
message("Selected scenarios: ", nrow(SCENARIOS))
message("Shard: ", SHARD_ID, "/", N_SHARDS)
message("Parallel workers: ", if (USE_PARALLEL) N_CORES else 1L)
message("============================================================")

if (RUN_MODE == "run") {
  for (i in seq_len(nrow(SCENARIOS))) {
    run_scenario_shard(SCENARIOS[i, ])
  }
  if (N_SHARDS == 1L) {
    result <- aggregate_results()
    print(result[, c(
      "scenario_id", "scenario_label", "completed_outer_trials",
      "proposed_rejection", "randomization_rejection"
    )], row.names = FALSE)
  } else {
    message(
      "Shard completed. After all shards finish, run with ",
      "PWRT_MODE=aggregate and the same PWRT_N_SHARDS value."
    )
  }
} else {
  result <- aggregate_results()
  print(result[, c(
    "scenario_id", "scenario_label", "completed_outer_trials",
    "proposed_rejection", "randomization_rejection"
  )], row.names = FALSE)
}
