#!/usr/bin/env Rscript

# =============================================================================
# Add SIGA-R to the completed SWIFT DIRECT-inspired power illustration
# =============================================================================
#
# Scientific purpose
# ------------------
# The original case-study simulation already contains 100,000 outer trials and
# the corresponding reference randomization-test p-values. This script does
# NOT rerun the 4,999-path randomization test. Instead, it:
#   1. reads the completed case-study shard files;
#   2. deterministically regenerates the same outer trials from their original
#      replicate identifiers and seeds;
#   3. performs the pair-path allocation-only calibration once;
#   4. evaluates SIGA-R on the same outer trials; and
#   5. joins SIGA-R with the existing SIGA-S and randomization-test results.
#
# The script also verifies that the regenerated SIGA-S p-values reproduce the
# existing case-study p-values. Aggregation stops if this audit fails.
#
# Required companion file
# -----------------------
#   siga_pair_path_engine.R
# Place it beside this script or set SWIFT_PAIR_ENGINE.
#
# Modes
# -----
#   SWIFT_R_MODE=calibrate : create/validate the reusable pair-path calibration
#   SWIFT_R_MODE=run       : run/resume one SIGA-R shard
#   SWIFT_R_MODE=aggregate : aggregate shards and create manuscript outputs
#
# Default locations reproduce the original case-study program:
#   project root:
#     ~/Dropbox/Tex/03_Post_PhD/029_exact_minimization/02_program
#   existing results:
#     swift_direct_siga_case_study_pureR_output
#   existing SIGA-S calibration:
#     swift_direct_siga_calibration_cache_pureR
#
# This extension evaluates the existing equal-response power scenario only.
# A new weak-null simulation is not required to demonstrate the practical
# distinction between SIGA-S and SIGA-R in this illustration.
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

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
  "SWIFT_PAIR_ENGINE",
  unset = file.path(BASE_DIR, "siga_pair_path_engine.R")
))
if (!file.exists(PAIR_ENGINE)) {
  stop(
    "Required pair-path engine not found: ", PAIR_ENGINE,
    "\nPlace siga_pair_path_engine.R beside this script or set SWIFT_PAIR_ENGINE.",
    call. = FALSE
  )
}
source(PAIR_ENGINE, local = FALSE)

# -----------------------------------------------------------------------------
# Environment and safe file helpers
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

ensure_directory <- function(path, attempts = 5L, wait_seconds = 1) {
  path <- path.expand(path)
  for (attempt in seq_len(attempts)) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (dir.exists(path)) {
      return(normalizePath(path, winslash = "/", mustWork = TRUE))
    }
    if (attempt < attempts) Sys.sleep(wait_seconds)
  }
  stop("Could not create or access directory: ", path, call. = FALSE)
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
      stop("Could not write: ", path, call. = FALSE)
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

# -----------------------------------------------------------------------------
# Paths and simulation settings
# -----------------------------------------------------------------------------

PROJECT_DIR <- path.expand(Sys.getenv(
  "PWRT_PROJECT_DIR",
  unset = "~/Dropbox/Tex/03_Post_PhD/029_exact_minimization/02_program"
))
OLD_OUTPUT_DIR <- path.expand(Sys.getenv(
  "SWIFT_OLD_OUTPUT_DIR",
  unset = file.path(PROJECT_DIR, "swift_direct_siga_case_study_pureR_output")
))
OLD_CALIBRATION_DIR <- path.expand(Sys.getenv(
  "SWIFT_OLD_CALIBRATION_DIR",
  unset = file.path(PROJECT_DIR, "swift_direct_siga_calibration_cache_pureR")
))
OUTPUT_DIR <- path.expand(Sys.getenv(
  "SWIFT_R_OUTPUT_DIR",
  unset = file.path(PROJECT_DIR, "swift_direct_siga_r_case_study_output")
))

if (!dir.exists(OLD_OUTPUT_DIR)) {
  stop("Existing case-study output directory not found: ", OLD_OUTPUT_DIR, call. = FALSE)
}
if (!dir.exists(OLD_CALIBRATION_DIR)) {
  stop("Existing SIGA-S calibration directory not found: ", OLD_CALIBRATION_DIR, call. = FALSE)
}
OUTPUT_DIR <- ensure_directory(OUTPUT_DIR)
SHARD_DIR <- ensure_directory(file.path(OUTPUT_DIR, "siga_r_shards"))
CALIBRATION_DIR <- ensure_directory(file.path(OUTPUT_DIR, "pair_calibration"))
PUBLICATION_DIR <- ensure_directory(file.path(OUTPUT_DIR, "publication"))

RUN_MODE <- tolower(trimws(Sys.getenv("SWIFT_R_MODE", unset = "run")))
if (!RUN_MODE %in% c("calibrate", "run", "aggregate")) {
  stop("SWIFT_R_MODE must be calibrate, run, or aggregate.", call. = FALSE)
}

N_OUTER <- env_integer("PWRT_N_OUTER", 100000L)
N_PAIR_CALIBRATION <- env_integer("PWRT_N_PAIR_CALIBRATION", 100000L)
PAIR_CALIBRATION_BATCH <- env_integer("PWRT_CALIBRATION_BATCH", 1000L)
N_SHARDS <- env_integer("SWIFT_R_N_SHARDS", 1L)
SHARD_ID <- env_integer("SWIFT_R_SHARD_ID", 1L)
OUTER_BATCH <- env_integer("SWIFT_R_OUTER_BATCH", 25L)
P_BIASED_COIN <- env_numeric("PWRT_P_BIASED_COIN", 0.80)
BASE_SEED <- env_integer("PWRT_SEED", 20260724L)
EPSILON_EXPONENT <- env_numeric("PWRT_EPSILON_EXPONENT", 1.0)
REPRODUCTION_TOLERANCE <- env_numeric("SWIFT_R_REPRO_TOL", 1e-10)
ALLOW_PARTIAL <- identical(Sys.getenv("SWIFT_R_ALLOW_PARTIAL", unset = "0"), "1")

stopifnot(
  N_OUTER >= 1L,
  N_PAIR_CALIBRATION >= 2L,
  PAIR_CALIBRATION_BATCH >= 1L,
  N_SHARDS >= 1L,
  SHARD_ID >= 1L,
  SHARD_ID <= N_SHARDS,
  OUTER_BATCH >= 1L,
  P_BIASED_COIN > 0.5,
  P_BIASED_COIN < 1,
  EPSILON_EXPONENT > 0,
  REPRODUCTION_TOLERANCE > 0
)

# -----------------------------------------------------------------------------
# Exact settings copied from the completed SWIFT DIRECT-inspired simulation
# -----------------------------------------------------------------------------

SCENARIO_ID <- 1L
TOTAL_N <- 404L
FACTOR_COUNT <- 5L
NONINFERIORITY_BOUNDARY <- -0.12
ALPHA <- 0.05
TARGET_CONTROL_RISK <- 0.622
BINARY_BETA <- c(-0.55, -0.25, -0.40, -0.25, -0.60)
BINARY_INTERACTION <- 0.00

AGE_MEDIAN <- 72.5
AGE_Q1 <- 64.5
AGE_Q3 <- 81.0
AGE_SD <- (AGE_Q3 - AGE_Q1) / (qnorm(0.75) - qnorm(0.25))
AGE_LOWER <- 18
AGE_UPPER <- 100
NIHSS_MEDIAN <- 17.0
NIHSS_Q1 <- (201 * 13 + 207 * 12) / 408
NIHSS_Q3 <- 20.0
NIHSS_SD <- (NIHSS_Q3 - NIHSS_Q1) / (qnorm(0.75) - qnorm(0.25))
NIHSS_LOWER <- 5
NIHSS_UPPER <- 30

truncated_tail_probability <- function(cut, mean, sd, lower, upper) {
  lo <- pnorm(lower, mean = mean, sd = sd)
  hi <- pnorm(upper, mean = mean, sd = sd)
  cc <- pnorm(cut, mean = mean, sd = sd)
  (hi - cc) / (hi - lo)
}

P_NIHSS_GT17 <- truncated_tail_probability(
  17.5, NIHSS_MEDIAN, NIHSS_SD, NIHSS_LOWER, NIHSS_UPPER
)
P_AGE_GE70 <- truncated_tail_probability(
  69.5, AGE_MEDIAN, AGE_SD, AGE_LOWER, AGE_UPPER
)
P_ICA_RELATED <- 117 / 408
P_TANDEM <- 63 / 408
ASPECTS_VALUES <- 4:10
ASPECTS_PROB <- c(0.01, 0.03, 0.08, 0.18, 0.30, 0.25, 0.15)
P_ASPECTS_4_7 <- sum(ASPECTS_PROB[ASPECTS_VALUES <= 7])
FACTOR_PROB <- c(
  P_NIHSS_GT17,
  P_AGE_GE70,
  P_ICA_RELATED,
  P_TANDEM,
  P_ASPECTS_4_7
)
names(FACTOR_PROB) <- c(
  "NIHSS >17",
  "Age >=70",
  "ICA-related occlusion",
  "Tandem lesion",
  "ASPECTS 4--7"
)

PATTERNS <- all_binary_patterns(FACTOR_COUNT)
PROFILE_PROB <- profile_prob_independent(FACTOR_PROB, PATTERNS)
WEIGHTS <- rep(1, FACTOR_COUNT + 1L)

# Reproduce the original logistic outcome model exactly. Under the power
# scenario, theta is set to zero, so p1(x) = p0(x) for every profile.
calibrate_original_outcome_model <- function() {
  nonintercept <- drop(PATTERNS %*% BINARY_BETA)
  if (FACTOR_COUNT >= 2L) {
    nonintercept <- nonintercept +
      BINARY_INTERACTION * PATTERNS[, 1L] * PATTERNS[, 2L]
  }
  mean_control <- function(alpha0) {
    sum(PROFILE_PROB * plogis(alpha0 + nonintercept))
  }
  intercept <- uniroot(
    function(a) mean_control(a) - TARGET_CONTROL_RISK,
    interval = c(-30, 30),
    tol = 1e-12
  )$root
  p0 <- plogis(intercept + nonintercept)
  list(
    intercept = intercept,
    beta = BINARY_BETA,
    interaction = BINARY_INTERACTION,
    theta = 0,
    p0 = p0,
    p1 = p0,
    achieved_control_risk = sum(PROFILE_PROB * p0),
    achieved_treatment_risk = sum(PROFILE_PROB * p0),
    achieved_risk_difference = 0
  )
}
OUTCOME_MODEL <- calibrate_original_outcome_model()
D_BOUNDARY <- OUTCOME_MODEL$p1 - OUTCOME_MODEL$p0 - NONINFERIORITY_BOUNDARY
stopifnot(
  length(D_BOUNDARY) == nrow(PATTERNS),
  max(abs(D_BOUNDARY - 0.12)) < 1e-12
)

# -----------------------------------------------------------------------------
# Reproducible regeneration of the completed outer trials
# -----------------------------------------------------------------------------

seed_value <- function(scenario_id, replicate_id, stream_id) {
  modulus <- 2147483646
  value <- as.double(BASE_SEED) +
    1000003 * as.double(stream_id) +
    104729 * as.double(scenario_id) +
    1009 * as.double(replicate_id)
  as.integer(value %% modulus + 1)
}

pair_calibration_seed <- function() {
  modulus <- 2147483646
  value <- as.double(BASE_SEED) + 17000003 +
    1009 * TOTAL_N + 104729 * FACTOR_COUNT
  as.integer(value %% modulus + 1)
}

generate_factors_original <- function(n) {
  X <- matrix(0L, nrow = n, ncol = FACTOR_COUNT)
  for (j in seq_len(FACTOR_COUNT)) {
    X[, j] <- rbinom(n, size = 1L, prob = FACTOR_PROB[j])
  }
  colnames(X) <- paste0("X", seq_len(FACTOR_COUNT))
  X
}

generate_outcome_original <- function(X, A) {
  eta0 <- OUTCOME_MODEL$intercept + drop(X %*% OUTCOME_MODEL$beta)
  if (FACTOR_COUNT >= 2L) {
    eta0 <- eta0 + OUTCOME_MODEL$interaction * X[, 1L] * X[, 2L]
  }
  p0 <- plogis(eta0)
  # theta is exactly zero in the completed power scenario.
  prob <- p0
  rbinom(length(A), size = 1L, prob = prob)
}

regenerate_trial <- function(replicate_id) {
  set.seed(seed_value(SCENARIO_ID, replicate_id, 1L))
  X <- generate_factors_original(TOTAL_N)
  z <- ps_assign_R(
    X = X,
    pbc = P_BIASED_COIN,
    weights_ = WEIGHTS,
    seed = seed_value(SCENARIO_ID, replicate_id, 2L)
  )
  A <- as.integer((z + 1L) / 2L)
  set.seed(seed_value(SCENARIO_ID, replicate_id, 3L))
  y <- generate_outcome_original(X, A)
  list(X = X, z = as.integer(z), A = A, y = y)
}

make_case_study_scores <- function(y, A, X) {
  W <- as.numeric(y) - NONINFERIORITY_BOUNDARY * as.numeric(A)
  unadjusted <- W - mean(W)
  C <- cbind(`(Intercept)` = 1, X)
  adjusted <- qr.resid(qr(C), W)
  adjusted <- adjusted - mean(adjusted)
  list(unadjusted = unadjusted, adjusted = adjusted)
}

# -----------------------------------------------------------------------------
# Existing results and calibration
# -----------------------------------------------------------------------------

existing_result_columns <- function() {
  c(
    "replicate",
    "p_proposed_unadjusted", "p_proposed_adjusted",
    "p_reference_unadjusted", "p_reference_adjusted",
    "reject_proposed_unadjusted", "reject_proposed_adjusted",
    "reject_reference_unadjusted", "reject_reference_adjusted"
  )
}

read_existing_results_from_directory <- function(directory) {
  files <- list.files(
    directory,
    pattern = "^shard_[0-9]+_of_[0-9]+[.]csv$",
    full.names = TRUE
  )
  if (!length(files)) {
    stop("No existing case-study shard CSV files found in ", directory, ".", call. = FALSE)
  }
  dat <- do.call(rbind, lapply(files, function(path) {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }))
  required <- existing_result_columns()
  missing <- setdiff(required, names(dat))
  if (length(missing)) {
    stop(
      "Existing shards in ", directory, " are missing columns: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  dat$replicate <- as.integer(dat$replicate)
  dat <- dat[!duplicated(dat$replicate), , drop = FALSE]
  dat <- dat[order(dat$replicate), , drop = FALSE]
  rownames(dat) <- NULL
  dat
}

has_complete_existing_coverage <- function(dat) {
  nrow(dat) >= N_OUTER &&
    all(seq_len(N_OUTER) %in% as.integer(dat$replicate))
}

same_existing_results <- function(x, y) {
  required <- existing_result_columns()
  common <- intersect(as.integer(x$replicate), as.integer(y$replicate))
  if (!length(common)) return(FALSE)
  x2 <- x[match(common, x$replicate), required, drop = FALSE]
  y2 <- y[match(common, y$replicate), required, drop = FALSE]
  isTRUE(all.equal(x2, y2, tolerance = 0, check.attributes = FALSE))
}

format_existing_candidate <- function(directory, dat = NULL, error = NULL) {
  if (!is.null(error)) {
    return(paste0("  - ", directory, " [unreadable: ", error, "]"))
  }
  range_text <- if (nrow(dat)) {
    paste0(min(dat$replicate), "--", max(dat$replicate))
  } else {
    "empty"
  }
  paste0(
    "  - ", directory,
    " [", nrow(dat), " unique replicates; range ", range_text,
    if (has_complete_existing_coverage(dat)) "; covers requested replicates" else "; incomplete",
    "]"
  )
}

read_existing_results <- function() {
  override <- trimws(Sys.getenv("SWIFT_OLD_SCENARIO_DIR", unset = ""))
  if (nzchar(override)) {
    directory <- path.expand(override)
    if (!dir.exists(directory)) {
      stop("SWIFT_OLD_SCENARIO_DIR does not exist: ", directory, call. = FALSE)
    }
    dat <- read_existing_results_from_directory(directory)
    if (!ALLOW_PARTIAL && !has_complete_existing_coverage(dat)) {
      stop(
        "The directory selected by SWIFT_OLD_SCENARIO_DIR does not cover replicates 1--",
        N_OUTER, ". Found ", nrow(dat), " unique replicates.", call. = FALSE
      )
    }
    message("Using existing scenario directory specified by SWIFT_OLD_SCENARIO_DIR: ", directory)
    return(dat)
  }

  root <- file.path(OLD_OUTPUT_DIR, "scenario_shards")
  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl("^scenario_01_", basename(dirs))]
  if (!length(dirs)) {
    stop("No existing scenario_01 directory found under ", root, ".", call. = FALSE)
  }

  candidates <- lapply(dirs, function(directory) {
    tryCatch(
      list(directory = directory,
           data = read_existing_results_from_directory(directory),
           error = NULL),
      error = function(e) list(directory = directory, data = NULL,
                               error = conditionMessage(e))
    )
  })
  readable <- which(vapply(candidates, function(x) is.null(x$error), logical(1)))
  if (!length(readable)) {
    details <- vapply(candidates, function(x) {
      format_existing_candidate(x$directory, error = x$error)
    }, character(1))
    stop(
      "None of the scenario_01 directories could be read:\n",
      paste(details, collapse = "\n"), call. = FALSE
    )
  }

  complete <- readable[vapply(candidates[readable], function(x) {
    has_complete_existing_coverage(x$data)
  }, logical(1))]

  if (length(complete) == 1L) {
    chosen <- candidates[[complete]]
    message("Using the only existing scenario directory covering the requested replicates: ", chosen$directory)
    return(chosen$data)
  }

  if (length(complete) > 1L) {
    reference <- candidates[[complete[1L]]]$data
    identical_complete <- if (length(complete) > 1L) {
      vapply(complete[-1L], function(index) {
        same_existing_results(reference, candidates[[index]]$data)
      }, logical(1))
    } else logical()
    if (all(identical_complete)) {
      complete_dirs <- vapply(candidates[complete], `[[`, character(1), "directory")
      mtimes <- file.info(complete_dirs)$mtime
      chosen_index <- complete[which.max(mtimes)]
      chosen <- candidates[[chosen_index]]
      message(
        "Found ", length(complete),
        " suitable scenario_01 directories with identical required results; using: ",
        chosen$directory
      )
      return(chosen$data)
    }

    details <- vapply(candidates, function(x) {
      format_existing_candidate(x$directory, x$data, x$error)
    }, character(1))
    stop(
      "Multiple suitable scenario_01 directories contain different results. ",
      "Set SWIFT_OLD_SCENARIO_DIR explicitly to the intended directory.\n",
      paste(details, collapse = "\n"), call. = FALSE
    )
  }

  if (ALLOW_PARTIAL) {
    counts <- vapply(candidates[readable], function(x) nrow(x$data), integer(1))
    chosen_index <- readable[which.max(counts)]
    chosen <- candidates[[chosen_index]]
    message(
      "No scenario_01 directory covers all requested replicates; SWIFT_R_ALLOW_PARTIAL=1, so using ",
      "the readable candidate with the most replicates: ", chosen$directory
    )
    return(chosen$data)
  }

  details <- vapply(candidates, function(x) {
    format_existing_candidate(x$directory, x$data, x$error)
  }, character(1))
  stop(
    "No scenario_01 directory covers replicates 1--", N_OUTER, ".\n",
    paste(details, collapse = "\n"),
    "\nSet SWIFT_OLD_SCENARIO_DIR explicitly after confirming the intended directory, ",
    "or set SWIFT_R_ALLOW_PARTIAL=1 for a deliberate partial run.",
    call. = FALSE
  )
}

old_calibration_file <- function() {
  prob_tag <- paste(formatC(FACTOR_PROB, format = "f", digits = 3), collapse = "-")
  expected <- file.path(
    OLD_CALIBRATION_DIR,
    sprintf(
      "swift_direct_siga_pureR_K%d_n%d_p%s_pbc%.3f_B%d.rds",
      FACTOR_COUNT, TOTAL_N, prob_tag, P_BIASED_COIN,
      env_integer("PWRT_N_CALIBRATION", 100000L)
    )
  )
  if (file.exists(expected)) return(expected)
  candidates <- list.files(
    OLD_CALIBRATION_DIR,
    pattern = sprintf("^swift_direct_siga_pureR_K%d_n%d_.*[.]rds$", FACTOR_COUNT, TOTAL_N),
    full.names = TRUE
  )
  if (length(candidates) != 1L) {
    stop(
      "Could not uniquely identify the existing SIGA-S calibration. Expected: ",
      expected, "; candidates found: ", length(candidates), call. = FALSE
    )
  }
  warning("Using the only matching existing SIGA-S calibration: ", candidates[1L])
  candidates[1L]
}

load_old_calibration <- function() {
  cal <- readRDS(old_calibration_file())
  required <- c("gamma", "n", "K", "J")
  missing <- required[!vapply(required, function(x) !is.null(cal[[x]]), logical(1L))]
  if (length(missing)) {
    stop("Existing SIGA-S calibration lacks: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (cal$n != TOTAL_N || cal$K != FACTOR_COUNT || cal$J != 2L ^ FACTOR_COUNT) {
    stop("Existing SIGA-S calibration dimensions do not match the case study.", call. = FALSE)
  }
  cal$gamma <- nearest_psd(cal$gamma)
  cal
}

pair_calibration_file <- function() {
  file.path(
    CALIBRATION_DIR,
    sprintf(
      "swift_direct_pair_cal_K%d_n%d_pbc%.3f_B%d_seed%d.rds",
      FACTOR_COUNT, TOTAL_N, P_BIASED_COIN,
      N_PAIR_CALIBRATION, pair_calibration_seed()
    )
  )
}

get_pair_calibration <- function(create = TRUE) {
  path <- pair_calibration_file()
  if (file.exists(path)) return(readRDS(path))
  if (!isTRUE(create)) stop("Pair calibration has not been created: ", path, call. = FALSE)

  message(
    "Creating pair-path calibration: K=", FACTOR_COUNT,
    ", n=", TOTAL_N,
    ", B=", N_PAIR_CALIBRATION
  )
  cal <- calibrate_pair_path_design_R(
    B0 = N_PAIR_CALIBRATION,
    n = TOTAL_N,
    patterns = PATTERNS,
    profile_prob = PROFILE_PROB,
    pbc = P_BIASED_COIN,
    weights_ = WEIGHTS,
    seed = pair_calibration_seed(),
    batch_size = PAIR_CALIBRATION_BATCH,
    progress = TRUE
  )
  tmp <- tempfile(pattern = "swift_pair_cal_", tmpdir = CALIBRATION_DIR, fileext = ".rds")
  saveRDS(cal, tmp)
  if (!file.rename(tmp, path)) {
    if (!file.copy(tmp, path, overwrite = TRUE)) {
      unlink(tmp)
      stop("Could not save pair calibration: ", path, call. = FALSE)
    }
    unlink(tmp)
  }
  cal
}

make_hybrid_calibration <- function() {
  old <- load_old_calibration()
  pair <- get_pair_calibration(create = FALSE)
  pair$gamma <- old$gamma
  pair$mean_u <- old$mean_u %||% rep(0, pair$J)
  pair$n <- TOTAL_N
  pair$K <- FACTOR_COUNT
  pair$J <- 2L ^ FACTOR_COUNT
  pair$patterns <- PATTERNS
  pair$profile_prob <- PROFILE_PROB
  pair$pbc <- P_BIASED_COIN
  pair$weights <- WEIGHTS
  class(pair) <- c("siga_pair_calibration", "list")
  pair
}

# -----------------------------------------------------------------------------
# SIGA-R evaluation for one regenerated outer trial
# -----------------------------------------------------------------------------

randomization_variance_with_safeguard <- function(sampling_result, d, calibration) {
  d <- as.numeric(d)
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

evaluate_score <- function(score, trial, calibration) {
  statistic <- 0.5 * sum(trial$z * score)
  sampling <- siga_sampling_variance(score, trial$X, calibration)
  randomization <- randomization_variance_with_safeguard(
    sampling_result = sampling,
    d = D_BOUNDARY,
    calibration = calibration
  )
  z_s <- statistic / sqrt(sampling$variance)
  z_r <- statistic / sqrt(randomization$variance)
  list(
    statistic = statistic,
    variance_s = sampling$variance,
    variance_r_raw = randomization$variance_raw,
    variance_r = randomization$variance,
    correction = randomization$correction,
    ratio_raw = randomization$variance_raw / sampling$variance,
    safeguarded = randomization$safeguarded,
    p_s = pnorm(z_s, lower.tail = FALSE),
    p_r = pnorm(z_r, lower.tail = FALSE),
    z_s = z_s,
    z_r = z_r
  )
}

run_extension_replicate <- function(replicate_id, existing_row, calibration) {
  generation_start <- proc.time()[3L]
  trial <- regenerate_trial(replicate_id)
  generation_seconds <- proc.time()[3L] - generation_start

  score_start <- proc.time()[3L]
  scores <- make_case_study_scores(trial$y, trial$A, trial$X)
  score_seconds <- proc.time()[3L] - score_start

  analysis_start <- proc.time()[3L]
  unadjusted <- evaluate_score(scores$unadjusted, trial, calibration)
  adjusted <- evaluate_score(scores$adjusted, trial, calibration)
  analysis_seconds <- proc.time()[3L] - analysis_start

  data.frame(
    replicate = replicate_id,
    p_s_existing_unadjusted = existing_row$p_proposed_unadjusted,
    p_s_regenerated_unadjusted = unadjusted$p_s,
    p_r_unadjusted = unadjusted$p_r,
    p_rt_unadjusted = existing_row$p_reference_unadjusted,
    reject_s_existing_unadjusted = as.logical(existing_row$reject_proposed_unadjusted),
    reject_s_regenerated_unadjusted = unadjusted$p_s <= ALPHA,
    reject_r_unadjusted = unadjusted$p_r <= ALPHA,
    reject_rt_unadjusted = as.logical(existing_row$reject_reference_unadjusted),
    variance_s_unadjusted = unadjusted$variance_s,
    variance_r_raw_unadjusted = unadjusted$variance_r_raw,
    variance_r_unadjusted = unadjusted$variance_r,
    ratio_r_over_s_unadjusted = unadjusted$ratio_raw,
    safeguard_unadjusted = unadjusted$safeguarded,
    p_s_existing_adjusted = existing_row$p_proposed_adjusted,
    p_s_regenerated_adjusted = adjusted$p_s,
    p_r_adjusted = adjusted$p_r,
    p_rt_adjusted = existing_row$p_reference_adjusted,
    reject_s_existing_adjusted = as.logical(existing_row$reject_proposed_adjusted),
    reject_s_regenerated_adjusted = adjusted$p_s <= ALPHA,
    reject_r_adjusted = adjusted$p_r <= ALPHA,
    reject_rt_adjusted = as.logical(existing_row$reject_reference_adjusted),
    variance_s_adjusted = adjusted$variance_s,
    variance_r_raw_adjusted = adjusted$variance_r_raw,
    variance_r_adjusted = adjusted$variance_r,
    ratio_r_over_s_adjusted = adjusted$ratio_raw,
    safeguard_adjusted = adjusted$safeguarded,
    data_regeneration_seconds = generation_seconds,
    score_construction_seconds = score_seconds,
    siga_r_analysis_seconds = analysis_seconds,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Sharded execution
# -----------------------------------------------------------------------------

extension_shard_file <- function(shard_id = SHARD_ID, n_shards = N_SHARDS) {
  file.path(
    SHARD_DIR,
    sprintf("siga_r_shard_%04d_of_%04d.csv", shard_id, n_shards)
  )
}

read_completed_extension <- function(path) {
  if (!file.exists(path)) return(integer())
  dat <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(dat) || !"replicate" %in% names(dat)) {
    stop("Could not read extension checkpoint: ", path, call. = FALSE)
  }
  unique(as.integer(dat$replicate))
}

run_extension_shard <- function() {
  existing <- read_existing_results()
  calibration <- make_hybrid_calibration()
  assigned <- seq.int(from = SHARD_ID, to = N_OUTER, by = N_SHARDS)
  assigned <- intersect(assigned, as.integer(existing$replicate))
  path <- extension_shard_file()
  completed <- read_completed_extension(path)
  pending <- setdiff(assigned, completed)

  message(
    "SIGA-R shard ", SHARD_ID, "/", N_SHARDS,
    "; assigned=", length(assigned),
    "; completed=", length(intersect(assigned, completed)),
    "; pending=", length(pending)
  )
  if (!length(pending)) return(invisible(NULL))

  existing_index <- match(pending, existing$replicate)
  batches <- split(pending, ceiling(seq_along(pending) / OUTER_BATCH))
  session_start <- proc.time()[3L]

  for (batch_index in seq_along(batches)) {
    ids <- batches[[batch_index]]
    rows <- match(ids, existing$replicate)
    batch <- do.call(rbind, lapply(seq_along(ids), function(k) {
      run_extension_replicate(
        replicate_id = ids[k],
        existing_row = existing[rows[k], , drop = FALSE],
        calibration = calibration
      )
    }))
    append_csv(batch, path)
    if (batch_index %% 10L == 0L || batch_index == length(batches)) {
      message(
        "  shard ", SHARD_ID,
        ": batch ", batch_index, "/", length(batches),
        "; elapsed ", format_elapsed(proc.time()[3L] - session_start)
      )
    }
  }
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# Aggregation and manuscript outputs
# -----------------------------------------------------------------------------

wilson_interval_local <- function(x, n, conf.level = 0.95) {
  z <- qnorm(1 - (1 - conf.level) / 2)
  p <- x / n
  den <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(estimate = p, lower = max(0, center - half), upper = min(1, center + half))
}

read_all_extension_shards <- function() {
  files <- list.files(
    SHARD_DIR,
    pattern = "^siga_r_shard_[0-9]+_of_[0-9]+[.]csv$",
    full.names = TRUE
  )
  if (!length(files)) stop("No SIGA-R extension shard files found.", call. = FALSE)
  dat <- do.call(rbind, lapply(files, function(path) {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }))
  dat <- dat[!duplicated(dat$replicate), , drop = FALSE]
  dat <- dat[order(dat$replicate), , drop = FALSE]
  if (!ALLOW_PARTIAL) {
    if (nrow(dat) != N_OUTER || !identical(as.integer(dat$replicate), seq_len(N_OUTER))) {
      stop(
        "SIGA-R extension contains ", nrow(dat),
        " unique replicates; expected exactly 1--", N_OUTER, ".",
        call. = FALSE
      )
    }
  }
  dat
}

summarize_analysis <- function(dat, analysis) {
  suffix <- if (analysis == "Unadjusted") "unadjusted" else "adjusted"
  reject_s <- as.logical(dat[[paste0("reject_s_existing_", suffix)]])
  reject_r <- as.logical(dat[[paste0("reject_r_", suffix)]])
  reject_rt <- as.logical(dat[[paste0("reject_rt_", suffix)]])
  p_s_old <- as.numeric(dat[[paste0("p_s_existing_", suffix)]])
  p_s_new <- as.numeric(dat[[paste0("p_s_regenerated_", suffix)]])
  ratio <- as.numeric(dat[[paste0("ratio_r_over_s_", suffix)]])
  safeguard <- as.logical(dat[[paste0("safeguard_", suffix)]])
  n <- nrow(dat)
  ci_s <- wilson_interval_local(sum(reject_s), n)
  ci_r <- wilson_interval_local(sum(reject_r), n)
  ci_rt <- wilson_interval_local(sum(reject_rt), n)
  paired <- as.numeric(reject_r) - as.numeric(reject_rt)
  paired_se <- if (n > 1L) sd(paired) / sqrt(n) else NA_real_

  data.frame(
    analysis = analysis,
    n_outer = n,
    siga_s_power = unname(ci_s["estimate"]),
    siga_s_lower_95 = unname(ci_s["lower"]),
    siga_s_upper_95 = unname(ci_s["upper"]),
    siga_r_power = unname(ci_r["estimate"]),
    siga_r_lower_95 = unname(ci_r["lower"]),
    siga_r_upper_95 = unname(ci_r["upper"]),
    rt_power = unname(ci_rt["estimate"]),
    rt_lower_95 = unname(ci_rt["lower"]),
    rt_upper_95 = unname(ci_rt["upper"]),
    siga_r_minus_rt = mean(paired),
    siga_r_minus_rt_se = paired_se,
    siga_r_minus_rt_lower_95 = mean(paired) - qnorm(0.975) * paired_se,
    siga_r_minus_rt_upper_95 = mean(paired) + qnorm(0.975) * paired_se,
    mean_variance_ratio = mean(ratio),
    q025_variance_ratio = unname(quantile(ratio, 0.025, type = 8)),
    q975_variance_ratio = unname(quantile(ratio, 0.975, type = 8)),
    safeguard_rate = mean(safeguard),
    max_abs_siga_s_pvalue_reproduction_error = max(abs(p_s_new - p_s_old)),
    stringsAsFactors = FALSE
  )
}

write_results_table <- function(summary, path) {
  fmt <- function(x) sprintf("%.2f", 100 * x)
  rows <- vapply(seq_len(nrow(summary)), function(i) {
    r <- summary[i, ]
    paste0(
      r$analysis, " & ",
      fmt(r$siga_s_power), " & ",
      fmt(r$siga_r_power), " & ",
      fmt(r$rt_power), " & ",
      sprintf("%+.2f", 100 * r$siga_r_minus_rt), " & ",
      sprintf("%.3f", r$mean_variance_ratio), " \\\\"
    )
  }, character(1L))
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Estimated non-inferiority power in the SWIFT DIRECT-inspired illustration.}",
    "\\label{tab:swift-direct-siga-summary}",
    "\\begin{threeparttable}",
    "\\small",
    "\\setlength{\\tabcolsep}{5.5pt}",
    "\\renewcommand{\\arraystretch}{1.10}",
    "\\begin{tabular}{lrrrrr}",
    "\\toprule",
    "Score & SIGA-S (\\%) & SIGA-R (\\%) & RT (\\%) & R$-$RT (pp) & Mean $\\widehat\\rho_n$ \\\\ ",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{tablenotes}[flushleft]",
    "\\footnotesize",
    "\\item RT denotes the reference fixed-score randomization test. R$-$RT is the paired difference between SIGA-R and RT rejection probabilities in percentage points. The diagnostic $\\widehat\\rho_n$ is the raw SIGA-R variance divided by the SIGA-S variance.",
    "\\item All three procedures were evaluated on the same 100,000 outer trials. SIGA-S targets repeated-sampling power for the marginal risk difference, whereas SIGA-R targets the rejection probability of the prespecified reference randomization test.",
    "\\end{tablenotes}",
    "\\end{threeparttable}",
    "\\end{table}",
    ""
  )
  writeLines(lines, path, useBytes = TRUE)
}

write_timing_table <- function(dat, pair_calibration, path) {
  old_summary_path <- file.path(OLD_OUTPUT_DIR, "swift_direct_siga_case_study_summary.csv")
  old <- if (file.exists(old_summary_path)) {
    read.csv(old_summary_path, stringsAsFactors = FALSE)
  } else {
    NULL
  }
  siga_r_analysis <- sum(as.numeric(dat$siga_r_analysis_seconds)) / 60
  pair_minutes <- pair_calibration$elapsed_seconds / 60
  if (is.null(old) || nrow(old) != 1L) {
    lines <- c(
      "Pair calibration and SIGA-R analysis timing (minutes):",
      sprintf("pair calibration = %.4f", pair_minutes),
      sprintf("SIGA-R analysis = %.4f", siga_r_analysis),
      sprintf("SIGA-R total = %.4f", pair_minutes + siga_r_analysis)
    )
    writeLines(lines, sub("[.]tex$", ".txt", path), useBytes = TRUE)
    return(invisible(NULL))
  }
  s_cal <- old$calibration_minutes
  s_analysis <- old$proposed_analysis_minutes
  rt_analysis <- old$reference_analysis_minutes
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Computation time in the SWIFT DIRECT-inspired illustration.}",
    "\\label{tab:swift-direct-siga-timing}",
    "\\begin{threeparttable}",
    "\\small",
    "\\begin{tabular}{lrrr}",
    "\\toprule",
    "Computational component & SIGA-S & SIGA-R & RT \\\\ ",
    "\\midrule",
    sprintf("One-time allocation calibration & %.2f & %.2f & -- \\\\ ", s_cal, pair_minutes),
    sprintf("Complete unadjusted and adjusted analyses & %.2f & %.2f & %.2f \\\\ ", s_analysis, siga_r_analysis, rt_analysis),
    sprintf("Total method-specific computation time & %.2f & %.2f & %.2f \\\\ ", s_cal + s_analysis, pair_minutes + siga_r_analysis, rt_analysis),
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{tablenotes}[flushleft]",
    "\\footnotesize",
    "\\item Times are minutes for all 100,000 outer trials. Data generation and score construction, which are common to the procedures, are excluded. The SIGA-R calibration is the reusable pair-path calibration.",
    "\\end{tablenotes}",
    "\\end{threeparttable}",
    "\\end{table}",
    ""
  )
  writeLines(lines, path, useBytes = TRUE)
}

aggregate_extension <- function() {
  dat <- read_all_extension_shards()
  summary <- rbind(
    summarize_analysis(dat, "Unadjusted"),
    summarize_analysis(dat, "Baseline-adjusted")
  )

  max_reproduction_error <- max(summary$max_abs_siga_s_pvalue_reproduction_error)
  if (!is.finite(max_reproduction_error) || max_reproduction_error > REPRODUCTION_TOLERANCE) {
    audit_path <- file.path(PUBLICATION_DIR, "swift_direct_siga_r_reproduction_failure.csv")
    safe_write_csv(summary, audit_path)
    stop(
      "The regenerated SIGA-S p-values did not reproduce the existing results. ",
      "Maximum absolute error = ", format(max_reproduction_error, scientific = TRUE),
      "; tolerance = ", format(REPRODUCTION_TOLERANCE, scientific = TRUE),
      ". Inspect ", audit_path, call. = FALSE
    )
  }

  summary_path <- file.path(PUBLICATION_DIR, "swift_direct_siga_r_summary.csv")
  detail_path <- file.path(PUBLICATION_DIR, "swift_direct_siga_r_detail.csv")
  table_path <- file.path(PUBLICATION_DIR, "swift_direct_siga_r_table.tex")
  timing_path <- file.path(PUBLICATION_DIR, "swift_direct_siga_r_timing.tex")
  paragraph_path <- file.path(PUBLICATION_DIR, "swift_direct_siga_r_results_paragraph.tex")
  audit_path <- file.path(PUBLICATION_DIR, "swift_direct_siga_r_audit.csv")

  safe_write_csv(summary, summary_path)
  safe_write_csv(dat, detail_path)
  write_results_table(summary, table_path)
  pair <- get_pair_calibration(create = FALSE)
  write_timing_table(dat, pair, timing_path)

  u <- summary[summary$analysis == "Unadjusted", ]
  a <- summary[summary$analysis == "Baseline-adjusted", ]
  paragraph <- paste0(
    "At the planned total sample size of 404, estimated non-inferiority power for the unadjusted score was ",
    sprintf("%.2f", 100 * u$siga_s_power), "\\% for SIGA-S, ",
    sprintf("%.2f", 100 * u$siga_r_power), "\\% for SIGA-R and ",
    sprintf("%.2f", 100 * u$rt_power), "\\% for the reference randomization test. ",
    "With baseline adjustment, the corresponding estimates were ",
    sprintf("%.2f", 100 * a$siga_s_power), "\\%, ",
    sprintf("%.2f", 100 * a$siga_r_power), "\\% and ",
    sprintf("%.2f", 100 * a$rt_power), "\\%, respectively. ",
    "The mean raw variance ratios $\\widehat\\rho_n$ were ",
    sprintf("%.3f", u$mean_variance_ratio), " and ",
    sprintf("%.3f", a$mean_variance_ratio),
    " for the unadjusted and adjusted scores."
  )
  writeLines(paragraph, paragraph_path, useBytes = TRUE)

  old_cal <- load_old_calibration()
  gamma_difference <- max(abs(old_cal$gamma - pair$gamma))
  audit <- data.frame(
    item = c(
      "completed outer trials",
      "maximum absolute SIGA-S p-value reproduction error",
      "maximum absolute difference between old and pair-calibration gamma entries",
      "pair-calibration replicates",
      "pair-calibration elapsed minutes",
      "SIGA-R safeguard rate, unadjusted",
      "SIGA-R safeguard rate, adjusted"
    ),
    value = c(
      nrow(dat),
      max_reproduction_error,
      gamma_difference,
      pair$B0,
      pair$elapsed_seconds / 60,
      u$safeguard_rate,
      a$safeguard_rate
    ),
    stringsAsFactors = FALSE
  )
  safe_write_csv(audit, audit_path)

  message("Aggregation complete.")
  message("Summary: ", summary_path)
  message("Table: ", table_path)
  message("Timing: ", timing_path)
  message("Paragraph: ", paragraph_path)
  invisible(summary)
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

message("SWIFT DIRECT SIGA-R case-study extension")
message("  mode                       : ", RUN_MODE)
message("  existing output            : ", OLD_OUTPUT_DIR)
message("  extension output           : ", OUTPUT_DIR)
message("  outer trials               : ", N_OUTER)
message("  pair calibration paths     : ", N_PAIR_CALIBRATION)
message("  biased-coin probability    : ", P_BIASED_COIN)
message("  base seed                  : ", BASE_SEED)
message("  epsilon_n                  : n^(-", EPSILON_EXPONENT, ")")

if (RUN_MODE == "calibrate") {
  invisible(load_old_calibration())
  invisible(get_pair_calibration(create = TRUE))
  message("Calibration complete.")
} else if (RUN_MODE == "run") {
  if (!file.exists(pair_calibration_file())) {
    stop("Run SWIFT_R_MODE=calibrate before launching shards.", call. = FALSE)
  }
  run_extension_shard()
} else {
  aggregate_extension()
}
