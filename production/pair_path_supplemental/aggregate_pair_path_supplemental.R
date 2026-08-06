#!/usr/bin/env Rscript

# Re-aggregate the completed 12-scenario / 120-shard supplemental SIGA run and create the portrait Supplementary Table 5.
# This script is self-contained: it does not source the simulation engine and
# does not rerun any outer trial or any 4,999-path randomization test.

options(stringsAsFactors = FALSE, warn = 1)

N_OUTER <- 100000L
N_RERANDOMIZATIONS <- 4999L
N_CALIBRATION <- 100000L
N_DIRECTION_CALIBRATION <- 200000L
N_SHARDS <- 120L
BASE_SEED <- 20260805L
RUN_VERSION <- "v20260805_supplemental_v2"
INCLUDE_ORACLE_N <- TRUE
RT_SHAPE_DIAGNOSTICS <- TRUE
CALIBRATION_SCHEMA_VERSION <- "paircal_v2_main_design_sizes"
VALIDATION_ABS_TOLERANCE <- 0.003
NOMINAL_ABS_TOLERANCE <- 0.004
VARIANCE_REL_TOLERANCE <- 0.12

script_directory <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

BASE_DIR <- script_directory()
run_tag <- paste0(
  "profile_manuscript_set_supplemental_",
  "M100000_B4999_B0100000_Bdir200000_eps1p00_",
  "oracle1_shape1_seed20260805_",
  "ver_v20260805_supplemental_v2_S120"
)
OUTPUT_DIR <- path.expand(Sys.getenv("PWRT_RUN_DIR", unset = ""))
if (!nzchar(OUTPUT_DIR)) {
  project_dir <- path.expand(Sys.getenv(
    "PWRT_PROJECT_DIR",
    unset = file.path(dirname(BASE_DIR), "pair_path_supplemental")
  ))
  output_root <- path.expand(Sys.getenv(
    "PWRT_OUTPUT_DIR",
    unset = file.path(project_dir, "master")
  ))
  OUTPUT_DIR <- file.path(output_root, run_tag)
}
SCENARIO_DIR <- file.path(OUTPUT_DIR, "scenario_shards")

if (!dir.exists(SCENARIO_DIR)) {
  stop("Scenario-shard directory is missing: ", SCENARIO_DIR, call. = FALSE)
}

safe_write_csv <- function(data, path) {
  parent <- dirname(path)
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create output directory: ", parent, call. = FALSE)
  }
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

to_logical <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x) || is.integer(x)) return(ifelse(is.na(x), NA, x != 0))
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[y %in% c("false", "f", "0", "no", "n")] <- FALSE
  out
}

flag01 <- function(x) {
  y <- to_logical(x)
  as.integer(y)
}

safe_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

safe_sum <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else sum(x)
}

safe_sd <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else sd(x)
}

safe_quantile <- function(x, probability) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else unname(quantile(x, probability, names = FALSE))
}

first_value <- function(dat, name, default = NA) {
  if (!name %in% names(dat) || !length(dat[[name]])) return(default)
  dat[[name]][1L]
}

require_columns <- function(dat, columns, context) {
  missing <- setdiff(columns, names(dat))
  if (length(missing)) {
    stop(
      context, " is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

wilson_interval <- function(x, n, conf.level = 0.95) {
  if (!is.finite(n) || n <= 0) {
    return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  z <- qnorm(1 - (1 - conf.level) / 2)
  phat <- x / n
  denominator <- 1 + z^2 / n
  center <- (phat + z^2 / (2 * n)) / denominator
  half <- z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2)) / denominator
  c(
    estimate = phat,
    lower = max(0, center - half),
    upper = min(1, center + half)
  )
}

rejection_summary <- function(x) {
  x <- to_logical(x)
  x <- x[!is.na(x)]
  if (!length(x)) {
    return(c(probability = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  ci <- wilson_interval(sum(x), length(x))
  c(
    probability = unname(ci["estimate"]),
    lower = unname(ci["lower"]),
    upper = unname(ci["upper"])
  )
}

paired_difference_interval <- function(x, y, conf.level = 0.95) {
  x <- as.numeric(to_logical(x))
  y <- as.numeric(to_logical(y))
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (!length(x)) {
    return(c(estimate = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  delta <- x - y
  estimate <- mean(delta)
  se <- if (length(delta) > 1L) sd(delta) / sqrt(length(delta)) else NA_real_
  z <- qnorm(1 - (1 - conf.level) / 2)
  c(
    estimate = estimate,
    se = se,
    lower = estimate - z * se,
    upper = estimate + z * se
  )
}

safe_ratio <- function(numerator, denominator) {
  out <- numerator / denominator
  out[!is.finite(out)] <- NA_real_
  out
}

central_moment_summary <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 4L) {
    return(c(mean = NA_real_, variance = NA_real_, skewness = NA_real_, excess_kurtosis = NA_real_))
  }
  center <- x - mean(x)
  second <- mean(center^2)
  if (!is.finite(second) || second <= 0) {
    return(c(mean = mean(x), variance = var(x), skewness = NA_real_, excess_kurtosis = NA_real_))
  }
  c(
    mean = mean(x),
    variance = var(x),
    skewness = mean(center^3) / second^(3 / 2),
    excess_kurtosis = mean(center^4) / second^2 - 3
  )
}

read_scenario_details <- function(directory) {
  files <- list.files(
    directory,
    pattern = "^shard_[0-9]+_of_120[.]csv$",
    full.names = TRUE
  )
  if (length(files) != N_SHARDS) {
    stop(
      basename(directory), " contains ", length(files),
      " shard files; expected ", N_SHARDS, ".",
      call. = FALSE
    )
  }

  pieces <- lapply(files, function(path) {
    dat <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    dat$source_file <- basename(path)
    dat
  })
  dat <- do.call(rbind, pieces)

  required_metadata <- c(
    "scenario_id", "configured_outer_trials", "rerandomizations_per_trial",
    "allocation_calibration_paths", "direction_calibration_paths",
    "epsilon_exponent", "configured_shards", "base_seed", "run_version",
    "include_oracle_n", "rt_shape_diagnostics", "calibration_schema_version",
    "machine_id", "replicate"
  )
  require_columns(dat, required_metadata, basename(directory))

  expected_numeric <- list(
    configured_outer_trials = N_OUTER,
    rerandomizations_per_trial = N_RERANDOMIZATIONS,
    allocation_calibration_paths = N_CALIBRATION,
    direction_calibration_paths = N_DIRECTION_CALIBRATION,
    epsilon_exponent = 1.0,
    configured_shards = N_SHARDS,
    base_seed = BASE_SEED
  )
  for (field in names(expected_numeric)) {
    observed <- unique(as.numeric(dat[[field]]))
    observed <- observed[is.finite(observed)]
    if (length(observed) != 1L || !isTRUE(all.equal(observed, as.numeric(expected_numeric[[field]])))) {
      stop(
        "Incompatible metadata in ", basename(directory), ": ", field,
        "=", paste(observed, collapse = ","),
        "; expected ", expected_numeric[[field]],
        call. = FALSE
      )
    }
  }

  observed_oracle <- unique(flag01(dat$include_oracle_n))
  observed_oracle <- observed_oracle[is.finite(observed_oracle)]
  if (length(observed_oracle) != 1L || observed_oracle != 1L) {
    stop("include_oracle_n is not consistently enabled in ", basename(directory), call. = FALSE)
  }
  observed_shape <- unique(flag01(dat$rt_shape_diagnostics))
  observed_shape <- observed_shape[is.finite(observed_shape)]
  if (length(observed_shape) != 1L || observed_shape != 1L) {
    stop("rt_shape_diagnostics is not consistently enabled in ", basename(directory), call. = FALSE)
  }

  observed_version <- unique(as.character(dat$run_version))
  if (length(observed_version) != 1L || !identical(observed_version, RUN_VERSION)) {
    stop(
      "Incompatible run_version in ", basename(directory), ": ",
      paste(observed_version, collapse = ","),
      call. = FALSE
    )
  }
  observed_schema <- unique(as.character(dat$calibration_schema_version))
  if (length(observed_schema) != 1L || !identical(observed_schema, CALIBRATION_SCHEMA_VERSION)) {
    stop(
      "Incompatible calibration_schema_version in ", basename(directory), ": ",
      paste(observed_schema, collapse = ","),
      call. = FALSE
    )
  }

  dat$replicate <- as.integer(dat$replicate)
  dat <- dat[order(dat$replicate), , drop = FALSE]

  duplicated_id <- duplicated(dat$replicate)
  if (any(duplicated_id)) {
    duplicate_values <- unique(dat$replicate[duplicated_id])
    compare_columns <- setdiff(names(dat), "source_file")
    for (id in duplicate_values) {
      block <- dat[dat$replicate == id, compare_columns, drop = FALSE]
      if (nrow(unique(block)) > 1L) {
        stop(
          "Conflicting duplicate replicate ", id,
          " in ", basename(directory),
          call. = FALSE
        )
      }
    }
    dat <- dat[!duplicated(dat$replicate), , drop = FALSE]
  }

  if (nrow(dat) != N_OUTER) {
    stop(
      basename(directory), " contains ", nrow(dat),
      " unique replicates; expected ", N_OUTER, ".",
      call. = FALSE
    )
  }
  expected_ids <- seq_len(N_OUTER)
  if (!identical(dat$replicate, expected_ids)) {
    missing_ids <- setdiff(expected_ids, dat$replicate)
    stop(
      "Replicate sequence is incomplete in ", basename(directory),
      if (length(missing_ids)) paste0("; first missing replicate: ", missing_ids[1L]) else "",
      call. = FALSE
    )
  }

  dat
}

summarize_analysis <- function(dat, analysis, alpha) {
  suffix <- if (analysis == "unadjusted") "unadjusted" else "adjusted"
  needed <- paste0(
    c(
      "reject_siga_s_", "reject_siga_r_", "reject_oracle_n_", "reject_rt_",
      "statistic_", "variance_s_", "variance_r_", "variance_rt_",
      "p_siga_s_", "p_siga_r_", "p_oracle_n_", "p_rt_",
      "randomization_skewness_", "randomization_excess_kurtosis_",
      "kappa_raw_", "kappa_truncated_", "safeguard_active_"
    ),
    suffix
  )
  require_columns(
    dat,
    c(
      needed, "empty_joint_strata", "min_positive_joint_stratum_count"
    ),
    paste0("scenario ", first_value(dat, "scenario_code", "?"), "/", analysis)
  )

  reject_s <- dat[[paste0("reject_siga_s_", suffix)]]
  reject_r <- dat[[paste0("reject_siga_r_", suffix)]]
  reject_o <- dat[[paste0("reject_oracle_n_", suffix)]]
  reject_rt <- dat[[paste0("reject_rt_", suffix)]]

  r_s <- rejection_summary(reject_s)
  r_r <- rejection_summary(reject_r)
  r_o <- rejection_summary(reject_o)
  r_rt <- rejection_summary(reject_rt)

  diff_s_rt <- paired_difference_interval(reject_s, reject_rt)
  diff_r_rt <- paired_difference_interval(reject_r, reject_rt)
  diff_o_rt <- paired_difference_interval(reject_o, reject_rt)
  diff_r_o <- paired_difference_interval(reject_r, reject_o)

  statistic <- as.numeric(dat[[paste0("statistic_", suffix)]])
  v_s <- as.numeric(dat[[paste0("variance_s_", suffix)]])
  v_r <- as.numeric(dat[[paste0("variance_r_", suffix)]])
  v_rt <- as.numeric(dat[[paste0("variance_rt_", suffix)]])
  p_s <- as.numeric(dat[[paste0("p_siga_s_", suffix)]])
  p_r <- as.numeric(dat[[paste0("p_siga_r_", suffix)]])
  p_o <- as.numeric(dat[[paste0("p_oracle_n_", suffix)]])
  p_rt <- as.numeric(dat[[paste0("p_rt_", suffix)]])

  ratio_r_s <- safe_ratio(v_r, v_s)
  ratio_s_rt <- safe_ratio(v_s, v_rt)
  ratio_r_rt <- safe_ratio(v_r, v_rt)

  stat_moments <- central_moment_summary(statistic)
  empirical_var_stat <- unname(stat_moments["variance"])
  mean_v_s <- safe_mean(v_s)
  mean_v_r <- safe_mean(v_r)
  mean_v_rt <- safe_mean(v_rt)
  target_ratio <- mean_v_rt / empirical_var_stat
  estimated_ratio <- mean_v_r / mean_v_s
  zcrit <- qnorm(1 - alpha)
  predicted_rt_from_target_ratio <- if (is.finite(target_ratio) && target_ratio > 0) {
    1 - pnorm(zcrit * sqrt(target_ratio))
  } else {
    NA_real_
  }
  trial_prediction <- 1 - pnorm(zcrit * sqrt(ratio_r_s))

  z_s <- statistic / sqrt(v_s)
  z_r <- statistic / sqrt(v_r)
  z_s_moments <- central_moment_summary(z_s)
  z_r_moments <- central_moment_summary(z_r)

  rt_skewness <- as.numeric(dat[[paste0("randomization_skewness_", suffix)]])
  rt_kurtosis <- as.numeric(dat[[paste0("randomization_excess_kurtosis_", suffix)]])

  se_for_gap <- c(diff_r_o["se"], diff_o_rt["se"])
  se_for_gap <- se_for_gap[is.finite(se_for_gap)]
  gap_tolerance <- max(
    VALIDATION_ABS_TOLERANCE,
    if (length(se_for_gap)) 3 * max(se_for_gap) else 0
  )
  variance_gap <- is.finite(diff_r_o["estimate"]) && abs(diff_r_o["estimate"]) > gap_tolerance
  tail_gap <- is.finite(diff_o_rt["estimate"]) && abs(diff_o_rt["estimate"]) > gap_tolerance
  gap_source <- if (variance_gap && tail_gap) {
    "mixed_variance_and_tail"
  } else if (variance_gap) {
    "variance_approximation"
  } else if (tail_gap) {
    "non_gaussian_or_discrete_tail"
  } else {
    "first_order_alignment"
  }

  data.frame(
    analysis = analysis,
    siga_s = unname(r_s["probability"]),
    siga_s_lower_95 = unname(r_s["lower"]),
    siga_s_upper_95 = unname(r_s["upper"]),
    siga_r = unname(r_r["probability"]),
    siga_r_lower_95 = unname(r_r["lower"]),
    siga_r_upper_95 = unname(r_r["upper"]),
    oracle_n = unname(r_o["probability"]),
    oracle_n_lower_95 = unname(r_o["lower"]),
    oracle_n_upper_95 = unname(r_o["upper"]),
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
    difference_oracle_n_minus_rt = unname(diff_o_rt["estimate"]),
    difference_oracle_n_minus_rt_se = unname(diff_o_rt["se"]),
    difference_oracle_n_minus_rt_lower_95 = unname(diff_o_rt["lower"]),
    difference_oracle_n_minus_rt_upper_95 = unname(diff_o_rt["upper"]),
    difference_siga_r_minus_oracle_n = unname(diff_r_o["estimate"]),
    difference_siga_r_minus_oracle_n_se = unname(diff_r_o["se"]),
    difference_siga_r_minus_oracle_n_lower_95 = unname(diff_r_o["lower"]),
    difference_siga_r_minus_oracle_n_upper_95 = unname(diff_r_o["upper"]),
    empirical_variance_statistic = empirical_var_stat,
    mean_variance_s = mean_v_s,
    mean_variance_r = mean_v_r,
    mean_variance_rt = mean_v_rt,
    relative_bias_variance_s = mean_v_s / empirical_var_stat - 1,
    relative_bias_variance_r_vs_rt = mean_v_r / mean_v_rt - 1,
    target_variance_ratio_rt_over_sampling = target_ratio,
    estimated_variance_ratio_r_over_s = estimated_ratio,
    mean_variance_ratio_r_over_s = safe_mean(ratio_r_s),
    sd_variance_ratio_r_over_s = safe_sd(ratio_r_s),
    q025_variance_ratio_r_over_s = safe_quantile(ratio_r_s, 0.025),
    median_variance_ratio_r_over_s = safe_quantile(ratio_r_s, 0.5),
    q975_variance_ratio_r_over_s = safe_quantile(ratio_r_s, 0.975),
    mean_variance_ratio_s_over_rt = safe_mean(ratio_s_rt),
    mean_variance_ratio_r_over_rt = safe_mean(ratio_r_rt),
    predicted_rt_size_from_target_ratio = predicted_rt_from_target_ratio,
    mean_trial_level_normal_prediction = safe_mean(trial_prediction),
    mean_abs_p_difference_siga_s_rt = safe_mean(abs(p_s - p_rt)),
    mean_abs_p_difference_siga_r_rt = safe_mean(abs(p_r - p_rt)),
    mean_abs_p_difference_oracle_n_rt = safe_mean(abs(p_o - p_rt)),
    p95_abs_p_difference_siga_s_rt = safe_quantile(abs(p_s - p_rt), 0.95),
    p95_abs_p_difference_siga_r_rt = safe_quantile(abs(p_r - p_rt), 0.95),
    p95_abs_p_difference_oracle_n_rt = safe_quantile(abs(p_o - p_rt), 0.95),
    statistic_mean = unname(stat_moments["mean"]),
    statistic_skewness = unname(stat_moments["skewness"]),
    statistic_excess_kurtosis = unname(stat_moments["excess_kurtosis"]),
    z_s_mean = unname(z_s_moments["mean"]),
    z_s_variance = unname(z_s_moments["variance"]),
    z_s_skewness = unname(z_s_moments["skewness"]),
    z_s_excess_kurtosis = unname(z_s_moments["excess_kurtosis"]),
    z_r_mean = unname(z_r_moments["mean"]),
    z_r_variance = unname(z_r_moments["variance"]),
    z_r_skewness = unname(z_r_moments["skewness"]),
    z_r_excess_kurtosis = unname(z_r_moments["excess_kurtosis"]),
    mean_rt_randomization_skewness = safe_mean(rt_skewness),
    mean_rt_randomization_excess_kurtosis = safe_mean(rt_kurtosis),
    correlation_zs2_rho = suppressWarnings(cor(z_s^2, ratio_r_s, use = "complete.obs")),
    mean_empty_joint_strata = safe_mean(dat$empty_joint_strata),
    probability_any_empty_joint_stratum = safe_mean(as.numeric(dat$empty_joint_strata > 0)),
    mean_min_positive_joint_stratum_count = safe_mean(dat$min_positive_joint_stratum_count),
    kappa_raw_mean = safe_mean(dat[[paste0("kappa_raw_", suffix)]]),
    kappa_truncation_rate = safe_mean(as.numeric(to_logical(dat[[paste0("kappa_truncated_", suffix)]]))),
    safeguard_rate = safe_mean(as.numeric(to_logical(dat[[paste0("safeguard_active_", suffix)]]))),
    gap_explanation = gap_source,
    stringsAsFactors = FALSE
  )
}

summarize_scenario <- function(dat) {
  alpha <- as.numeric(first_value(dat, "alpha", 0.025))
  by_analysis <- rbind(
    summarize_analysis(dat, "unadjusted", alpha),
    summarize_analysis(dat, "adjusted", alpha)
  )

  machine_ids <- sort(unique(as.character(dat$machine_id)))
  machine_ids <- machine_ids[nzchar(machine_ids)]

  fixed <- data.frame(
    scenario_id = as.integer(first_value(dat, "scenario_id")),
    scenario_code = as.character(first_value(dat, "scenario_code")),
    scenario_label = as.character(first_value(dat, "scenario_label")),
    scenario_class = as.character(first_value(dat, "scenario_class")),
    set = as.character(first_value(dat, "set")),
    design_id = as.integer(first_value(dat, "design_id")),
    design_label = as.character(first_value(dat, "design_label")),
    factor_count = as.integer(first_value(dat, "factor_count")),
    target_per_group = as.integer(first_value(dat, "target_per_group")),
    total_n = as.integer(first_value(dat, "total_n")),
    pbc = as.numeric(first_value(dat, "pbc")),
    profile_type = as.character(first_value(dat, "profile_type")),
    outcome_type = as.character(first_value(dat, "outcome_type")),
    model_type = as.character(first_value(dat, "model_type")),
    direction = as.character(first_value(dat, "direction")),
    boundary = as.numeric(first_value(dat, "boundary")),
    alpha = alpha,
    achieved_effect = as.numeric(first_value(dat, "achieved_effect")),
    pair_ratio_model = as.numeric(first_value(dat, "pair_ratio_model")),
    pair_scalar_gap_model = as.numeric(first_value(dat, "pair_scalar_gap_model")),
    n_outer = nrow(dat),
    rerandomizations_per_trial = N_RERANDOMIZATIONS,
    allocation_calibration_paths = N_CALIBRATION,
    direction_calibration_paths = N_DIRECTION_CALIBRATION,
    base_seed = BASE_SEED,
    run_version = RUN_VERSION,
    include_oracle_n = INCLUDE_ORACLE_N,
    rt_shape_diagnostics = RT_SHAPE_DIAGNOSTICS,
    machine_count = length(machine_ids),
    machine_ids = paste(machine_ids, collapse = ";"),
    mean_realized_pair_scalar_gap = safe_mean(dat$pair_scalar_gap_realized),
    mean_observed_treated = safe_mean(dat$observed_treated),
    data_generation_minutes = safe_sum(dat$data_generation_seconds) / 60,
    score_construction_minutes = safe_sum(dat$score_construction_seconds) / 60,
    siga_analysis_minutes = safe_sum(dat$siga_analysis_seconds) / 60,
    rt_analysis_minutes = safe_sum(dat$rt_analysis_seconds) / 60,
    stringsAsFactors = FALSE
  )

  cbind(fixed[rep(1L, nrow(by_analysis)), , drop = FALSE], by_analysis)
}

write_latex_tables <- function(summary, output_directory) {
  fmt <- function(x, digits = 2) {
    ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "--")
  }
  clean <- function(x) gsub("_", "-", as.character(x), fixed = TRUE)
  lines <- c(
    "\\begingroup",
    "\\fontsize{8.5pt}{10.0pt}\\selectfont",
    "\\setlength{\\tabcolsep}{2.0pt}",
    "\\renewcommand{\\arraystretch}{1.06}",
    "\\setlength{\\LTleft}{\\fill}",
    "\\setlength{\\LTright}{\\fill}",
    paste0(
      "\\begin{longtable}{@{}>{\\centering\\arraybackslash}m{28pt}",
      ">{\\centering\\arraybackslash}m{18pt}",
      ">{\\centering\\arraybackslash}m{34pt}",
      ">{\\raggedright\\arraybackslash}m{48pt}",
      ">{\\raggedright\\arraybackslash}m{66pt}",
      ">{\\raggedright\\arraybackslash}m{52pt}",
      ">{\\centering\\arraybackslash}m{40pt}",
      "*{3}{>{\\centering\\arraybackslash}m{31pt}}@{}}"
    ),
    "\\caption{Rejection probabilities in the targeted pair-path sensitivity analysis.}",
    "\\label{tab:supp-pair-path-sensitivity-results}\\\\",
    "\\toprule",
    paste0(
      "Code & $F$ & \\shortstack{$n$/\\\\group} & Outcome & Model & Analysis & ",
      "\\shortstack{Mean\\\\$\\widehat\\rho_n$} & SIGA-S & SIGA-R & RT \\\\"
    ),
    "\\midrule",
    "\\endfirsthead",
    "\\multicolumn{10}{c}{\\tablename\\ \\thetable{} -- continued}\\\\",
    "\\toprule",
    paste0(
      "Code & $F$ & \\shortstack{$n$/\\\\group} & Outcome & Model & Analysis & ",
      "\\shortstack{Mean\\\\$\\widehat\\rho_n$} & SIGA-S & SIGA-R & RT \\\\"
    ),
    "\\midrule",
    "\\endhead",
    "\\midrule",
    "\\multicolumn{10}{r}{Continued on next page}\\\\",
    "\\endfoot",
    "\\bottomrule",
    "\\endlastfoot"
  )
  analysis_order <- match(summary$analysis, c("unadjusted", "adjusted"))
  summary <- summary[order(summary$scenario_id, analysis_order), , drop = FALSE]
  for (i in seq_len(nrow(summary))) {
    row <- summary[i, ]
    lines <- c(lines, paste0(
      row$scenario_code, " & ", row$factor_count, " & ", row$target_per_group,
      " & ", tools::toTitleCase(row$outcome_type),
      " & ", clean(row$model_type),
      " & ", tools::toTitleCase(row$analysis),
      " & ", fmt(row$mean_variance_ratio_r_over_s, 3),
      " & ", fmt(100 * row$siga_s, 2),
      " & ", fmt(100 * row$siga_r, 2),
      " & ", fmt(100 * row$rt, 2), " \\\\"
    ))
  }
  lines <- c(
    lines,
    "\\end{longtable}",
    paste0(
      "\\noindent{\\footnotesize Rejection probabilities are percentages at a nominal one-sided level of 2.5\\%. ",
      "Codes S01--S08 denote practically heterogeneous scenarios under $p_{\\rm bc}=0.80$, ",
      "whereas S09--S12 denote deliberately strong pair-path scenarios under $p_{\\rm bc}=0.95$.}"
    ),
    "\\endgroup"
  )
  path <- file.path(output_directory, "pair_path_supplemental_table.tex")
  writeLines(lines, path)
  invisible(path)
}

apply_validation_gates <- function(summary) {
  rows <- lapply(seq_len(nrow(summary)), function(i) {
    row <- summary[i, ]
    se_candidates <- c(
      row$difference_siga_r_minus_rt_se,
      row$difference_siga_r_minus_oracle_n_se,
      row$difference_oracle_n_minus_rt_se
    )
    se_candidates <- se_candidates[is.finite(se_candidates)]
    paired_tolerance <- max(
      VALIDATION_ABS_TOLERANCE,
      if (length(se_candidates)) 3 * max(se_candidates) else 0
    )
    r_close_rt <- abs(row$difference_siga_r_minus_rt) <= paired_tolerance
    r_close_oracle <- abs(row$difference_siga_r_minus_oracle_n) <= paired_tolerance
    oracle_close_rt <- abs(row$difference_oracle_n_minus_rt) <= paired_tolerance
    rt_agreement_or_explained <- r_close_rt || r_close_oracle
    s_nominal <- abs(row$siga_s - row$alpha) <= NOMINAL_ABS_TOLERANCE
    variance_r_close <- abs(row$relative_bias_variance_r_vs_rt) <= VARIANCE_REL_TOLERANCE
    variance_s_close <- abs(row$relative_bias_variance_s) <= VARIANCE_REL_TOLERANCE
    safeguard_ok <- row$safeguard_rate <= 0.001

    direction_ok <- TRUE
    if (
      row$scenario_class == "strong_pair_path" &&
      row$pair_ratio_model > 1.05
    ) {
      direction_ok <- row$siga_r <= row$siga_s && row$rt <= row$siga_s
    }

    primary_pass <- if (row$set == "supplemental") {
      s_nominal && variance_r_close && variance_s_close && safeguard_ok &&
        direction_ok && rt_agreement_or_explained
    } else {
      NA
    }

    data.frame(
      scenario_id = row$scenario_id,
      scenario_code = row$scenario_code,
      scenario_class = row$scenario_class,
      analysis = row$analysis,
      paired_tolerance = paired_tolerance,
      r_close_rt = r_close_rt,
      r_close_oracle = r_close_oracle,
      oracle_close_rt = oracle_close_rt,
      rt_agreement_or_explained = rt_agreement_or_explained,
      s_nominal = s_nominal,
      variance_r_close_to_rt = variance_r_close,
      variance_s_close_to_sampling = variance_s_close,
      safeguard_ok = safeguard_ok,
      direction_ok = direction_ok,
      gap_explanation = row$gap_explanation,
      primary_pass = primary_pass,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

shard_files <- list.files(
  SCENARIO_DIR,
  pattern = "^shard_[0-9]+_of_120[.]csv$",
  recursive = TRUE,
  full.names = TRUE
)
if (length(shard_files) != 1440L) {
  stop("Expected 1440 scenario-shard files, but found ", length(shard_files), ".", call. = FALSE)
}
scenario_directories <- sort(unique(dirname(shard_files)))
if (length(scenario_directories) != 12L) {
  stop("Expected 12 scenario directories, but found ", length(scenario_directories), ".", call. = FALSE)
}
cat("Confirmed 12 scenarios and 1440 completed scenario-shard files.\n")

output_names <- c(
  "pair_path_supplemental_summary.csv",
  "pair_path_supplemental_gates.csv",
  "pair_path_supplemental_table.tex"
)
existing_outputs <- file.path(OUTPUT_DIR, output_names)
existing_outputs <- existing_outputs[file.exists(existing_outputs)]
if (length(existing_outputs)) {
  backup_dir <- file.path(
    OUTPUT_DIR,
    paste0("backup_before_reaggregation_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )
  if (!dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE) && !dir.exists(backup_dir)) {
    stop("Could not create backup directory: ", backup_dir, call. = FALSE)
  }
  copied <- file.copy(existing_outputs, backup_dir, overwrite = FALSE)
  if (!all(copied)) stop("Could not back up all existing summary/table files.", call. = FALSE)
  cat("Previous summary/table files backed up to:\n", backup_dir, "\n", sep = "")
}

all_summary <- vector("list", length(scenario_directories))
for (i in seq_along(scenario_directories)) {
  directory <- scenario_directories[i]
  cat(
    sprintf(
      "[%02d/12] Reading and aggregating %s\n",
      i, basename(directory)
    )
  )
  dat <- read_scenario_details(directory)
  all_summary[[i]] <- summarize_scenario(dat)
  rm(dat)
  invisible(gc(FALSE))
}

summary <- do.call(rbind, all_summary)
summary <- summary[order(summary$scenario_id, summary$analysis), , drop = FALSE]
validation <- apply_validation_gates(summary)

summary_path <- file.path(OUTPUT_DIR, "pair_path_supplemental_summary.csv")
gates_path <- file.path(OUTPUT_DIR, "pair_path_supplemental_gates.csv")
safe_write_csv(summary, summary_path)
safe_write_csv(validation, gates_path)
write_latex_tables(summary, OUTPUT_DIR)

expected_outputs <- file.path(OUTPUT_DIR, output_names)
missing_outputs <- expected_outputs[!file.exists(expected_outputs)]
if (length(missing_outputs)) {
  stop(
    "Aggregation finished, but expected output files are missing:\n",
    paste(missing_outputs, collapse = "\n"),
    call. = FALSE
  )
}

failed <- validation$primary_pass %in% FALSE
cat("\nRe-aggregation completed. No outer simulation or RT calculation was rerun.\n")
cat("Updated files:\n", paste(expected_outputs, collapse = "\n"), "\n", sep = "")
if (any(failed)) {
  cat(
    "\nSupplemental validation gates marked as failed: ",
    paste0(validation$scenario_code[failed], "/", validation$analysis[failed], collapse = ", "),
    "\n",
    sep = ""
  )
} else {
  cat("\nAll supplemental primary validation gates passed.\n")
}
