# Reproducible workflows for the full-grid benchmarks and
# and SWIFT DIRECT-inspired simulation.
# Source siga_core.R and siga_study_definitions.R first.

analyse_boundary_scores <- function(trial,
                                    boundary,
                                    calibration,
                                    d,
                                    siga_target = c("S", "R"),
                                    alternative = c("two.sided", "greater", "less"),
                                    use_lnm_for_unadjusted = FALSE,
                                    epsilon = NULL) {
  siga_target <- match.arg(siga_target)
  alternative <- match.arg(alternative)
  scores <- make_boundary_scores(trial$y, trial$A, trial$X, boundary)
  statistic_unadjusted <- 0.5 * sum(trial$z * scores$unadjusted)
  statistic_adjusted <- 0.5 * sum(trial$z * scores$adjusted)
  vs_u <- siga_sampling_variance(scores$unadjusted, trial$X, calibration)
  vs_a <- siga_sampling_variance(scores$adjusted, trial$X, calibration)
  if (siga_target == "S") {
    variance_u <- vs_u$variance
    variance_a <- vs_a$variance
    p_u <- if (use_lnm_for_unadjusted) {
      lattice_normal_mixture_pvalue(
        trial$y, trial$A, statistic_unadjusted, variance_u,
        calibration$treated_count_prob, alternative
      )
    } else {
      gaussian_pvalue(statistic_unadjusted, variance_u, alternative)$p
    }
    p_a <- gaussian_pvalue(statistic_adjusted, variance_a, alternative)$p
    vr_u <- vr_a <- NULL
  } else {
    vr_u <- siga_randomization_variance(vs_u, d, calibration, epsilon)
    vr_a <- siga_randomization_variance(vs_a, d, calibration, epsilon)
    variance_u <- vr_u$variance
    variance_a <- vr_a$variance
    p_u <- gaussian_pvalue(statistic_unadjusted, variance_u, alternative)$p
    p_a <- gaussian_pvalue(statistic_adjusted, variance_a, alternative)$p
  }
  list(
    scores = scores,
    statistic_unadjusted = statistic_unadjusted,
    statistic_adjusted = statistic_adjusted,
    p_unadjusted = p_u,
    p_adjusted = p_a,
    sampling_unadjusted = vs_u,
    sampling_adjusted = vs_a,
    randomization_unadjusted = vr_u,
    randomization_adjusted = vr_a
  )
}

run_full_grid_replicate <- function(scenario,
                                    design,
                                    calibration,
                                    siga_target,
                                    B_randomization,
                                    base_seed,
                                    replicate_id) {
  model <- if (scenario$outcome == "continuous") {
    make_continuous_benchmark_model(design, scenario$effect)
  } else {
    make_binary_benchmark_model(design, scenario$effect)
  }
  trial <- generate_benchmark_trial(
    design, model,
    seed_profile = seed_value(base_seed, scenario$scenario_id, replicate_id, 1L),
    seed_allocation = seed_value(base_seed, scenario$scenario_id, replicate_id, 2L),
    seed_outcome = seed_value(base_seed, scenario$scenario_id, replicate_id, 3L)
  )
  epsilon <- 1 / design$total_n

  if (scenario$objective == "Superiority") {
    boundary <- 0
    result <- analyse_boundary_scores(
      trial, boundary, calibration, model$d_at(boundary), siga_target,
      alternative = "two.sided",
      use_lnm_for_unadjusted = siga_target == "S" && scenario$outcome == "binary",
      epsilon = epsilon
    )
    rt <- fixed_score_randomization_test(
      trial$X,
      cbind(unadjusted = result$scores$unadjusted,
            adjusted = result$scores$adjusted),
      trial$z, B = B_randomization, pbc = design$pbc,
      weights = design$weights,
      seed = seed_value(base_seed, scenario$scenario_id, replicate_id, 4L)
    )
    p_rt_u <- rt$two_sided[1L]
    p_rt_a <- rt$two_sided[2L]
    boundary_pair <- c(boundary, NA_real_)
    results_lower <- result
    results_upper <- NULL
  } else if (scenario$objective == "Non-inferiority") {
    boundary <- scenario$lower_boundary
    result <- analyse_boundary_scores(
      trial, boundary, calibration, model$d_at(boundary), siga_target,
      alternative = "greater", use_lnm_for_unadjusted = FALSE,
      epsilon = epsilon
    )
    rt <- fixed_score_randomization_test(
      trial$X,
      cbind(unadjusted = result$scores$unadjusted,
            adjusted = result$scores$adjusted),
      trial$z, B = B_randomization, pbc = design$pbc,
      weights = design$weights,
      seed = seed_value(base_seed, scenario$scenario_id, replicate_id, 4L)
    )
    p_rt_u <- rt$greater[1L]
    p_rt_a <- rt$greater[2L]
    boundary_pair <- c(boundary, NA_real_)
    results_lower <- result
    results_upper <- NULL
  } else {
    lower <- scenario$lower_boundary
    upper <- scenario$upper_boundary
    lower_result <- analyse_boundary_scores(
      trial, lower, calibration, model$d_at(lower), siga_target,
      alternative = "greater", use_lnm_for_unadjusted = FALSE,
      epsilon = epsilon
    )
    upper_result <- analyse_boundary_scores(
      trial, upper, calibration, model$d_at(upper), siga_target,
      alternative = "less", use_lnm_for_unadjusted = FALSE,
      epsilon = epsilon
    )
    rt <- fixed_score_randomization_test(
      trial$X,
      cbind(
        unadjusted_lower = lower_result$scores$unadjusted,
        unadjusted_upper = upper_result$scores$unadjusted,
        adjusted_lower = lower_result$scores$adjusted,
        adjusted_upper = upper_result$scores$adjusted
      ),
      trial$z, B = B_randomization, pbc = design$pbc,
      weights = design$weights,
      seed = seed_value(base_seed, scenario$scenario_id, replicate_id, 4L)
    )
    lower_result$p_unadjusted <- max(lower_result$p_unadjusted, upper_result$p_unadjusted)
    lower_result$p_adjusted <- max(lower_result$p_adjusted, upper_result$p_adjusted)
    p_rt_u <- max(rt$greater[1L], rt$less[2L])
    p_rt_a <- max(rt$greater[3L], rt$less[4L])
    boundary_pair <- c(lower, upper)
    results_lower <- lower_result
    results_upper <- upper_result
  }

  randomization_summary <- function(lower_object, upper_object, analysis) {
    if (siga_target == "S") {
      return(c(ratio_lower = NA_real_, ratio_upper = NA_real_, safeguard = FALSE))
    }
    slot <- if (analysis == "unadjusted") "randomization_unadjusted" else "randomization_adjusted"
    first <- lower_object[[slot]]
    second <- if (is.null(upper_object)) NULL else upper_object[[slot]]
    c(
      ratio_lower = first$variance_ratio,
      ratio_upper = if (is.null(second)) NA_real_ else second$variance_ratio,
      safeguard = first$safeguard_active || (!is.null(second) && second$safeguard_active)
    )
  }
  ratio_u <- randomization_summary(results_lower, results_upper, "unadjusted")
  ratio_a <- randomization_summary(results_lower, results_upper, "adjusted")

  data.frame(
    benchmark = if (siga_target == "S") "sampling_targeted" else "randomization_targeted",
    scenario_id = scenario$scenario_id,
    replicate = replicate_id,
    outcome = scenario$outcome,
    design_name = scenario$design_name,
    factor_count = scenario$factor_count,
    n_per_group = scenario$n_per_group,
    objective = scenario$objective,
    measure = scenario$measure,
    effect = scenario$effect,
    lower_boundary = boundary_pair[1L],
    upper_boundary = boundary_pair[2L],
    alpha = scenario$alpha,
    p_siga_unadjusted = results_lower$p_unadjusted,
    p_siga_adjusted = results_lower$p_adjusted,
    p_rt_unadjusted = p_rt_u,
    p_rt_adjusted = p_rt_a,
    reject_siga_unadjusted = results_lower$p_unadjusted <= scenario$alpha,
    reject_siga_adjusted = results_lower$p_adjusted <= scenario$alpha,
    reject_rt_unadjusted = p_rt_u <= scenario$alpha,
    reject_rt_adjusted = p_rt_a <= scenario$alpha,
    ratio_lower_unadjusted = ratio_u["ratio_lower"],
    ratio_upper_unadjusted = ratio_u["ratio_upper"],
    ratio_lower_adjusted = ratio_a["ratio_lower"],
    ratio_upper_adjusted = ratio_a["ratio_upper"],
    safeguard_unadjusted = as.logical(ratio_u["safeguard"]),
    safeguard_adjusted = as.logical(ratio_a["safeguard"]),
    stringsAsFactors = FALSE
  )
}

calibration_cache_name <- function(output_dir, target, design, B0, seed) {
  file.path(
    ensure_directory(file.path(output_dir, "calibration")),
    paste0("calibration_", target, "_", design$design_label,
           "_B0", B0, "_seed", normalise_seed(seed), ".rds")
  )
}

get_full_grid_calibration <- function(output_dir,
                                      target,
                                      design,
                                      B0,
                                      seed,
                                      batch_size = 1000L,
                                      progress = TRUE) {
  path <- calibration_cache_name(output_dir, target, design, B0, seed)
  if (file.exists(path)) return(readRDS(path))
  calibration <- if (target == "S") {
    calibrate_one_path_design(
      B0, design$total_n, design$patterns, design$profile_prob,
      design$pbc, design$weights, seed, batch_size, progress
    )
  } else {
    calibrate_three_path_design(
      B0, design$total_n, design$patterns, design$profile_prob,
      design$pbc, design$weights, seed, batch_size, progress
    )
  }
  saveRDS(calibration, path)
  calibration
}

run_full_grid_benchmark <- function(target = c("S", "R"),
                                    n_outer = 100000L,
                                    B_randomization = 4999L,
                                    B_calibration = 100000L,
                                    base_seed = NULL,
                                    output_dir = "output/full_grid",
                                    n_shards = 1L,
                                    shard_id = 1L,
                                    scenario_ids = NULL,
                                    batch_size = 10L,
                                    calibration_batch = 1000L,
                                    progress = TRUE) {
  target <- match.arg(target)
  n_outer <- as.integer(n_outer)
  B_randomization <- as.integer(B_randomization)
  B_calibration <- as.integer(B_calibration)
  n_shards <- as.integer(n_shards)
  shard_id <- as.integer(shard_id)
  batch_size <- as.integer(batch_size)
  base_seed <- as.integer(base_seed %||% if (target == "S") 20260725L else 20260801L)
  if (n_outer < 1L || B_randomization < 1L || B_calibration < 2L ||
      n_shards < 1L || shard_id < 1L || shard_id > n_shards) {
    stop("Invalid benchmark run settings.", call. = FALSE)
  }
  output_dir <- ensure_directory(output_dir)
  scenarios <- make_full_grid_scenarios()
  if (!is.null(scenario_ids)) {
    scenarios <- scenarios[scenarios$scenario_id %in% as.integer(scenario_ids), , drop = FALSE]
  }
  designs <- make_full_grid_designs()
  safe_write_csv(scenarios, file.path(output_dir, paste0("scenario_definitions_", target, ".csv")))

  calibrations <- list()
  for (design_name in unique(scenarios$design_name)) {
    design <- designs[[design_name]]
    calibration_seed <- seed_value(base_seed, design$design_id, 0L, 90L)
    calibrations[[design_name]] <- get_full_grid_calibration(
      output_dir, target, design, B_calibration, calibration_seed,
      calibration_batch, progress
    )
  }

  result_dir <- ensure_directory(file.path(output_dir, paste0("shards_", target)))
  for (row in seq_len(nrow(scenarios))) {
    scenario <- scenarios[row, ]
    design <- designs[[scenario$design_name]]
    calibration <- calibrations[[scenario$design_name]]
    result_path <- file.path(
      result_dir,
      sprintf("scenario_%03d_shard_%03d_of_%03d.csv",
              scenario$scenario_id, shard_id, n_shards)
    )
    completed <- integer()
    if (file.exists(result_path) && file.info(result_path)$size > 0) {
      existing <- tryCatch(read.csv(result_path), error = function(e) NULL)
      if (!is.null(existing) && "replicate" %in% names(existing)) {
        completed <- unique(as.integer(existing$replicate))
      }
    }
    assigned <- seq.int(shard_id, n_outer, by = n_shards)
    pending <- setdiff(assigned, completed)
    if (progress) message("scenario ", scenario$scenario_id, ": pending ", length(pending))
    if (!length(pending)) next
    for (from in seq.int(1L, length(pending), by = batch_size)) {
      ids <- pending[from:min(length(pending), from + batch_size - 1L)]
      batch <- do.call(rbind, lapply(ids, function(replicate_id) {
        run_full_grid_replicate(
          scenario, design, calibration, target,
          B_randomization, base_seed, replicate_id
        )
      }))
      append_csv(batch, result_path)
      if (progress) message("  wrote through replicate ", max(ids))
    }
  }
  invisible(output_dir)
}

summarize_full_grid <- function(output_dir,
                                target = c("S", "R"),
                                n_outer = NULL,
                                allow_partial = FALSE) {
  target <- match.arg(target)
  scenario_file <- file.path(output_dir, paste0("scenario_definitions_", target, ".csv"))
  scenarios <- read.csv(scenario_file, stringsAsFactors = FALSE)
  files <- list.files(file.path(output_dir, paste0("shards_", target)),
                      pattern = "\\.csv$", full.names = TRUE)
  if (!length(files)) stop("No shard files were found.", call. = FALSE)
  details <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  details <- details[!duplicated(details[c("scenario_id", "replicate")]), , drop = FALSE]
  summaries <- list()
  index <- 0L
  for (i in seq_len(nrow(scenarios))) {
    scenario <- scenarios[i, ]
    dat <- details[details$scenario_id == scenario$scenario_id, , drop = FALSE]
    if (!allow_partial && !is.null(n_outer) && nrow(dat) != n_outer) {
      stop("Scenario ", scenario$scenario_id, " has ", nrow(dat),
           " rows; expected ", n_outer, ".", call. = FALSE)
    }
    for (analysis in c("unadjusted", "adjusted")) {
      index <- index + 1L
      siga <- dat[[paste0("reject_siga_", analysis)]]
      rt <- dat[[paste0("reject_rt_", analysis)]]
      s_interval <- wilson_interval(sum(siga), length(siga))
      rt_interval <- wilson_interval(sum(rt), length(rt))
      difference <- paired_difference_interval(siga, rt)
      ratio_columns <- c(paste0("ratio_lower_", analysis), paste0("ratio_upper_", analysis))
      ratio_values <- unlist(dat[ratio_columns], use.names = FALSE)
      ratio_values <- ratio_values[is.finite(ratio_values)]
      summaries[[index]] <- data.frame(
        scenario_id = scenario$scenario_id,
        outcome = scenario$outcome,
        factor_count = scenario$factor_count,
        n_per_group = scenario$n_per_group,
        objective = scenario$objective,
        measure = scenario$measure,
        effect = scenario$effect,
        analysis = analysis,
        siga_percent = 100 * s_interval["estimate"],
        siga_lower_95_percent = 100 * s_interval["lower"],
        siga_upper_95_percent = 100 * s_interval["upper"],
        rt_percent = 100 * rt_interval["estimate"],
        rt_lower_95_percent = 100 * rt_interval["lower"],
        rt_upper_95_percent = 100 * rt_interval["upper"],
        difference_pp = 100 * difference["estimate"],
        difference_lower_95_pp = 100 * difference["lower"],
        difference_upper_95_pp = 100 * difference["upper"],
        mean_variance_ratio = if (length(ratio_values)) mean(ratio_values) else NA_real_,
        safeguard_frequency = mean(dat[[paste0("safeguard_", analysis)]], na.rm = TRUE),
        n_replicates = nrow(dat),
        stringsAsFactors = FALSE
      )
    }
  }
  summary <- do.call(rbind, summaries)
  safe_write_csv(details, file.path(output_dir, paste0("full_grid_details_", target, ".csv")))
  safe_write_csv(summary, file.path(output_dir, paste0("full_grid_summary_", target, ".csv")))
  summary
}

# The exact targeted pair-path supplemental workflow is delegated to
# production/pair_path_supplemental by scripts/03_run_pair_path_supplemental.R.

run_swift_direct_inspired <- function(n_outer = 100000L,
                                      B_randomization = 4999L,
                                      B_calibration = 100000L,
                                      base_seed = 20260724L,
                                      output_dir = "output/swift_direct_inspired",
                                      n_shards = 1L,
                                      shard_id = 1L,
                                      batch_size = 10L,
                                      calibration_batch = 1000L,
                                      progress = TRUE) {
  output_dir <- ensure_directory(output_dir)
  design <- make_swift_direct_design()
  model <- make_swift_direct_model(design)
  calibration_file <- file.path(output_dir, paste0("calibration_B0", B_calibration, ".rds"))
  calibration <- if (file.exists(calibration_file)) readRDS(calibration_file) else {
    object <- calibrate_three_path_design(
      B_calibration, design$total_n, design$patterns, design$profile_prob,
      design$pbc, design$weights, seed_value(base_seed, 1L, 0L, 90L),
      calibration_batch, progress
    )
    saveRDS(object, calibration_file)
    object
  }
  result_path <- file.path(output_dir, sprintf(
    "swift_shard_%03d_of_%03d.csv", shard_id, n_shards
  ))
  assigned <- seq.int(shard_id, n_outer, by = n_shards)
  completed <- integer()
  if (file.exists(result_path) && file.info(result_path)$size > 0) {
    existing <- tryCatch(read.csv(result_path), error = function(e) NULL)
    if (!is.null(existing)) completed <- unique(existing$replicate)
  }
  pending <- setdiff(assigned, completed)
  if (!length(pending)) return(invisible(output_dir))
  for (from in seq.int(1L, length(pending), by = batch_size)) {
    ids <- pending[from:min(length(pending), from + batch_size - 1L)]
    batch <- do.call(rbind, lapply(ids, function(replicate_id) {
      trial <- generate_swift_trial(
        design, model,
        seed_value(base_seed, 1L, replicate_id, 1L),
        seed_value(base_seed, 1L, replicate_id, 2L),
        seed_value(base_seed, 1L, replicate_id, 3L)
      )
      boundary <- design$boundary
      d <- model$d_at(boundary)
      scores <- make_boundary_scores(trial$y, trial$A, trial$X, boundary)
      analyse <- function(score) {
        statistic <- 0.5 * sum(trial$z * score)
        vs <- siga_sampling_variance(score, trial$X, calibration)
        vr <- siga_randomization_variance(vs, d, calibration, 1 / design$total_n)
        c(
          p_s = gaussian_pvalue(statistic, vs$variance, "greater")$p,
          p_r = gaussian_pvalue(statistic, vr$variance, "greater")$p,
          rho = vr$variance_ratio,
          safeguard = vr$safeguard_active
        )
      }
      unadjusted <- analyse(scores$unadjusted)
      adjusted <- analyse(scores$adjusted)
      rt <- fixed_score_randomization_test(
        trial$X, cbind(scores$unadjusted, scores$adjusted), trial$z,
        B_randomization, design$pbc, design$weights,
        seed_value(base_seed, 1L, replicate_id, 4L)
      )
      data.frame(
        replicate = replicate_id,
        p_siga_s_unadjusted = unadjusted["p_s"],
        p_siga_r_unadjusted = unadjusted["p_r"],
        p_rt_unadjusted = rt$greater[1L],
        p_siga_s_adjusted = adjusted["p_s"],
        p_siga_r_adjusted = adjusted["p_r"],
        p_rt_adjusted = rt$greater[2L],
        rho_unadjusted = unadjusted["rho"],
        rho_adjusted = adjusted["rho"],
        safeguard_unadjusted = as.logical(unadjusted["safeguard"]),
        safeguard_adjusted = as.logical(adjusted["safeguard"]),
        stringsAsFactors = FALSE
      )
    }))
    append_csv(batch, result_path)
    if (progress) message("SWIFT: wrote through replicate ", max(ids))
  }
  invisible(output_dir)
}

summarize_swift_direct <- function(output_dir, alpha = 0.05) {
  files <- list.files(output_dir, pattern = "swift_shard_.*\\.csv$", full.names = TRUE)
  details <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  details <- details[!duplicated(details$replicate), , drop = FALSE]
  rows <- lapply(c("unadjusted", "adjusted"), function(analysis) {
    s <- details[[paste0("p_siga_s_", analysis)]] <= alpha
    r <- details[[paste0("p_siga_r_", analysis)]] <= alpha
    rt <- details[[paste0("p_rt_", analysis)]] <= alpha
    data.frame(
      analysis = analysis,
      siga_s_percent = 100 * mean(s),
      siga_r_percent = 100 * mean(r),
      rt_percent = 100 * mean(rt),
      siga_r_minus_rt_pp = 100 * mean(r - rt),
      mean_rho = mean(details[[paste0("rho_", analysis)]], na.rm = TRUE),
      safeguard_frequency = mean(details[[paste0("safeguard_", analysis)]], na.rm = TRUE),
      n_replicates = nrow(details),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  safe_write_csv(details, file.path(output_dir, "swift_direct_details.csv"))
  safe_write_csv(summary, file.path(output_dir, "swift_direct_summary.csv"))
  summary
}
