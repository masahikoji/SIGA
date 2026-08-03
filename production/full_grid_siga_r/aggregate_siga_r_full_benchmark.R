#!/usr/bin/env Rscript

# =============================================================================
# Aggregate the comprehensive SIGA-R benchmark and create publication outputs
# =============================================================================
#
# Required input
# --------------
# The directory created by simulation_siga_r_full_benchmark.R, containing:
#   configuration/full_scenario_grid.csv
#   scenario_shards/scenario_*/shard_*.csv
#
# Default manuscript requirement
# ------------------------------
# Exactly 100,000 unique outer replicates for each of all 56 scenarios.
# Set PWRT_ALLOW_PARTIAL=1 only for smoke-test aggregation.
# =============================================================================

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
OUTPUT_ROOT <- path.expand(Sys.getenv(
  "PWRT_OUTPUT_DIR",
  unset = file.path(BASE_DIR, "siga_r_full_benchmark_output")
))
if (!dir.exists(OUTPUT_ROOT)) stop("Output base directory not found: ", OUTPUT_ROOT, call. = FALSE)
OUTPUT_ROOT <- normalizePath(OUTPUT_ROOT, winslash = "/", mustWork = TRUE)

EXPECTED_OUTER <- suppressWarnings(as.integer(Sys.getenv("PWRT_N_OUTER", unset = "100000")))
EXPECTED_RERAND <- suppressWarnings(as.integer(Sys.getenv("PWRT_N_RERAND", unset = "4999")))
EXPECTED_CALIBRATION <- suppressWarnings(as.integer(Sys.getenv("PWRT_N_CALIBRATION", unset = "100000")))
EXPECTED_SEED <- suppressWarnings(as.integer(Sys.getenv("PWRT_SEED", unset = "20260801")))
EXPECTED_PBC <- suppressWarnings(as.numeric(Sys.getenv("PWRT_P_BIASED_COIN", unset = "0.80")))
EXPECTED_EPSILON_EXPONENT <- suppressWarnings(as.numeric(Sys.getenv("PWRT_EPSILON_EXPONENT", unset = "1.0")))
if (is.na(EXPECTED_OUTER) || EXPECTED_OUTER < 1L) stop("PWRT_N_OUTER is invalid.")
if (is.na(EXPECTED_RERAND) || EXPECTED_RERAND < 1L) stop("PWRT_N_RERAND is invalid.")
if (is.na(EXPECTED_CALIBRATION) || EXPECTED_CALIBRATION < 2L) stop("PWRT_N_CALIBRATION is invalid.")
ALLOW_PARTIAL <- identical(Sys.getenv("PWRT_ALLOW_PARTIAL", unset = "0"), "1")

RUN_TAG <- paste0(
  "M", EXPECTED_OUTER,
  "_B", EXPECTED_RERAND,
  "_Bpsi", EXPECTED_CALIBRATION,
  "_pbc", gsub("[.]", "p", formatC(EXPECTED_PBC, format = "f", digits = 3)),
  "_eps", gsub("[.]", "p", formatC(EXPECTED_EPSILON_EXPONENT, format = "f", digits = 2)),
  "_seed", EXPECTED_SEED
)
RUN_DIR <- file.path(OUTPUT_ROOT, RUN_TAG)
if (!dir.exists(RUN_DIR)) stop("Run directory not found: ", RUN_DIR, call. = FALSE)
RUN_DIR <- normalizePath(RUN_DIR, winslash = "/", mustWork = TRUE)

CONFIG_DIR <- file.path(RUN_DIR, "configuration")
SCENARIO_DIR <- file.path(RUN_DIR, "scenario_shards")
PUBLICATION_DIR <- file.path(RUN_DIR, "publication_tables")
dir.create(PUBLICATION_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(PUBLICATION_DIR)) stop("Could not create publication directory.")

scenario_grid_path <- file.path(CONFIG_DIR, "full_scenario_grid.csv")
if (!file.exists(scenario_grid_path)) {
  stop("Scenario grid not found: ", scenario_grid_path, call. = FALSE)
}
SCENARIOS <- read.csv(scenario_grid_path, stringsAsFactors = FALSE)
if (!ALLOW_PARTIAL && nrow(SCENARIOS) != 56L) {
  stop("The manuscript grid must contain 56 scenarios; found ", nrow(SCENARIOS), ".")
}

wilson_interval <- function(x, n, conf.level = 0.95) {
  z <- qnorm(1 - (1 - conf.level) / 2)
  p <- x / n
  den <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(estimate = p, lower = max(0, center - half), upper = min(1, center + half))
}

safe_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

safe_quantile <- function(x, prob) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  unname(quantile(x, probs = prob, names = FALSE, type = 8))
}

read_scenario_detail <- function(scenario) {
  prefix <- sprintf("scenario_%02d_", scenario$scenario_id)
  dirs <- list.dirs(SCENARIO_DIR, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[startsWith(basename(dirs), prefix)]
  if (length(dirs) != 1L) {
    stop(
      "Expected one directory for scenario ", scenario$scenario_id,
      "; found ", length(dirs), call. = FALSE
    )
  }
  files <- list.files(
    dirs,
    pattern = "^shard_[0-9]+_of_[0-9]+[.]csv$",
    full.names = TRUE
  )
  if (!length(files)) stop("No shard CSV files for scenario ", scenario$scenario_id)
  detail <- do.call(rbind, lapply(files, function(path) {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }))
  detail <- detail[detail$scenario_id == scenario$scenario_id, , drop = FALSE]
  detail <- detail[!duplicated(detail$replicate), , drop = FALSE]
  detail <- detail[order(detail$replicate), , drop = FALSE]

  required_configuration <- c(
    "configured_outer_trials", "rerandomizations_per_trial",
    "pair_calibration_replicates", "epsilon_exponent"
  )
  missing_configuration <- setdiff(required_configuration, names(detail))
  if (length(missing_configuration)) {
    stop(
      "Scenario ", scenario$scenario_id,
      " is missing configuration columns: ",
      paste(missing_configuration, collapse = ", "), call. = FALSE
    )
  }
  unique_outer <- unique(as.integer(detail$configured_outer_trials))
  unique_rerand <- unique(as.integer(detail$rerandomizations_per_trial))
  unique_calibration <- unique(as.integer(detail$pair_calibration_replicates))
  unique_epsilon <- unique(round(as.numeric(detail$epsilon_exponent), 12))
  if (length(unique_outer) != 1L || unique_outer != EXPECTED_OUTER ||
      length(unique_rerand) != 1L || unique_rerand != EXPECTED_RERAND ||
      length(unique_calibration) != 1L || unique_calibration != EXPECTED_CALIBRATION ||
      length(unique_epsilon) != 1L || abs(unique_epsilon - EXPECTED_EPSILON_EXPONENT) > 1e-12) {
    stop(
      "Scenario ", scenario$scenario_id,
      " contains simulation settings inconsistent with the requested run tag.",
      call. = FALSE
    )
  }

  if (!ALLOW_PARTIAL) {
    if (nrow(detail) != EXPECTED_OUTER) {
      stop(
        "Scenario ", scenario$scenario_id, " has ", nrow(detail),
        " unique replicates; expected ", EXPECTED_OUTER, call. = FALSE
      )
    }
    if (!identical(as.integer(detail$replicate), seq_len(EXPECTED_OUTER))) {
      missing <- setdiff(seq_len(EXPECTED_OUTER), detail$replicate)
      stop(
        "Scenario ", scenario$scenario_id,
        " has nonconsecutive replicate IDs. First missing: ",
        paste(head(missing, 10L), collapse = ","), call. = FALSE
      )
    }
  }
  detail
}

summarize_analysis <- function(detail, scenario, analysis) {
  suffix <- if (analysis == "Unadjusted") "unadjusted" else "adjusted"
  col_s <- paste0("reject_s_", suffix)
  col_r <- paste0("reject_r_", suffix)
  col_rt <- paste0("reject_rt_", suffix)
  p_s <- paste0("p_s_", suffix)
  p_r <- paste0("p_r_", suffix)
  p_rt <- paste0("p_rt_", suffix)

  required <- c(col_s, col_r, col_rt, p_s, p_r, p_rt)
  missing <- setdiff(required, names(detail))
  if (length(missing)) {
    stop(
      "Scenario ", scenario$scenario_id, " is missing columns: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }

  rs <- as.integer(detail[[col_s]])
  rr <- as.integer(detail[[col_r]])
  rt <- as.integer(detail[[col_rt]])
  n <- length(rt)
  ci_s <- wilson_interval(sum(rs), n)
  ci_r <- wilson_interval(sum(rr), n)
  ci_rt <- wilson_interval(sum(rt), n)

  paired_diff <- rr - rt
  diff_est <- mean(paired_diff)
  diff_se <- if (n > 1L) sd(paired_diff) / sqrt(n) else NA_real_
  diff_lo <- diff_est - qnorm(0.975) * diff_se
  diff_hi <- diff_est + qnorm(0.975) * diff_se

  ratio_rs_1 <- detail[[paste0("ratio_r_over_s_", suffix, "_1")]]
  ratio_rs_2 <- detail[[paste0("ratio_r_over_s_", suffix, "_2")]]
  ratio_rrt_1 <- detail[[paste0("ratio_r_over_rt_", suffix, "_1")]]
  ratio_rrt_2 <- detail[[paste0("ratio_r_over_rt_", suffix, "_2")]]
  safeguard_1 <- detail[[paste0("safeguard_", suffix, "_1")]]
  safeguard_2 <- detail[[paste0("safeguard_", suffix, "_2")]]

  data.frame(
    scenario_id = scenario$scenario_id,
    scenario_label = scenario$scenario_label,
    outcome = scenario$outcome,
    analysis = analysis,
    design_id = scenario$design_id,
    factor_count = scenario$factor_count,
    target_per_group = scenario$target_per_group,
    total_n = scenario$total_n,
    objective = scenario$objective,
    scenario_role = scenario$scenario_role,
    true_effect = scenario$true_effect,
    alpha = scenario$alpha,
    completed_outer_trials = n,
    siga_s_rejection = unname(ci_s["estimate"]),
    siga_s_lower_95 = unname(ci_s["lower"]),
    siga_s_upper_95 = unname(ci_s["upper"]),
    siga_r_rejection = unname(ci_r["estimate"]),
    siga_r_lower_95 = unname(ci_r["lower"]),
    siga_r_upper_95 = unname(ci_r["upper"]),
    rt_rejection = unname(ci_rt["estimate"]),
    rt_lower_95 = unname(ci_rt["lower"]),
    rt_upper_95 = unname(ci_rt["upper"]),
    siga_r_minus_rt = diff_est,
    siga_r_minus_rt_se = diff_se,
    siga_r_minus_rt_lower_95 = diff_lo,
    siga_r_minus_rt_upper_95 = diff_hi,
    siga_s_minus_rt = mean(rs - rt),
    pvalue_mae_siga_r_vs_rt = mean(abs(detail[[p_r]] - detail[[p_rt]])),
    pvalue_rmse_siga_r_vs_rt = sqrt(mean((detail[[p_r]] - detail[[p_rt]])^2)),
    pvalue_correlation_siga_r_vs_rt = cor(detail[[p_r]], detail[[p_rt]]),
    mean_ratio_r_over_s_1 = safe_mean(ratio_rs_1),
    q025_ratio_r_over_s_1 = safe_quantile(ratio_rs_1, 0.025),
    q975_ratio_r_over_s_1 = safe_quantile(ratio_rs_1, 0.975),
    mean_ratio_r_over_s_2 = safe_mean(ratio_rs_2),
    q025_ratio_r_over_s_2 = safe_quantile(ratio_rs_2, 0.025),
    q975_ratio_r_over_s_2 = safe_quantile(ratio_rs_2, 0.975),
    mean_ratio_r_over_rt_1 = safe_mean(ratio_rrt_1),
    q025_ratio_r_over_rt_1 = safe_quantile(ratio_rrt_1, 0.025),
    q975_ratio_r_over_rt_1 = safe_quantile(ratio_rrt_1, 0.975),
    mean_ratio_r_over_rt_2 = safe_mean(ratio_rrt_2),
    q025_ratio_r_over_rt_2 = safe_quantile(ratio_rrt_2, 0.025),
    q975_ratio_r_over_rt_2 = safe_quantile(ratio_rrt_2, 0.975),
    safeguard_rate_1 = safe_mean(as.numeric(safeguard_1)),
    safeguard_rate_2 = safe_mean(as.numeric(safeguard_2)),
    stringsAsFactors = FALSE
  )
}

summary_rows <- vector("list", 2L * nrow(SCENARIOS))
pos <- 0L
for (i in seq_len(nrow(SCENARIOS))) {
  scenario <- SCENARIOS[i, ]
  message("Aggregating scenario ", scenario$scenario_id, "/", nrow(SCENARIOS))
  detail <- read_scenario_detail(scenario)
  for (analysis in c("Unadjusted", "Adjusted")) {
    pos <- pos + 1L
    summary_rows[[pos]] <- summarize_analysis(detail, scenario, analysis)
  }
  rm(detail)
  invisible(gc(FALSE))
}

SUMMARY <- do.call(rbind, summary_rows)
SUMMARY <- SUMMARY[order(SUMMARY$scenario_id, SUMMARY$analysis), , drop = FALSE]
summary_csv <- file.path(PUBLICATION_DIR, "siga_r_full_benchmark_summary.csv")
write.csv(SUMMARY, summary_csv, row.names = FALSE, na = "")

# -----------------------------------------------------------------------------
# Publication formatting
# -----------------------------------------------------------------------------

objective_label <- function(x) {
  switch(x,
    superiority = "Superiority",
    noninferiority = "Non-inferiority",
    equivalence = "Equivalence",
    x
  )
}

measure_label <- function(x) {
  switch(x,
    type1 = "Type I error",
    power = "Power",
    type1_lower = "Error at lower limit",
    type1_upper = "Error at upper limit",
    x
  )
}

format_effect <- function(x) sprintf("%+.3f", x)

format_ratio_pair <- function(row) {
  a <- row$mean_ratio_r_over_s_1
  b <- row$mean_ratio_r_over_s_2
  if (is.finite(b)) sprintf("%.3f/%.3f", a, b) else sprintf("%.3f", a)
}

write_longtable <- function(data, path, caption, label) {
  lines <- c(
    "\\begin{longtable}{rrllrrrrr}",
    paste0("\\caption{", caption, "}\\label{", label, "}\\\\"),
    "\\toprule",
    "Factors & $n$/group & Objective & Measure & Effect & SIGA-R (\\%) & RT (\\%) & Difference (pp) & $V_{\\rm R}/V_{\\rm S}$ \\\\",
    "\\midrule",
    "\\endfirsthead",
    "\\multicolumn{9}{c}{\\tablename\\ \\thetable{} continued}\\\\",
    "\\toprule",
    "Factors & $n$/group & Objective & Measure & Effect & SIGA-R (\\%) & RT (\\%) & Difference (pp) & $V_{\\rm R}/V_{\\rm S}$ \\\\",
    "\\midrule",
    "\\endhead",
    "\\midrule",
    "\\multicolumn{9}{r}{Continued on next page}\\\\",
    "\\endfoot",
    "\\bottomrule",
    "\\endlastfoot"
  )
  for (i in seq_len(nrow(data))) {
    row <- data[i, ]
    lines <- c(lines, paste0(
      row$factor_count, " & ",
      row$target_per_group, " & ",
      objective_label(row$objective), " & ",
      measure_label(row$scenario_role), " & ",
      format_effect(row$true_effect), " & ",
      sprintf("%.2f", 100 * row$siga_r_rejection), " & ",
      sprintf("%.2f", 100 * row$rt_rejection), " & ",
      sprintf("%+.2f", 100 * row$siga_r_minus_rt), " & ",
      format_ratio_pair(row), " \\\\"
    ))
  }
  lines <- c(
    lines,
    "\\end{longtable}",
    "\\noindent\\footnotesize RT denotes the reference fixed-score randomization test. For equivalence scenarios, the two variance ratios correspond to the lower and upper boundaries. SIGA-R and RT were evaluated on the same outer trials. Differences should be calculated from unrounded estimates.\\normalsize",
    ""
  )
  writeLines(lines, path, useBytes = TRUE)
}

for (outcome in c("continuous", "binary")) {
  for (analysis in c("Unadjusted", "Adjusted")) {
    dat <- SUMMARY[SUMMARY$outcome == outcome & SUMMARY$analysis == analysis, , drop = FALSE]
    title_outcome <- if (outcome == "continuous") "continuous outcomes" else "binary outcomes"
    title_analysis <- tolower(analysis)
    file_tag <- paste(outcome, tolower(analysis), sep = "_")
    write_longtable(
      dat,
      file.path(PUBLICATION_DIR, paste0("supp_siga_r_", file_tag, ".tex")),
      paste0(
        "Randomization-targeted operating characteristics for ", title_outcome,
        " using the ", title_analysis, " analysis."
      ),
      paste0("tab:siga-r-", outcome, "-", tolower(analysis))
    )
  }
}

# Compact main-text summary, without presenting SIGA-S and SIGA-R as competitors.
compact_rows <- list()
pos <- 0L
for (outcome in c("continuous", "binary")) {
  for (analysis in c("Unadjusted", "Adjusted")) {
    dat <- SUMMARY[SUMMARY$outcome == outcome & SUMMARY$analysis == analysis, , drop = FALSE]
    power <- dat[dat$scenario_role == "power", , drop = FALSE]
    ratios <- c(dat$mean_ratio_r_over_s_1, dat$mean_ratio_r_over_s_2)
    ratios <- ratios[is.finite(ratios)]
    pos <- pos + 1L
    compact_rows[[pos]] <- data.frame(
      Outcome = if (outcome == "continuous") "Continuous" else "Binary",
      Analysis = analysis,
      Scenarios = nrow(dat),
      MaxAbsDifferenceAllPP = 100 * max(abs(dat$siga_r_minus_rt)),
      MaxAbsDifferencePowerPP = 100 * max(abs(power$siga_r_minus_rt)),
      MinMeanVarianceRatio = min(ratios),
      MaxMeanVarianceRatio = max(ratios),
      MaxSafeguardRate = max(c(dat$safeguard_rate_1, dat$safeguard_rate_2), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}
COMPACT <- do.call(rbind, compact_rows)
write.csv(
  COMPACT,
  file.path(PUBLICATION_DIR, "siga_r_full_benchmark_compact_summary.csv"),
  row.names = FALSE
)

main_lines <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Summary of the comprehensive validation of SIGA-R against the reference randomization test.}",
  "\\label{tab:siga-r-comprehensive-summary}",
  "\\begin{threeparttable}",
  "\\small",
  "\\setlength{\\tabcolsep}{5.5pt}",
  "\\renewcommand{\\arraystretch}{1.10}",
  "\\begin{tabular}{llrrrr}",
  "\\toprule",
  "Outcome & Analysis & Scenarios & Max. difference & Max. power difference & Mean $V_{\\rm R}/V_{\\rm S}$ range \\\\",
  " & & & (pp) & (pp) & \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(COMPACT))) {
  row <- COMPACT[i, ]
  main_lines <- c(main_lines, paste0(
    row$Outcome, " & ", row$Analysis, " & ", row$Scenarios, " & ",
    sprintf("%.2f", row$MaxAbsDifferenceAllPP), " & ",
    sprintf("%.2f", row$MaxAbsDifferencePowerPP), " & ",
    sprintf("%.3f--%.3f", row$MinMeanVarianceRatio, row$MaxMeanVarianceRatio),
    " \\\\"
  ))
}
main_lines <- c(
  main_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}[flushleft]",
  "\\footnotesize",
  "\\item Differences are absolute differences between SIGA-R and the reference fixed-score randomization test, calculated from unrounded rejection probabilities. The power comparison includes the superiority, non-inferiority, and equivalence power scenarios. SIGA-S is not included because this table evaluates fidelity to the conditional randomization target rather than ranks the two SIGA procedures as competing tests.",
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}",
  ""
)
writeLines(
  main_lines,
  file.path(PUBLICATION_DIR, "main_siga_r_comprehensive_summary.tex"),
  useBytes = TRUE
)

# Automatically generated results paragraph with exact values.
max_all <- max(abs(SUMMARY$siga_r_minus_rt))
max_power <- max(abs(SUMMARY$siga_r_minus_rt[SUMMARY$scenario_role == "power"]))
max_type1 <- max(abs(SUMMARY$siga_r_minus_rt[SUMMARY$scenario_role != "power"]))
all_ratios <- c(SUMMARY$mean_ratio_r_over_s_1, SUMMARY$mean_ratio_r_over_s_2)
all_ratios <- all_ratios[is.finite(all_ratios)]
max_safeguard <- max(c(SUMMARY$safeguard_rate_1, SUMMARY$safeguard_rate_2), na.rm = TRUE)

results_text <- paste0(
  "Across all 56 benchmark scenarios and both the unadjusted and adjusted analyses, ",
  "the maximum absolute difference between SIGA-R and the reference randomization test was ",
  sprintf("%.2f", 100 * max_all), " percentage points. The corresponding maxima were ",
  sprintf("%.2f", 100 * max_type1), " percentage points for null-boundary operating characteristics and ",
  sprintf("%.2f", 100 * max_power), " percentage points for power. The scenario-level mean variance ratio ",
  "$\\widehat V_{\\rm R}/\\widehat V_{\\rm S}$ ranged from ",
  sprintf("%.3f", min(all_ratios)), " to ", sprintf("%.3f", max(all_ratios)), ". ",
  "The numerical safeguard was activated in at most ",
  sprintf("%.4f", 100 * max_safeguard), "\\% of outer trials in any scenario-boundary-analysis combination."
)
writeLines(
  results_text,
  file.path(PUBLICATION_DIR, "siga_r_full_benchmark_results_paragraph.tex"),
  useBytes = TRUE
)

# Configuration audit file.
audit <- data.frame(
  item = c(
    "scenarios", "analysis_rows", "expected_outer_per_scenario",
    "minimum_completed_outer", "maximum_completed_outer",
    "maximum_absolute_siga_r_minus_rt_pp",
    "maximum_power_absolute_siga_r_minus_rt_pp",
    "maximum_safeguard_rate"
  ),
  value = c(
    length(unique(SUMMARY$scenario_id)), nrow(SUMMARY), EXPECTED_OUTER,
    min(SUMMARY$completed_outer_trials), max(SUMMARY$completed_outer_trials),
    100 * max_all, 100 * max_power, max_safeguard
  ),
  stringsAsFactors = FALSE
)
write.csv(audit, file.path(PUBLICATION_DIR, "aggregation_audit.csv"), row.names = FALSE)

message("Aggregation complete.")
message("Summary CSV: ", summary_csv)
message("Publication files: ", PUBLICATION_DIR)
