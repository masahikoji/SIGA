#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || length(args) > 4L) {
  stop(
    paste0(
      "Usage: Rscript create_pair_path_publication_tables.R ",
      "SUMMARY_CSV SCENARIO_DEFINITIONS_CSV DESIGN_DIAGNOSTICS_CSV [OUTPUT_DIR]"
    ),
    call. = FALSE
  )
}

summary_path <- normalizePath(args[1L], mustWork = TRUE)
scenario_path <- normalizePath(args[2L], mustWork = TRUE)
diagnostic_path <- normalizePath(args[3L], mustWork = TRUE)
output_dir <- if (length(args) == 4L) args[4L] else dirname(summary_path)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## Main-manuscript table: selected scenarios and the adjusted analysis only.
## The complete 40-row table is written separately for the Supplement.
main_scenarios <- c(1L, 8L, 9L, 15L, 16L, 19L, 20L)
main_analysis <- "adjusted"

summary <- read.csv(summary_path, stringsAsFactors = FALSE, check.names = FALSE)
scenarios <- read.csv(scenario_path, stringsAsFactors = FALSE, check.names = FALSE)
diagnostics <- read.csv(diagnostic_path, stringsAsFactors = FALSE, check.names = FALSE)

required_summary <- c(
  "scenario_id", "set", "design_id", "outcome_type", "model_type",
  "direction", "analysis", "mean_variance_ratio_r_over_s",
  "predicted_rt_size_from_mean_ratio", "siga_s", "siga_r", "rt",
  "difference_siga_s_minus_rt", "difference_siga_r_minus_rt"
)
required_scenarios <- c(
  "scenario_id", "set", "design_id", "outcome_type", "model_type",
  "direction", "boundary", "pair_ratio"
)
required_diagnostics <- c(
  "design_id", "factor_count", "total_n", "pbc", "profile_type",
  "calibration_paths", "calibration_seconds", "min_generalized_ratio",
  "max_generalized_ratio", "max_abs_cross_D01", "max_abs_cross_D02",
  "max_abs_cross_D12", "max_abs_cross_G0102",
  "frobenius_cross_G0102", "max_abs_factor_correlation"
)

check_columns <- function(dat, required, label) {
  missing <- setdiff(required, names(dat))
  if (length(missing)) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}
check_columns(summary, required_summary, "Summary CSV")
check_columns(scenarios, required_scenarios, "Scenario-definition CSV")
check_columns(diagnostics, required_diagnostics, "Design-diagnostic CSV")

latex_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "--"
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x
}

clean_label <- function(x) {
  x <- gsub("_", "-", as.character(x), fixed = TRUE)
  latex_escape(x)
}

fmt <- function(x, digits = 3L) {
  x <- as.numeric(x)
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "--")
}

pct <- function(x, digits = 2L) {
  x <- as.numeric(x)
  ifelse(is.finite(x), formatC(100 * x, format = "f", digits = digits), "--")
}

sci <- function(x, digits = 2L) {
  x <- as.numeric(x)
  ifelse(is.finite(x), formatC(x, format = "e", digits = digits), "--")
}

design_code <- function(x) paste0("D", as.integer(x))
analysis_rank <- match(summary$analysis, c("adjusted", "unadjusted"))
summary <- summary[order(summary$scenario_id, analysis_rank), , drop = FALSE]
scenarios <- scenarios[order(scenarios$scenario_id), , drop = FALSE]
diagnostics <- diagnostics[order(diagnostics$design_id), , drop = FALSE]

## -------------------------------------------------------------------------
## Main-manuscript selected table
## -------------------------------------------------------------------------
main <- summary[
  summary$scenario_id %in% main_scenarios & summary$analysis == main_analysis,
  , drop = FALSE
]
main <- main[match(main_scenarios, main$scenario_id, nomatch = 0L), , drop = FALSE]

main_lines <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  paste0(
    "\\caption{Selected pair-path validation results for the adjusted analysis. ",
    "Rejection probabilities and differences are percentages; complete results are reported in the Supplement.}"
  ),
  "\\label{tab:pair-path-validation-main}",
  "\\small",
  "\\setlength{\\tabcolsep}{3.4pt}",
  "\\begin{tabular}{rcclrrrrrr}",
  "\\toprule",
  paste0(
    "ID & Design & Outcome & Setting & $V_{\\mathrm R}/V_{\\mathrm S}$ & ",
    "Pred. RT & SIGA-S & SIGA-R & RT & SIGA-R$-$RT \\\\"
  ),
  "\\midrule"
)
for (i in seq_len(nrow(main))) {
  row <- main[i, ]
  setting <- paste(clean_label(row$model_type), clean_label(row$direction), sep = "; ")
  main_lines <- c(
    main_lines,
    paste0(
      row$scenario_id, " & ",
      design_code(row$design_id), " & ",
      clean_label(row$outcome_type), " & ",
      setting, " & ",
      fmt(row$mean_variance_ratio_r_over_s, 3L), " & ",
      pct(row$predicted_rt_size_from_mean_ratio, 2L), " & ",
      pct(row$siga_s, 2L), " & ",
      pct(row$siga_r, 2L), " & ",
      pct(row$rt, 2L), " & ",
      pct(row$difference_siga_r_minus_rt, 2L), " \\\\"
    )
  )
}
main_lines <- c(main_lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(main_lines, file.path(output_dir, "pair_path_main_selected_table.tex"))

## -------------------------------------------------------------------------
## Supplement: complete 40-row longtable
## -------------------------------------------------------------------------
full_lines <- c(
  "\\begin{landscape}",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{2.6pt}",
  "\\renewcommand{\\arraystretch}{1.06}",
  "\\begin{longtable}{rcllllrrrrrrr}",
  paste0(
    "\\caption{Complete pair-path theory-validation results. Rejection probabilities, ",
    "predictions and differences are percentages.}",
    "\\label{tab:supp-pair-path-validation}\\\\"
  ),
  "\\toprule",
  paste0(
    "ID & Design & Outcome & Model & Direction & Analysis & ",
    "$V_{\\mathrm R}/V_{\\mathrm S}$ & Pred. RT & SIGA-S & SIGA-R & RT & ",
    "SIGA-S$-$RT & SIGA-R$-$RT \\\\"
  ),
  "\\midrule",
  "\\endfirsthead",
  "",
  "\\multicolumn{13}{c}{\\tablename\\ \\thetable{} -- continued}\\\\",
  "\\toprule",
  paste0(
    "ID & Design & Outcome & Model & Direction & Analysis & ",
    "$V_{\\mathrm R}/V_{\\mathrm S}$ & Pred. RT & SIGA-S & SIGA-R & RT & ",
    "SIGA-S$-$RT & SIGA-R$-$RT \\\\"
  ),
  "\\midrule",
  "\\endhead",
  "",
  "\\midrule",
  "\\multicolumn{13}{r}{Continued on next page}\\\\",
  "\\endfoot",
  "",
  "\\bottomrule",
  "\\endlastfoot"
)
for (i in seq_len(nrow(summary))) {
  row <- summary[i, ]
  full_lines <- c(
    full_lines,
    paste0(
      row$scenario_id, " & ",
      design_code(row$design_id), " & ",
      clean_label(row$outcome_type), " & ",
      clean_label(row$model_type), " & ",
      clean_label(row$direction), " & ",
      clean_label(row$analysis), " & ",
      fmt(row$mean_variance_ratio_r_over_s, 3L), " & ",
      pct(row$predicted_rt_size_from_mean_ratio, 2L), " & ",
      pct(row$siga_s, 2L), " & ",
      pct(row$siga_r, 2L), " & ",
      pct(row$rt, 2L), " & ",
      pct(row$difference_siga_s_minus_rt, 2L), " & ",
      pct(row$difference_siga_r_minus_rt, 2L), " \\\\"
    )
  )
}
full_lines <- c(full_lines, "\\end{longtable}", "\\end{landscape}")
writeLines(full_lines, file.path(output_dir, "pair_path_supplement_full_table.tex"))

## -------------------------------------------------------------------------
## Supplement: scenario definitions
## -------------------------------------------------------------------------
scenario_lines <- c(
  "\\begin{landscape}",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{3.0pt}",
  "\\renewcommand{\\arraystretch}{1.06}",
  "\\begin{longtable}{rclllllrr}",
  paste0(
    "\\caption{Definitions of the pair-path validation scenarios.}",
    "\\label{tab:supp-pair-path-scenarios}\\\\"
  ),
  "\\toprule",
  paste0(
    "ID & Set & Design & Outcome & Model & Direction & Boundary & ",
    "Design ratio & $\\alpha$ \\\\"
  ),
  "\\midrule",
  "\\endfirsthead",
  "",
  "\\multicolumn{9}{c}{\\tablename\\ \\thetable{} -- continued}\\\\",
  "\\toprule",
  paste0(
    "ID & Set & Design & Outcome & Model & Direction & Boundary & ",
    "Design ratio & $\\alpha$ \\\\"
  ),
  "\\midrule",
  "\\endhead",
  "",
  "\\bottomrule",
  "\\endlastfoot"
)
for (i in seq_len(nrow(scenarios))) {
  row <- scenarios[i, ]
  alpha_value <- if ("alpha" %in% names(scenarios)) row$alpha else NA_real_
  scenario_lines <- c(
    scenario_lines,
    paste0(
      row$scenario_id, " & ",
      clean_label(row$set), " & ",
      design_code(row$design_id), " & ",
      clean_label(row$outcome_type), " & ",
      clean_label(row$model_type), " & ",
      clean_label(row$direction), " & ",
      fmt(row$boundary, 2L), " & ",
      fmt(row$pair_ratio, 3L), " & ",
      fmt(alpha_value, 3L), " \\\\"
    )
  )
}
scenario_lines <- c(scenario_lines, "\\end{longtable}", "\\end{landscape}")
writeLines(scenario_lines, file.path(output_dir, "pair_path_supplement_scenario_table.tex"))

## -------------------------------------------------------------------------
## Supplement: allocation calibration diagnostics
## -------------------------------------------------------------------------
diagnostics$max_abs_cross_D <- apply(
  diagnostics[, c("max_abs_cross_D01", "max_abs_cross_D02", "max_abs_cross_D12")],
  1L, max, na.rm = TRUE
)

diagnostic_lines <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  paste0(
    "\\caption{Three-copy allocation-calibration diagnostics. ",
    "The cross-copy quantities should be close to zero.}"
  ),
  "\\label{tab:supp-pair-path-diagnostics}",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{3.0pt}",
  "\\begin{tabular}{rlllrrrrrr}",
  "\\toprule",
  paste0(
    "Design & $(K,n,p_{\\rm bc})$ & Profile & $B_0$ & Min ratio & Max ratio & ",
    "Max cross-$D$ & Max cross-$G$ & Frobenius cross-$G$ & Max factor corr. \\\\"
  ),
  "\\midrule"
)
for (i in seq_len(nrow(diagnostics))) {
  row <- diagnostics[i, ]
  design_tuple <- paste0(
    "$(", row$factor_count, ",", row$total_n, ",", fmt(row$pbc, 2L), ")$"
  )
  diagnostic_lines <- c(
    diagnostic_lines,
    paste0(
      design_code(row$design_id), " & ",
      design_tuple, " & ",
      clean_label(row$profile_type), " & ",
      formatC(as.integer(row$calibration_paths), format = "d"), " & ",
      fmt(row$min_generalized_ratio, 3L), " & ",
      fmt(row$max_generalized_ratio, 3L), " & ",
      sci(row$max_abs_cross_D, 2L), " & ",
      sci(row$max_abs_cross_G0102, 2L), " & ",
      sci(row$frobenius_cross_G0102, 2L), " & ",
      fmt(row$max_abs_factor_correlation, 3L), " \\\\"
    )
  )
}
diagnostic_lines <- c(
  diagnostic_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(diagnostic_lines, file.path(output_dir, "pair_path_supplement_diagnostic_table.tex"))

outputs <- c(
  "pair_path_main_selected_table.tex",
  "pair_path_supplement_full_table.tex",
  "pair_path_supplement_scenario_table.tex",
  "pair_path_supplement_diagnostic_table.tex"
)
message("Publication tables written to: ", normalizePath(output_dir, mustWork = FALSE))
for (x in outputs) message("  - ", x)
