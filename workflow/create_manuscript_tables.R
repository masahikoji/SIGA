#!/usr/bin/env Rscript

# Create table-ready TeX fragments from manuscript_outputs. This script does
# not simulate data. It is run after build_manuscript_outputs.R.

options(stringsAsFactors = FALSE, warn = 1)

script_directory <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = FALSE)))
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
root_dir <- normalizePath(path.expand(Sys.getenv("SIGA_PROJECT_DIR", unset = file.path(script_directory(), ".."))), winslash = "/", mustWork = TRUE)
input_dir <- file.path(root_dir, "manuscript_outputs")
expected_dir <- file.path(root_dir, "expected")
out_dir <- file.path(input_dir, "tex_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) stop("Could not create: ", out_dir, call. = FALSE)

read_req <- function(path) {
  if (!file.exists(path)) stop("Required file not found: ", path, call. = FALSE)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
write_tex <- function(lines, name) writeLines(lines, file.path(out_dir, name), useBytes = TRUE)
escape_tex <- function(x) gsub("_", "\\\\_", as.character(x), fixed = TRUE)
fmt_effect <- function(x) sprintf("%+.3f", as.numeric(x))
fmt_pct <- function(x) sprintf("%.2f", 100 * as.numeric(x))
objective_label <- function(x) {
  z <- as.character(x)
  z[x == "superiority"] <- "Superiority"
  z[x == "noninferiority"] <- "Non-inferiority"
  z[x == "equivalence"] <- "Equivalence"
  z
}
role_label <- function(x) {
  z <- as.character(x)
  z[x == "type1"] <- "Type I error"
  z[x == "power"] <- "Power"
  z[x == "type1_lower"] <- "Error at $L$"
  z[x == "type1_upper"] <- "Error at $U$"
  z
}

continuous_parameters <- read_req(file.path(expected_dir, "reported_continuous_power_parameters.csv"))
binary_parameters <- read_req(file.path(expected_dir, "reported_binary_power_parameters.csv"))
operating <- read_req(file.path(input_dir, "combined_operating_characteristics.csv"))
timing <- read_req(file.path(input_dir, "timing_summary_complete.csv"))
calibration <- read_req(file.path(input_dir, "calibration_timing_repeats.csv"))
case <- read_req(file.path(input_dir, "swift_direct_case_study_summary.csv"))

# Web Table 1.
lines <- c(
  "\\begin{table}[!htbp]", "\\centering",
  "\\caption{Treatment effects used in the continuous-outcome power scenarios.}",
  "\\label{tab:power-effects}", "\\begin{tabular}{ccccc}", "\\toprule",
  "Factors & Size/group & Superiority $\\tau$ & Non-inferiority $\\tau$ & Equivalence $\\tau$ \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(continuous_parameters))) {
  x <- continuous_parameters[i, ]
  lines <- c(lines, sprintf("%d & %d & %.3f & %.3f & %.3f \\\\", x$factor_count, x$target_per_group,
                            x$superiority_effect, x$noninferiority_effect, x$equivalence_effect))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
write_tex(lines, "web_table_1_continuous_power_parameters.tex")

# Web Table 2.
lines <- c(
  "\\begin{table}[!htbp]", "\\centering",
  "\\caption{Marginal risk differences and equivalence limits used in the binary-outcome power scenarios.}",
  "\\label{tab:binary-effects}", "\\begin{tabular}{cccccc}", "\\toprule",
  "Factors & Size/group & Superiority $\\Delta$ & Non-inferiority $\\Delta$ & Equivalence limits & Equivalence $\\Delta$ \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(binary_parameters))) {
  x <- binary_parameters[i, ]
  limits <- sprintf("$(%.3f,%.3f)$", x$equivalence_lower, x$equivalence_upper)
  lines <- c(lines, sprintf("%d & %d & %.3f & %.3f & %s & %.3f \\\\", x$factor_count, x$target_per_group,
                            x$superiority_risk_difference, x$noninferiority_risk_difference,
                            limits, x$equivalence_risk_difference))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
write_tex(lines, "web_table_2_binary_power_parameters.tex")

objective_order <- c(superiority = 1L, noninferiority = 2L, equivalence = 3L)
role_order <- c(type1 = 1L, type1_lower = 1L, type1_upper = 2L, power = 3L)
outcome_order <- c(Continuous = 1L, Binary = 2L)
analysis_order <- c(Unadjusted = 1L, Adjusted = 2L)

write_operating <- function(k, filename) {
  dat <- operating[operating$factor_count == k, ]
  dat <- dat[order(outcome_order[dat$outcome], dat$target_per_group,
                   objective_order[dat$objective], role_order[dat$scenario_role]), ]
  word <- if (k == 2L) "two" else "five"
  lines <- c(
    "\\begin{landscape}", "\\begingroup", "\\singlespacing", "\\scriptsize",
    "\\setlength{\\tabcolsep}{3.2pt}", "\\renewcommand{\\arraystretch}{1.08}",
    "\\begin{longtable}{@{}lllllrrrr@{}}",
    sprintf("\\caption{Operating characteristics for %s-factor designs.}\\label{tab:operating-characteristics-%s-factor}\\\\", word, word),
    "\\toprule",
    "Outcome & Size & Objective & Measure & Effect & \\multicolumn{2}{c}{Unadjusted} & \\multicolumn{2}{c}{Adjusted} \\\\",
    "\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}",
    " & & & & & SIGA (\\%) & RT (\\%) & SIGA (\\%) & RT (\\%) \\\\",
    "\\midrule", "\\endfirsthead",
    "\\multicolumn{9}{c}{\\tablename\\ \\thetable{} -- continued}\\\\", "\\toprule",
    "Outcome & Size & Objective & Measure & Effect & \\multicolumn{2}{c}{Unadjusted} & \\multicolumn{2}{c}{Adjusted} \\\\",
    "\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}",
    " & & & & & SIGA (\\%) & RT (\\%) & SIGA (\\%) & RT (\\%) \\\\",
    "\\midrule", "\\endhead", "\\midrule", "\\multicolumn{9}{r}{Continued on next page}\\\\", "\\endfoot", "\\bottomrule", "\\endlastfoot"
  )
  previous <- ""
  for (i in seq_len(nrow(dat))) {
    group <- paste(dat$outcome[i], dat$target_per_group[i])
    if (nzchar(previous) && group != previous) lines <- c(lines, "\\addlinespace")
    previous <- group
    lines <- c(lines, sprintf("%s & %d & %s & %s & %s & %s & %s & %s & %s \\\\",
                              dat$outcome[i], dat$target_per_group[i], objective_label(dat$objective[i]),
                              role_label(dat$scenario_role[i]), fmt_effect(dat$true_effect[i]),
                              fmt_pct(dat$proposed_unadjusted[i]), fmt_pct(dat$reference_unadjusted[i]),
                              fmt_pct(dat$proposed_adjusted[i]), fmt_pct(dat$reference_adjusted[i])))
  }
  lines <- c(lines, "\\end{longtable}", "\\endgroup", "\\end{landscape}")
  write_tex(lines, filename)
}
write_operating(2L, "web_table_3_operating_characteristics_two_factors.tex")
write_operating(5L, "web_table_4_operating_characteristics_five_factors.tex")

write_timing <- function(k, filename) {
  dat <- timing[timing$factor_count == k, ]
  dat <- dat[order(outcome_order[dat$outcome], dat$target_per_group,
                   objective_order[dat$objective], analysis_order[dat$analysis]), ]
  word <- if (k == 2L) "two" else "five"
  lines <- c(
    "\\begin{table}[!htbp]", "\\centering",
    sprintf("\\caption{Projected single-core computation times for %s-factor designs.}", word),
    sprintf("\\label{tab:computation-time-%s-factor}", word),
    "\\small", "\\begin{tabular}{llrlrr}", "\\toprule",
    "Outcome & Analysis & $n$/group & Objective & SIGA projected time (min) & RT projected time (min) \\\\",
    "\\midrule"
  )
  previous <- ""
  for (i in seq_len(nrow(dat))) {
    group <- paste(dat$outcome[i], dat$target_per_group[i])
    if (nzchar(previous) && group != previous) lines <- c(lines, "\\addlinespace")
    previous <- group
    lines <- c(lines, sprintf("%s & %s & %d & %s & %.3f & %.2f \\\\", dat$outcome[i], dat$analysis[i],
                              dat$target_per_group[i], objective_label(dat$objective[i]),
                              dat$siga_projected_minutes_100k[i], dat$rt_projected_minutes_100k[i]))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  write_tex(lines, filename)
}
write_timing(2L, "web_table_5_timing_two_factors.tex")
write_timing(5L, "web_table_6_timing_five_factors.tex")

# Web Table 7.
cal_wide <- reshape(calibration[, c("factor_count", "target_per_group", "total_n", "allocation_calibration_paths", "calibration_repeat", "calibration_seconds")],
                    idvar = c("factor_count", "target_per_group", "total_n", "allocation_calibration_paths"),
                    timevar = "calibration_repeat", direction = "wide")
names(cal_wide) <- sub("calibration_seconds\\.", "run", names(cal_wide))
names(cal_wide)[names(cal_wide) == "allocation_calibration_paths"] <- "B0"
cal_wide$mean_seconds <- rowMeans(cal_wide[c("run1", "run2", "run3")])
cal_wide <- cal_wide[order(cal_wide$factor_count, cal_wide$target_per_group), ]
lines <- c(
  "\\begin{table}[!htbp]", "\\centering",
  "\\caption{Single-core computation times for one allocation-only SIGA calibration.}",
  "\\label{tab:calibration-time}", "\\small", "\\begin{tabular}{rrrrrrrr}", "\\toprule",
  "Factors & $n$/group & Total $n$ & $B_0$ & Run 1 (s) & Run 2 (s) & Run 3 (s) & Mean (s) \\\\", "\\midrule"
)
for (i in seq_len(nrow(cal_wide))) {
  x <- cal_wide[i, ]
  lines <- c(lines, sprintf("%d & %d & %d & %d & %.3f & %.3f & %.3f & %.3f \\\\",
                            x$factor_count, x$target_per_group, x$total_n, x$B0,
                            x$run1, x$run2, x$run3, x$mean_seconds))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
write_tex(lines, "web_table_7_calibration_timing.tex")

# Main manuscript illustration table.
x <- case[1L, ]
lines <- c(
  "\\begin{table}[!htbp]", "\\centering",
  "\\caption{Power and computation time in the SWIFT DIRECT-inspired illustration.}",
  "\\label{tab:swift-direct-siga-summary}", "\\begin{tabular}{lrr}", "\\toprule",
  " & SIGA & Randomization test \\\\", "\\midrule",
  "\\multicolumn{3}{l}{\\textit{Estimated non-inferiority power (\\%)}}\\\\",
  sprintf("Unadjusted score & %.2f & %.2f \\\\", 100 * x$proposed_unadjusted, 100 * x$reference_unadjusted),
  sprintf("Baseline-adjusted score & %.2f & %.2f \\\\", 100 * x$proposed_adjusted, 100 * x$reference_adjusted),
  "\\addlinespace", "\\multicolumn{3}{l}{\\textit{Method-specific computation time (min)}}\\\\",
  sprintf("Allocation calibration & %.2f & -- \\\\", x$calibration_minutes),
  sprintf("Analysis & %.2f & %.2f \\\\", x$proposed_analysis_minutes, x$reference_analysis_minutes),
  sprintf("Total & %.2f & %.2f \\\\", x$proposed_total_minutes, x$reference_analysis_minutes),
  "\\bottomrule", "\\end{tabular}", "\\end{table}"
)
write_tex(lines, "main_table_swift_direct_illustration.tex")

message("Table-ready TeX fragments written to: ", out_dir)
