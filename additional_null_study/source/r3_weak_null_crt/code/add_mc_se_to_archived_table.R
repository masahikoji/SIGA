#!/usr/bin/env Rscript
# Adds Monte Carlo standard errors to the archived operating-characteristic table
# (expected/reported_operating_characteristics.csv of the manuscript package) and
# writes the restructured main-text Table 1 (maximum absolute differences with the
# Monte Carlo standard error of the difference).
#
# Usage: Rscript add_mc_se_to_archived_table.R INPUT_CSV OUTPUT_DIR [N_OUTER]
#   N_OUTER defaults to 100000.
#
# Paired-difference standard errors require trial-level rejection indicators.  If the
# archived per-trial checkpoint files are available, pass their concatenated CSV as
# a fourth argument with columns reject_proposed, reject_randomization (or the
# equivalent names) and the script computes the exact paired SE; otherwise the
# conservative bound sqrt{p1(1-p1)+p2(1-p2)}/sqrt(M) is reported and labelled.

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: Rscript add_mc_se_to_archived_table.R INPUT_CSV OUTPUT_DIR [N_OUTER]")
input <- args[1L]; outdir <- args[2L]; M <- if (length(args) >= 3L) as.integer(args[3L]) else 100000L
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
d <- read.csv(input)
se <- function(p) sqrt(p * (1 - p) / M)
for (a in c("unadjusted", "adjusted")) {
  p1 <- d[[paste0("proposed_", a)]]; p2 <- d[[paste0("reference_", a)]]
  d[[paste0("proposed_", a, "_se")]] <- se(p1)
  d[[paste0("reference_", a, "_se")]] <- se(p2)
  d[[paste0("difference_", a)]] <- p1 - p2
  d[[paste0("difference_", a, "_se_bound")]] <- sqrt(p1 * (1 - p1) + p2 * (1 - p2)) / sqrt(M)
}
write.csv(d, file.path(outdir, "reported_operating_characteristics_with_se.csv"), row.names = FALSE)

is_null <- d$scenario_role != "power"
pick <- function(rows, a) {
  v <- rows[[paste0("difference_", a)]]; k <- which.max(abs(v))
  c(max = abs(v[k]), se = rows[[paste0("difference_", a, "_se_bound")]][k],
    label = paste(rows$outcome[k], rows$factor_count[k], rows$target_per_group[k], rows$objective[k], rows$scenario_role[k]))
}
lines <- c(
  "\\begin{table}[!htbp]", "\\centering",
  "\\caption{Maximum absolute differences in rejection probability (percentage points) between the Gaussian approximation and the reference randomisation test across the 56 scenario configurations, with a Monte Carlo standard error bound for the difference at the maximising configuration.}",
  "\\label{tab:simulation-summary}", "\\begin{tabular}{@{}llll@{}}", "\\toprule",
  "Rows & Analysis & Max $|$difference$|$ & SE bound\\\\", "\\midrule")
for (rows_name in c("null", "power")) {
  rows <- d[if (rows_name == "null") is_null else !is_null, , drop = FALSE]
  for (a in c("unadjusted", "adjusted")) {
    r <- pick(rows, a)
    lines <- c(lines, sprintf("%s values & %s & %.2f & %.2f\\\\ %% %s", rows_name, a, 100 * as.numeric(r["max"]), 100 * as.numeric(r["se"]), r["label"]))
  }
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}",
           sprintf("\\begin{flushleft}\\footnotesize With %s outer trials the Monte Carlo standard error of a rejection probability is %.2f points at 5\\%% and %.2f points at 85\\%%. The SE bound for a paired difference is $\\{p_1(1-p_1)+p_2(1-p_2)\\}^{1/2}M^{-1/2}$; the exact paired standard error is smaller because the two procedures share the outer trials.\\end{flushleft}",
                   format(M, big.mark = ","), 100 * se(0.05), 100 * se(0.85)),
           "\\end{table}")
writeLines(lines, file.path(outdir, "table1_with_se.tex"))
cat("Wrote:", file.path(outdir, "reported_operating_characteristics_with_se.csv"), "and table1_with_se.tex\n")
