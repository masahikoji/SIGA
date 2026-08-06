#!/usr/bin/env Rscript

# Create the portrait Supplementary Table 5 from the released aggregate summary.
# Oracle-N and method-minus-RT columns are deliberately excluded.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L) {
  stop("Usage: Rscript create_pair_path_supplemental_table.R SUMMARY_CSV [OUTPUT_TEX]", call. = FALSE)
}
summary_path <- normalizePath(args[1L], mustWork = TRUE)
output_path <- if (length(args) == 2L) args[2L] else file.path(dirname(summary_path), "pair_path_supplemental_table.tex")
summary <- read.csv(summary_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c(
  "scenario_id", "scenario_code", "factor_count", "target_per_group",
  "outcome_type", "model_type", "analysis", "mean_variance_ratio_r_over_s",
  "siga_s", "siga_r", "rt"
)
missing <- setdiff(required, names(summary))
if (length(missing)) stop("Summary is missing: ", paste(missing, collapse = ", "), call. = FALSE)
fmt <- function(x, digits = 2) ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "--")
clean <- function(x) gsub("_", "-", as.character(x), fixed = TRUE)
analysis_order <- match(summary$analysis, c("unadjusted", "adjusted"))
summary <- summary[order(summary$scenario_id, analysis_order), , drop = FALSE]
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
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
writeLines(lines, output_path)
message("Wrote ", normalizePath(output_path, winslash = "/", mustWork = TRUE))
