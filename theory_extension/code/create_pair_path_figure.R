#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L) {
  stop("Usage: Rscript create_pair_path_figure.R SUMMARY_CSV [OUTPUT_PDF]", call. = FALSE)
}
summary_path <- normalizePath(args[1L], mustWork = TRUE)
output_path <- if (length(args) == 2L) {
  args[2L]
} else {
  file.path(dirname(summary_path), "pair_path_theory_validation_figure.pdf")
}

dat <- read.csv(summary_path, stringsAsFactors = FALSE)
required <- c(
  "scenario_id", "analysis", "rt", "siga_s", "siga_r",
  "predicted_rt_size_from_mean_ratio", "mean_variance_ratio_r_over_s"
)
missing <- setdiff(required, names(dat))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)

pdf(output_path, width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.2, 1.0))

limits <- range(c(dat$predicted_rt_size_from_mean_ratio, dat$rt), finite = TRUE)
pad <- 0.05 * diff(limits)
if (!is.finite(pad) || pad == 0) pad <- 0.001
limits <- limits + c(-pad, pad)
plot(
  dat$predicted_rt_size_from_mean_ratio,
  dat$rt,
  xlim = limits,
  ylim = limits,
  xlab = "RT size predicted from mean V_R/V_S",
  ylab = "Empirical RT rejection probability",
  pch = ifelse(dat$analysis == "adjusted", 17, 1),
  main = "Variance-ratio prediction"
)
abline(0, 1, lty = 2)
text(
  dat$predicted_rt_size_from_mean_ratio,
  dat$rt,
  labels = dat$scenario_id,
  pos = 3,
  cex = 0.7
)

limits2 <- range(c(dat$rt, dat$siga_s, dat$siga_r), finite = TRUE)
pad2 <- 0.05 * diff(limits2)
if (!is.finite(pad2) || pad2 == 0) pad2 <- 0.001
limits2 <- limits2 + c(-pad2, pad2)
plot(
  dat$rt,
  dat$siga_r,
  xlim = limits2,
  ylim = limits2,
  xlab = "Empirical RT rejection probability",
  ylab = "Gaussian rejection probability",
  pch = 16,
  main = "SIGA-R and SIGA-S versus RT"
)
points(dat$rt, dat$siga_s, pch = 1)
abline(0, 1, lty = 2)
legend(
  "topleft",
  legend = c("SIGA-R", "SIGA-S"),
  pch = c(16, 1),
  bty = "n"
)
text(dat$rt, dat$siga_r, labels = dat$scenario_id, pos = 3, cex = 0.7)

dev.off()
message("Figure written to: ", output_path)
