#!/usr/bin/env Rscript
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])))), "R", "load_all.R"))

input <- Sys.getenv(
  "SIGA_STRESS_REPORTED",
  file.path(SIGA_ROOT, "data", "reported_results", "pair_path_stress_reported.csv")
)
output <- Sys.getenv(
  "SIGA_FIGURE_OUTPUT",
  file.path(SIGA_ROOT, "figures", "pair_path_validation_selected.pdf")
)
data <- read.csv(input, stringsAsFactors = FALSE)
selected <- c(1L, 3L, 6L, 10L, 19L, 20L)
data <- data[data$scenario_id %in% selected & data$analysis == "adjusted", , drop = FALSE]
data <- data[match(selected, data$scenario_id), , drop = FALSE]
labels <- c(
  "1: D1 continuous / homogeneous",
  "3: D1 continuous / realistic",
  "6: D1 binary / common log-odds",
  "10: D2 continuous / realistic",
  "19: D4 continuous / max-ratio stress",
  "20: D4 binary / max-ratio stress"
)
methods <- c("SIGA-S", "SIGA-R", "RT")
values <- cbind(data$siga_s_percent, data$siga_r_percent, data$rt_percent)
N <- 100000
wilson <- function(percent) {
  interval <- wilson_interval(round(percent * N / 100), N)
  100 * interval[c("lower", "upper")]
}
ensure_directory(dirname(output))
pdf(output, width = 8.8, height = 5.6, useDingbats = FALSE)
par(mar = c(4.5, 11, 1.5, 1), xaxs = "i")
y_base <- rev(seq(1, by = 4, length.out = length(selected)))
plot(NA, xlim = c(1.2, 3.15), ylim = range(y_base) + c(-2, 2),
     xlab = "Rejection probability (%)", ylab = "", yaxt = "n", bty = "l")
abline(v = 2.5, lty = 2)
abline(v = pretty(c(1.2, 3.15)), col = "grey90")
axis(2, at = y_base, labels = labels, las = 1, tick = FALSE)
point_symbols <- c(1, 0, 2)
line_types <- c(1, 1, 1)
for (i in seq_along(selected)) {
  for (j in seq_along(methods)) {
    y <- y_base[i] + (2 - j) * 0.65
    ci <- wilson(values[i, j])
    segments(ci[1L], y, ci[2L], y, lty = line_types[j])
    segments(ci[1L], y - 0.10, ci[1L], y + 0.10)
    segments(ci[2L], y - 0.10, ci[2L], y + 0.10)
    points(values[i, j], y, pch = point_symbols[j])
  }
}
legend("bottomright", legend = methods, pch = point_symbols, horiz = TRUE, bty = "n")
dev.off()
message("Wrote ", normalizePath(output, mustWork = TRUE))
