#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L) {
  stop(
    "Usage: Rscript create_pair_path_publication_figure.R SUMMARY_CSV [OUTPUT_PDF]",
    call. = FALSE
  )
}

summary_path <- normalizePath(args[1L], mustWork = TRUE)
output_path <- if (length(args) == 2L) {
  args[2L]
} else {
  file.path(dirname(summary_path), "pair_path_theory_validation_figure_publication.pdf")
}

## Scenarios labelled once in each panel. Edit this vector if a different
## subset should be highlighted in the manuscript figure.
label_scenarios <- c(2L, 8L, 9L, 15L, 16L, 19L, 20L)

## Label positions in base-R text(): 1=below, 2=left, 3=above, 4=right.
label_positions_panel_a <- c(
  `2` = 4L, `8` = 4L, `9` = 2L, `15` = 3L,
  `16` = 3L, `19` = 2L, `20` = 4L
)
label_positions_panel_b <- c(
  `2` = 4L, `8` = 3L, `9` = 2L, `15` = 4L,
  `16` = 3L, `19` = 2L, `20` = 4L
)

dat <- read.csv(summary_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c(
  "scenario_id", "analysis", "rt", "siga_s", "siga_r",
  "predicted_rt_size_from_mean_ratio", "mean_variance_ratio_r_over_s"
)
missing <- setdiff(required, names(dat))
if (length(missing)) {
  stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

numeric_columns <- c(
  "rt", "siga_s", "siga_r", "predicted_rt_size_from_mean_ratio",
  "mean_variance_ratio_r_over_s"
)
for (v in numeric_columns) dat[[v]] <- as.numeric(dat[[v]])

bad <- !complete.cases(dat[, c("scenario_id", "analysis", numeric_columns)])
if (any(bad)) {
  stop(
    "The summary contains ", sum(bad),
    " row(s) with missing or non-finite plotting values.",
    call. = FALSE
  )
}
if (!all(dat$analysis %in% c("adjusted", "unadjusted"))) {
  stop("analysis must contain only 'adjusted' and 'unadjusted'.", call. = FALSE)
}

safe_limits <- function(x, fraction = 0.06, minimum_pad = 0.0006) {
  z <- range(x, finite = TRUE)
  if (!all(is.finite(z))) stop("No finite plotting values were found.", call. = FALSE)
  pad <- fraction * diff(z)
  if (!is.finite(pad) || pad <= 0) pad <- minimum_pad
  z + c(-pad, pad)
}

point_shape <- ifelse(dat$analysis == "adjusted", 24, 21)
point_fill <- ifelse(dat$analysis == "adjusted", "black", "white")

pdf(
  output_path,
  width = 9.6,
  height = 4.7,
  family = "Helvetica",
  useDingbats = FALSE,
  onefile = TRUE
)

par(
  mfrow = c(1, 2),
  mar = c(4.4, 4.4, 2.3, 0.8),
  oma = c(1.2, 0.2, 0.1, 0.2),
  mgp = c(2.7, 0.8, 0),
  tcl = -0.25,
  las = 1,
  xaxs = "i",
  yaxs = "i"
)

## Panel A: variance-ratio prediction versus empirical RT size.
limits_a <- safe_limits(c(dat$predicted_rt_size_from_mean_ratio, dat$rt))
plot(
  dat$predicted_rt_size_from_mean_ratio,
  dat$rt,
  xlim = limits_a,
  ylim = limits_a,
  xlab = "Predicted RT rejection probability",
  ylab = "Empirical RT rejection probability",
  pch = point_shape,
  bg = point_fill,
  col = "black",
  cex = 0.88,
  main = "Variance-ratio prediction"
)
abline(0, 1, lty = 2, lwd = 1)
legend(
  "bottomright",
  legend = c("Unadjusted", "Adjusted"),
  pch = c(21, 24),
  pt.bg = c("white", "black"),
  col = "black",
  pt.cex = 0.9,
  cex = 0.78,
  bty = "n",
  inset = 0.01
)

centres_a <- aggregate(
  cbind(
    x = dat$predicted_rt_size_from_mean_ratio,
    y = dat$rt
  ),
  by = list(scenario_id = dat$scenario_id),
  FUN = mean
)
centres_a <- centres_a[centres_a$scenario_id %in% label_scenarios, , drop = FALSE]
for (i in seq_len(nrow(centres_a))) {
  sid <- as.character(centres_a$scenario_id[i])
  pos <- unname(label_positions_panel_a[sid])
  if (!length(pos) || is.na(pos)) pos <- 3L
  text(
    centres_a$x[i], centres_a$y[i], labels = sid,
    pos = pos, offset = 0.45, cex = 0.72, font = 2
  )
}

## Panel B: SIGA-R and SIGA-S versus empirical RT size.
limits_b <- safe_limits(c(dat$rt, dat$siga_s, dat$siga_r))
plot(
  NA_real_, NA_real_,
  xlim = limits_b,
  ylim = limits_b,
  xlab = "Empirical RT rejection probability",
  ylab = "Gaussian rejection probability",
  main = "SIGA-R and SIGA-S versus RT"
)
abline(0, 1, lty = 2, lwd = 1)

## Method is represented by fill; analysis is represented by shape.
points(
  dat$rt, dat$siga_r,
  pch = point_shape, bg = "black", col = "black", cex = 0.88
)
points(
  dat$rt, dat$siga_s,
  pch = point_shape, bg = "white", col = "black", cex = 0.88
)

legend(
  "topleft",
  legend = c(
    "SIGA-R, unadjusted", "SIGA-R, adjusted",
    "SIGA-S, unadjusted", "SIGA-S, adjusted"
  ),
  pch = c(21, 24, 21, 24),
  pt.bg = c("black", "black", "white", "white"),
  col = "black",
  pt.cex = 0.86,
  cex = 0.72,
  bty = "n",
  inset = 0.005
)

label_frame_b <- rbind(
  data.frame(scenario_id = dat$scenario_id, x = dat$rt, y = dat$siga_r),
  data.frame(scenario_id = dat$scenario_id, x = dat$rt, y = dat$siga_s)
)
centres_b <- aggregate(cbind(x, y) ~ scenario_id, data = label_frame_b, FUN = mean)
centres_b <- centres_b[centres_b$scenario_id %in% label_scenarios, , drop = FALSE]
for (i in seq_len(nrow(centres_b))) {
  sid <- as.character(centres_b$scenario_id[i])
  pos <- unname(label_positions_panel_b[sid])
  if (!length(pos) || is.na(pos)) pos <- 3L
  text(
    centres_b$x[i], centres_b$y[i], labels = sid,
    pos = pos, offset = 0.45, cex = 0.72, font = 2
  )
}

mtext(
  "Numbers identify selected stress scenarios; the complete scenario definitions are reported in the Supplement.",
  side = 1, outer = TRUE, line = 0.15, cex = 0.74
)

dev.off()
message("Publication figure written to: ", normalizePath(output_path, mustWork = FALSE))
