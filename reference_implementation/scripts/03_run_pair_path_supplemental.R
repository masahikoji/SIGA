#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/"))
repo <- normalizePath(file.path(here, "..", ".."), winslash = "/")
mode <- tolower(Sys.getenv("SIGA_MODE", "audit"))
script <- if (mode %in% c("audit", "calibrate", "run")) {
  file.path(repo, "production", "pair_path_supplemental", "simulation_pair_path_supplemental.R")
} else if (mode == "aggregate") {
  file.path(repo, "production", "pair_path_supplemental", "aggregate_pair_path_supplemental.R")
} else {
  stop("SIGA_MODE must be audit, calibrate, run, or aggregate.", call. = FALSE)
}
Sys.setenv(PWRT_MODE = if (mode == "aggregate") "run" else mode)
status <- system2(file.path(R.home("bin"), "Rscript"), script)
quit(status = status)
