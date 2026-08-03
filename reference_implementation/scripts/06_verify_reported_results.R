#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/"))
repo <- normalizePath(file.path(here, "..", ".."), winslash = "/")
target <- file.path(repo, "scripts", "06_verify_production_results.R")
status <- system2(file.path(R.home("bin"), "Rscript"), target)
quit(status = status)
