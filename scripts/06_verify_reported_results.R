#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/"))
source(file.path(here, "06_verify_production_results.R"), local = new.env(parent = globalenv()))
