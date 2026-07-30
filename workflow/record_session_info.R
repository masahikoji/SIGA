#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1L]] else "sessionInfo.txt"
con <- file(out, open = "wt")
on.exit(close(con), add = TRUE)
cat("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n", file = con, sep = "")
print(sessionInfo(), file = con)
cat("\nRNG settings:\n", file = con)
print(RNGkind(), file = con)
