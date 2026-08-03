#!/usr/bin/env Rscript
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])))), "R", "load_all.R"))

profile <- tolower(Sys.getenv("SIGA_PROFILE", "smoke"))
settings <- switch(profile,
  smoke = list(n_outer = 20L, B_randomization = 99L, B_calibration = 500L),
  pilot = list(n_outer = 5000L, B_randomization = 999L, B_calibration = 20000L),
  manuscript = SIGA_MANUSCRIPT$swift_direct_inspired,
  stop("SIGA_PROFILE must be smoke, pilot, or manuscript.", call. = FALSE)
)
mode <- tolower(Sys.getenv("SIGA_MODE", "run"))
output <- Sys.getenv("SIGA_OUTPUT_DIR", file.path(SIGA_ROOT, "output", "swift_direct_inspired"))
if (mode == "run") {
  run_swift_direct_inspired(
    n_outer = as.integer(Sys.getenv("SIGA_N_OUTER", settings$n_outer)),
    B_randomization = as.integer(Sys.getenv("SIGA_N_RERAND", settings$B_randomization)),
    B_calibration = as.integer(Sys.getenv("SIGA_N_CALIBRATION", settings$B_calibration)),
    base_seed = as.integer(Sys.getenv("SIGA_SEED", SIGA_MANUSCRIPT$swift_direct_inspired$base_seed)),
    output_dir = output,
    n_shards = as.integer(Sys.getenv("SIGA_N_SHARDS", "1")),
    shard_id = as.integer(Sys.getenv("SIGA_SHARD_ID", "1")),
    batch_size = as.integer(Sys.getenv("SIGA_OUTER_BATCH", "10")),
    calibration_batch = as.integer(Sys.getenv("SIGA_CALIBRATION_BATCH", "1000"))
  )
} else if (mode == "aggregate") {
  print(summarize_swift_direct(output))
} else stop("SIGA_MODE must be run or aggregate.", call. = FALSE)
