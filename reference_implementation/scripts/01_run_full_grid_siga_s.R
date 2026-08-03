#!/usr/bin/env Rscript
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])))), "R", "load_all.R"))

profile <- tolower(Sys.getenv("SIGA_PROFILE", "smoke"))
settings <- switch(profile,
  smoke = list(n_outer = 20L, B_randomization = 99L, B_calibration = 500L),
  pilot = list(n_outer = 5000L, B_randomization = 999L, B_calibration = 20000L),
  manuscript = SIGA_MANUSCRIPT$siga_s_continuous,
  stop("SIGA_PROFILE must be smoke, pilot, or manuscript.", call. = FALSE)
)
mode <- tolower(Sys.getenv("SIGA_MODE", "run"))
output <- Sys.getenv("SIGA_OUTPUT_DIR", file.path(SIGA_ROOT, "output", "full_grid_siga_s"))
scenario_text <- trimws(Sys.getenv("SIGA_SCENARIO_IDS", ""))
scenario_ids <- if (nzchar(scenario_text)) as.integer(strsplit(scenario_text, ",", fixed=TRUE)[[1L]]) else NULL
if (mode == "run") {
  run_full_grid_benchmark(
    target = "S",
    n_outer = as.integer(Sys.getenv("SIGA_N_OUTER", settings$n_outer)),
    B_randomization = as.integer(Sys.getenv("SIGA_N_RERAND", settings$B_randomization)),
    B_calibration = as.integer(Sys.getenv("SIGA_N_CALIBRATION", settings$B_calibration)),
    base_seed = as.integer(Sys.getenv("SIGA_SEED", SIGA_MANUSCRIPT$siga_s_continuous$base_seed)),
    output_dir = output,
    n_shards = as.integer(Sys.getenv("SIGA_N_SHARDS", "1")),
    shard_id = as.integer(Sys.getenv("SIGA_SHARD_ID", "1")),
    scenario_ids = scenario_ids,
    batch_size = as.integer(Sys.getenv("SIGA_OUTER_BATCH", "10")),
    calibration_batch = as.integer(Sys.getenv("SIGA_CALIBRATION_BATCH", "1000"))
  )
} else if (mode == "aggregate") {
  print(summarize_full_grid(
    output, target = "S",
    n_outer = as.integer(Sys.getenv("SIGA_N_OUTER", settings$n_outer)),
    allow_partial = identical(Sys.getenv("SIGA_ALLOW_PARTIAL", "0"), "1")
  ))
} else stop("SIGA_MODE must be run or aggregate.", call. = FALSE)
