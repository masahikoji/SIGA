# Production workflows

The files in this directory are the source snapshots used to generate the
aggregate results under `../data/production_results/raw/`. They are kept
separate from the portable implementation under
`../reference_implementation/`.

## Manuscript-scale settings

| Workflow | Outer trials | RT paths/trial | Calibration paths | Base seed |
|---|---:|---:|---:|---:|
| SIGA-S continuous, unadjusted and adjusted | 100,000 | 4,999 | 100,000 one-path | 20260725 |
| SIGA-S binary, large n | 100,000 | 4,999 | 100,000 one-path | 20260726 |
| SIGA-S binary, small n | 100,000 | 4,999 | 100,000 one-path | 20260728 |
| SIGA-R full grid | 100,000 | 4,999 | 100,000 three-path | 20260801 |
| Pair-path stress study | 100,000 | 4,999 | 100,000 three-path | 20260729 |
| SWIFT DIRECT-inspired SIGA-S/RT and SIGA-R extension | 100,000 | 4,999 | 100,000 one-path plus 100,000 three-path | 20260724 |
| Standardized timing benchmark | as recorded in output CSV | 4,999 | 100,000 | 20260729 |

The full-grid SIGA-S experiment was produced by separate continuous and binary
workflows; there is therefore no single SIGA-S base seed. The randomization-
targeted full grid was an independent rerun, so its RT percentages need not be
identical to those in the sampling-targeted full grid.

## Reported-result verification

From the repository root:

```bash
Rscript scripts/06_verify_production_results.R
```

This command reads the included production aggregate CSVs and verifies the
reported numerical summaries.

## Smoke and kernel checks

```bash
Rscript production/common/check_pair_path_engine.R
Rscript production/common/00_verify_rt_kernel_equivalence.R
Rscript reference_implementation/scripts/99_smoke_test.R
```

## Rerunning production workflows

All production scripts use environment variables for output paths and Monte
Carlo settings. Use a new local output directory; do not write concurrent jobs
from different machines into the same cloud-synchronized shard directory.

### SIGA-R full grid

```bash
cd production/full_grid_siga_r
PWRT_OUTPUT_DIR=/absolute/path/to/output PWRT_N_OUTER=100000 PWRT_N_RERAND=4999 PWRT_N_CALIBRATION=100000 PWRT_SEED=20260801 PWRT_WORKERS=24 bash run_siga_r_full_benchmark_m3_ultra.sh
```

### Pair-path stress study

Run each shard independently, then aggregate using the same `PWRT_OUTPUT_DIR`,
`PWRT_N_SHARDS`, profile, scenario set, and seed:

```bash
cd production/pair_path_stress
PWRT_PROFILE=manuscript PWRT_SCENARIO_SET=all PWRT_MODE=run PWRT_PROJECT_DIR=/absolute/path/to/pair_project PWRT_OUTPUT_DIR=/absolute/path/to/pair_output PWRT_N_SHARDS=40 PWRT_SHARD_ID=1 PWRT_SEED=20260729 Rscript simulation_pair_path_stress.R

PWRT_PROFILE=manuscript PWRT_SCENARIO_SET=all PWRT_MODE=aggregate PWRT_PROJECT_DIR=/absolute/path/to/pair_project PWRT_OUTPUT_DIR=/absolute/path/to/pair_output PWRT_N_SHARDS=40 PWRT_SEED=20260729 Rscript simulation_pair_path_stress.R
```

### SWIFT DIRECT-inspired extension

First run `simulation_swift_direct_siga_s_rt.R` to create the original SIGA-S
and RT output. Then point the extension to that output and explicitly select
`scenario_01_swift_direct_NI_equal_response_power` if both power and boundary
runs are present:

```bash
cd production/swift_direct
PWRT_PROJECT_DIR=/absolute/path/to/project SWIFT_OLD_SCENARIO_DIR=/absolute/path/to/scenario_01_swift_direct_NI_equal_response_power SWIFT_R_OUTPUT_DIR=/absolute/path/to/swift_r_output PWRT_SEED=20260724 SWIFT_R_WORKERS=24 bash run_swift_direct_siga_r_extension_m3_ultra.sh
```

### Binary large-n aggregation correction

The archived large-n binary script contains a base-R name-propagation issue in
its aggregation helper: `c(probability = ci["estimate"])` creates the name
`probability.estimate`, so later indexing by `"probability"` returns `NA`.
The active script `full_grid_siga_s/simulation_binary_large_n.R` applies the
minimal `unname()` correction. The archived original source and a unified diff
are retained in the same directory. The corrected 100,000-replicate aggregate
CSV is the result used in the reported comparison.

### SIGA-S workflows and timing

The scripts under `full_grid_siga_s/` and `timing/` document their environment
variables in their headers. Set `PWRT_PROJECT_DIR`, output paths, shard counts,
and the workflow-specific seed in the table above. The large raw shard outputs
are intentionally not included in this snapshot.
