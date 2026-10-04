# R3 simulation: weak-null heterogeneity and the variance-corrected reference test

This add-on is placed beside `theory_extension/` of the manuscript code supplement.
It does not modify any existing program. `code/siga_pair_path_engine.R` is an
unchanged copy of the theory-extension engine (allocation, three-copy calibration,
outcome models, SIGA-S, SIGA-R, reference test). The new functions are in
`code/r3_engine_additions.R`; the driver is `code/simulation_r3_weak_null_crt.R`.
Base R only.

## What it adds (referee point R3)

1. **Heterogeneous-effect weak-null scenarios in the four main-study designs**
   (K2 n=200, K2 n=1000, K5 n=400, K5 n=2000; p_bc = 0.80), for continuous and
   binary outcomes, at every null value of the main study (superiority 0,
   non-inferiority margin, both equivalence limits) and at the manuscript's fixed
   alternatives. Individual effects are `Delta + d_s + eta_i` with pi-weighted
   mean of `d` equal to zero, so the marginal effect equals the manuscript value.
2. **The variance-corrected reference test (CRT)**: the reference test with the
   critical value rescaled by `lambda = sqrt(V_R / V_S)`, in three versions:
   model-based `d(b)`, data-based `d_hat(b)` (stratum mean differences), and a
   noise-debiased data-based version. With `d = 0` the CRT is the reference test.
3. **Monte Carlo standard errors** for every rejection probability and for every
   paired difference (methods share the outer trials), plus the empirical variance
   of the regenerated statistics for comparison with `V_R`.
4. **Tables**: null-value and power tables with SEs, and a restructured Table 1
   (maximum absolute differences with the SE of the difference).
   `code/add_mc_se_to_archived_table.R` adds SEs to the archived main-study table.

Tie convention: a regenerated statistic tied with the observed statistic (|T*| = |T_n|, or
T* = T_n for one-sided tests) always counts as at least as extreme, as in the reference test.
Without this clause a rescaling factor lambda > 1 would exclude exact ties, which are frequent
for lattice-valued scores (binary outcomes), and the corrected test would become anti-conservative
in finite samples although lambda - 1 is negligible; the clause is asymptotically immaterial.

Methods on common outer trials: `rt`, `crt_model`, `crt_data`, `crt_data_db`,
`siga_s`, `siga_r`; analyses: `unadjusted`, `adjusted` (baseline-only residuals).
Objectives: superiority (two-sided, alpha .05), non-inferiority (one-sided, .025),
equivalence (two one-sided tests at .05, intersection-union). Equivalence uses the
same regenerated paths for both limits, as in the manuscript.

## Scenario sets (`R3_SCENARIO_SET`)

| set | rows | models |
|---|---|---|
| `null` (default) | 32: 4 designs x 2 outcomes x 4 null values | continuous: `Y(0)=mu_s+e`, `Y(1)=Y(0)+Delta+d_s+eta`, direction = centred first factor, `max|d_s|=0.5`, `sd(eta)=0.25`; binary: `p1_s = p0_s + Delta + d_s`, `max|d_s|=0.15` (shrunk if a probability would leave (0.02, 0.98); the shrink factor is recorded) |
| `power` | 24: 4 x 2 x 3 alternatives | same models at the manuscript alternatives |
| `interaction` | 32 | as `null` with the highest-order-interaction direction (the least favourable direction in the exact computations) |
| `manuscript` | 56 | the original models: homogeneous continuous effect; common-log-odds binary (reproduces the main study with CRT added) |
| `all` | 168 | everything (also interaction x power) |

Scenario ids are fixed in the 168-row grid (`r3_scenario_definitions.csv`), so
`R3_SCENARIO_IDS=2,8,14` selects by id regardless of the set.

## Compiled regeneration kernel (strongly recommended)

The reference test is the only expensive step: every outer trial regenerates `B` paths,
each sequential over the `n` participants, so one scenario at manuscript scale costs
`100,000 x 4,999 x n` sequential assignment steps. `code/rt_kernel.c` implements this loop
in C through base R's `.Call` interface (no Rcpp needed) and draws the uniforms from R's
own generator in the same order as the pure-R loop, so the regenerated statistics and all
p-values are **bit-identical** to the pure-R implementation (verified by
`check_r3_engine.R`). It is compiled automatically on first use (`R CMD SHLIB`) into
`code/build/<platform>_R<version>/` and used whenever a C compiler is available.

- macOS: install the Xcode command-line tools once (`xcode-select --install`).
- Linux: `build-essential` (gcc) and the R development headers (`r-base-dev`).
- `R3_KERNEL=R` forces the pure-R loop; `R3_KERNEL=C` fails loudly if compilation fails;
  `R3_BUILD_DIR` moves the build directory (e.g. outside a synchronised folder).

Measured per trial with `B = 4,999` (this is the whole reference test for one trial,
including all score columns; one core of a 2024 laptop-class CPU):

| design | pure R | C kernel |
|---|---|---|
| K2 n=200 | 0.27 s | 0.07 s |
| K2 n=1000 | 0.80 s | 0.14 s |
| K5 n=400 | 0.47 s | 0.07 s |
| K5 n=2000 | 2.4 s | 0.36 s |

On an Apple M3 Ultra the absolute times are roughly half of these. One 32-scenario set at
manuscript scale is therefore about 70 core-hours with the C kernel (about 3 hours on
24 cores) instead of about 400 core-hours in pure R.

## Environment variables

| variable | default | meaning |
|---|---|---|
| `R3_PROFILE` | `pilot` | `smoke` (20/99/500), `pilot` (2000/999/20000), `manuscript` (100000/4999/100000): outer trials / regenerated paths / calibration paths |
| `R3_MODE` | `run` | `calibrate`, `run`, `aggregate` |
| `R3_SCENARIO_SET` | `null` | see above |
| `R3_SCENARIO_IDS` | (all in set) | comma-separated ids |
| `R3_N_OUTER`, `R3_N_RERAND`, `R3_N_CALIBRATION` | profile | explicit sizes |
| `R3_N_SHARDS`, `R3_SHARD_ID` | 1, 1 | sharding (replicates `SHARD_ID, SHARD_ID+N_SHARDS, ...`) |
| `R3_CORES` | min(8, cores-1) | `parallel::mclapply` workers (Unix); set it explicitly on a large machine |
| `R3_OUTER_BATCH` | 10 x cores | trials per parallel batch and checkpoint; must be at least the core count, otherwise only `R3_OUTER_BATCH` cores are used |
| `R3_PROJECT_DIR` | `../r3_output` | calibration cache and outputs |
| `R3_PBC` | 0.80 | biased-coin probability for all designs (e.g. 0.95 for a stress run) |
| `R3_CONT_SD`, `R3_CONT_MAX_D`, `R3_CONT_ETA_SD` | 1.0, 0.5, 0.25 | continuous outcome SD, heterogeneity size, individual-effect SD |
| `R3_BIN_MAX_D` | 0.15 | binary heterogeneity size (risk-difference scale) |
| `R3_BINARY_CONTINUITY` | 0 | 1 = use the Appendix-D lattice continuity correction for unadjusted binary superiority (SIGA-S/-R only) |
| `R3_SKIP_DETAILS` | 0 | 1 = do not write the (large) trial-level details CSV at aggregation |
| `R3_ALLOW_PARTIAL` | 0 | 1 = aggregate incomplete runs |
| `R3_SEED` | 20261002 | base seed |
| `R3_KERNEL` | `auto` | `auto` / `C` / `R`: regeneration loop implementation (see above) |
| `R3_BUILD_DIR` | `code/build/...` | where the compiled kernel is placed |

The run tag (output sub-directory) records the sizes and the heterogeneity settings,
so runs with different settings never mix.

## Terminal instructions

```bash
cd r3_weak_null_crt

# 0. checks (seconds) and smoke test (about a minute)
bash workflow/run_checks.sh
bash workflow/run_smoke_test.sh

# 1. pilot on one machine (2,000 trials, 999 paths; null set; a few hours on 8 cores)
R3_CORES=8 bash workflow/run_pilot.sh null

# 2. manuscript scale, sharded over machines/cores: run every shard exactly once
#    (calibration for each design is created once and cached in R3_PROJECT_DIR)
bash workflow/run_calibration_only.sh manuscript
R3_CORES=8 bash workflow/run_full_shard.sh 40 1 null
R3_CORES=8 bash workflow/run_full_shard.sh 40 2 null
#   ...
R3_CORES=8 bash workflow/run_full_shard.sh 40 40 null

# 3. aggregate after all shards are complete
bash workflow/aggregate_full_results.sh 40 null

# 4. other sets (same pattern)
R3_CORES=8 bash workflow/run_full_shard.sh 40 1 power
R3_CORES=8 bash workflow/run_full_shard.sh 40 1 interaction
R3_CORES=8 bash workflow/run_full_shard.sh 40 1 manuscript

# 5. stress variant showing the correction at work (stronger coin, less noise)
R3_PBC=0.95 R3_CONT_SD=0.25 R3_CONT_MAX_D=1.0 R3_CONT_ETA_SD=0.10 \
  R3_PROJECT_DIR=$PWD/r3_output_stress R3_CORES=8 bash workflow/run_full_shard.sh 40 1 null

# 6. Monte Carlo standard errors for the archived main-study table
Rscript code/add_mc_se_to_archived_table.R \
  ../expected/reported_operating_characteristics.csv r3_output/archived_se 100000
```

Shards are restartable: completed replicates are read from the shard file and skipped.
Progress and completion: the log ends with `Shard 1 of 1 completed` when a run is finished;
`bash workflow/progress.sh LOGFILE` prints the scenarios finished and the latest progress line,
and `nohup bash workflow/notify_when_done.sh LOGFILE &` raises a macOS notification (and a
spoken message) when the run ends or if the R process disappears.
Use a separate `R3_PROJECT_DIR` per machine if there is no shared storage and copy the
`scenario_shards/` directories together before aggregating.

Approximate cost at manuscript scale: see the kernel table above. With the C kernel a
32-scenario set is about 70 core-hours (about 3 hours on 24 cores); in pure R about
400–500 core-hours. The pilot profile is about 1/250 of a manuscript-scale set.

## Outputs (in `R3_PROJECT_DIR/simulation_output/<run tag>/`)

- `r3_scenario_definitions.csv`: models, null values, `max|d|`, shrink factors, model pair ratios `d'Psi d / d'Pi d`.
- `r3_design_diagnostics.csv`: calibration summaries and min/max generalised ratios.
- `scenario_shards/scenario_<id>_<label>/shard_*.csv`: trial-level checkpoints.
- `r3_<set>_summary.csv`: per scenario and analysis: rejection probabilities with SE
  (`<method>`, `<method>_se`), paired differences with SE (`diff_<method>_minus_rt`,
  `diff_siga_s_minus_crt_data`), deviations from the nominal level for null rows,
  mean `lambda` for each CRT version, mean `V_R/V_S`, mean empirical-RT-variance
  ratios, mean error of `d_hat`, floor-activation rate, timing.
- `r3_<set>_details.csv`: all trial-level rows (omit with `R3_SKIP_DETAILS=1`).
- `r3_<set>_null_table.tex`, `r3_<set>_power_table.tex`: rejection probabilities (%) with SE.
- `r3_<set>_table1_summary.tex`: maximum absolute differences with the SE of the
  difference at the maximising scenario (restructured Table 1).

## Notes

- The data-based `d_hat(b)` uses within-stratum arm means; strata with an empty arm
  contribute zero. With 32 strata and n = 400 the plug-in quadratic form is inflated by
  estimation noise, which makes the plug-in CRT slightly conservative; `crt_data_db`
  subtracts the estimated noise contribution and is the recommended data-based version.
  The summary reports the mean of each `lambda` so the three versions can be compared.
- Under the sharp null (`manuscript` set, continuous rows) `lambda_data^2 = 1 + O_p(1/n)`,
  so CRT-data and RT agree up to that order; CRT-model equals RT exactly there.
- `mean_ratio_vrt_empirical_over_vr_model` close to one confirms that `V_R` tracks the
  conditional variance of the regenerated statistics; `mean_ratio_vrt_empirical_over_vs`
  is the realised `V_R/V_S` ratio.
- The default heterogeneity (`max|d|=0.5` with outcome SD 1) gives `V_R/V_S` of about
  1.005 in the K2 designs, i.e. the correction is small, which is itself a finding for the
  main study; the `interaction` set and the stress variant make it larger.
