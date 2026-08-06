# Targeted pair-path supplemental sensitivity analysis

This directory contains the frozen production driver, aggregation code and publication-table generator for the 12-scenario sensitivity analysis reported in Supplementary Appendix G and Supplementary Table 5.


## Production-source provenance

`simulation_pair_path_supplemental.R` is a semantically exact R deparse of the 105,077-byte source used for the completed multi-machine run. The original source had SHA-256 `35e0adb0da0b7b584fb093cee5f4869c2ad7e79242cd7b02c0c0ea2e19988ca6`. The only post-run change in that driver is the `unname()` correction in the Wilson rejection-summary helper; this changes aggregation column names but does not alter trial generation, allocation calibration, outcome generation, SIGA calculations or the 4,999-path reference tests.

For manuscript tables, use `aggregate_pair_path_supplemental.R` or `create_pair_path_supplemental_table.R`. The publication table is portrait and intentionally omits Oracle-N and method-minus-RT columns.

## Released manuscript configuration

- 12 weak-null scenarios and two analyses per scenario;
- 100,000 outer trials per scenario;
- 4,999 regenerated allocation paths per reference randomization test;
- 100,000 three-copy analysis-calibration replicates;
- an independent 200,000-replicate direction-selection calibration for designs D5 and D6;
- 120 shards;
- base seed 20260805;
- run version `v20260805_supplemental_v2`;
- SIGA-R safeguard `max(V_R_raw, V_S/n)`.

## Run shards

Use a local output root. Each machine must receive a disjoint shard identifier from 1 through 120.

```bash
PWRT_PROFILE=manuscript \
PWRT_MODE=run \
PWRT_SCENARIO_SET=supplemental \
PWRT_PROJECT_DIR=/absolute/path/to/pair_path_supplemental \
PWRT_OUTPUT_DIR=/absolute/path/to/pair_path_supplemental/master \
PWRT_N_SHARDS=120 \
PWRT_SHARD_ID=1 \
PWRT_SEED=20260805 \
PWRT_CORES=8 \
Rscript simulation_pair_path_supplemental.R
```

## Aggregate completed shards

```bash
PWRT_RUN_DIR=/absolute/path/to/profile_manuscript_set_supplemental_M100000_B4999_B0100000_Bdir200000_eps1p00_oracle1_shape1_seed20260805_ver_v20260805_supplemental_v2_S120 \
Rscript aggregate_pair_path_supplemental.R
```

The aggregation step verifies 12 scenario directories, 1,440 shard files, 100,000 unique replicates per scenario and the frozen metadata before writing the summary, gate file and portrait Supplementary Table 5. It does not rerun outer trials or the 4,999-path tests.

## Recreate only Supplementary Table 5

```bash
Rscript create_pair_path_supplemental_table.R \
  /path/to/pair_path_supplemental_summary.csv \
  /path/to/pair_path_supplemental_table.tex
```

The publication table intentionally excludes Oracle-N and method-minus-RT columns.
