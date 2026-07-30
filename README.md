# SIGA reproducibility repository

This repository contains the simulation programs and frozen reported-value files for the manuscript

> **Fast Approximation of Randomization-Test Power under Biased-Coin Minimization**  
> Masahiro Kojima

The code implements the broad SIGA-S versus fixed-score randomization-test benchmark, the SWIFT DIRECT-inspired illustration, timing calculations, and the targeted pair-path validation of SIGA-S and SIGA-R.

## Repository identifier

Replace these placeholders when the public repository and archive are created:

- Repository: `https://github.com/USERNAME/REPOSITORY`
- Release: `v1.0.0`
- Archived release DOI: `https://doi.org/10.5281/zenodo.XXXXXXX`

## Requirements

- R with `Rscript` available;
- base R only; no contributed R packages are required;
- Bash for the supplied workflow scripts;
- Linux or macOS is recommended for multicore execution; Windows uses sequential execution.

The full manuscript calculations are computationally intensive: the production profile uses 100,000 outer trials, 4,999 regenerated allocations per randomization test, and 100,000 allocation-only calibration replicates.

## Quick validation

Run the reduced production-engine checks:

```bash
bash workflow/run_smoke_tests.sh
```

Run the reduced pair-path validation:

```bash
cd theory_extension
bash workflow/run_smoke_test.sh
```

The smoke tests verify execution only and do not reproduce manuscript estimates.

## Broad manuscript benchmark

Run every shard exactly once:

```bash
bash workflow/run_full_shard.sh N_SHARDS SHARD_ID
```

After all shard files are available:

```bash
bash workflow/aggregate_full_results.sh N_SHARDS
```

Run the standardized timing benchmark separately on the target machine:

```bash
bash workflow/run_timing_benchmark.sh
```

Build and compare the manuscript-level outputs with the frozen values in `expected/`:

```bash
SIGA_VERIFY_TIMING=1 bash workflow/finalize_manuscript_outputs.sh
```

## Pair-path theory validation

The extension distinguishes:

- **SIGA-S**, which uses the trace-matched sampling-calibrated variance;
- **SIGA-R**, which adds the pair-path correction and targets the conditional fixed-score randomization distribution;
- **RT**, the regenerated-path fixed-score randomization test.

For a 40-shard production run:

```bash
cd theory_extension
bash workflow/run_full_shard.sh 40 1 core
# repeat with shard IDs 2,...,40
bash workflow/aggregate_full_results.sh 40 core

bash workflow/run_full_shard.sh 40 1 sensitivity
# repeat with shard IDs 2,...,40
bash workflow/aggregate_full_results.sh 40 sensitivity
```

After aggregation, verify displayed manuscript values:

```bash
Rscript workflow/verify_pair_path_reported_results.R \
  path/to/pair_path_theory_validation_summary.csv
```

The targeted study includes aligned controls, practical outcome models, heterogeneous-effect stress directions, correlated profiles, and a stronger biased coin. In the practically motivated models the two variance targets were nearly aligned. In deliberately amplified stress directions, SIGA-R reproduced the main predicted shift but was not uniformly closer to RT in every finite-sample scenario; the largest absolute SIGA-R--RT difference in the reported stress study was 0.68 percentage points. The asymptotic result should therefore not be interpreted as a uniform finite-sample error bound.

## Directory map

- `code/`: frozen production programs for the broad benchmark and case study;
- `workflow/`: broad benchmark, aggregation, timing, and verification scripts;
- `config/`: frozen simulation configuration;
- `expected/`: values reported in the manuscript and Supplement;
- `theory_extension/code/`: pair-path calibration and validation engine;
- `theory_extension/workflow/`: pair-path smoke, sharding, aggregation, and verification scripts;
- `figures/`: final pair-path figure and its source data;
- `MANIFEST.csv` and `SHA256SUMS`: file inventory and checksums.

## Reproducibility notes

The workflows create `sessionInfo.txt`; retain it with archived outputs. The exact R version used for the original broad simulation was not recorded in the source archive, so reported values are checked at displayed precision rather than by requiring bitwise-identical output across platforms. Random-number streams and manuscript-scale settings are fixed in the programs.

## Theoretical scope

The numerical implementation uses weighted absolute-imbalance biased-coin minimization. The manuscript's allocation-copy limit is stated under product-chain regularity. The simulations and calibration diagnostics do not replace an algorithm-specific proof of that condition for the weighted absolute-imbalance rule.

## License and citation

The code is released under the MIT License. See `CITATION.cff` for citation metadata after the repository placeholders have been replaced.
