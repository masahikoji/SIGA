# Pair-path theory-validation extension

This extension supplements the broad manuscript reproduction programs and does not alter them. It compares SIGA-S, SIGA-R, and the regenerated-path fixed-score randomization test under aligned, practical, and deliberately heterogeneous scenarios.

## Methods represented in the code

- **SIGA-S**: trace-matched sampling-calibrated variance;
- **SIGA-R**: SIGA-S plus `n d'(Psi - Pi_hat) d / 16`;
- **RT**: fixed-score conditional randomization test under the same biased-coin minimization rule.

The manuscript-scale defaults are 100,000 outer trials, 4,999 rerandomizations, and 100,000 three-copy calibration replicates.

## Run a smoke test

```bash
bash workflow/run_smoke_test.sh
```

## Run manuscript-scale shards

```bash
bash workflow/run_full_shard.sh 40 SHARD_ID core
bash workflow/aggregate_full_results.sh 40 core
bash workflow/run_full_shard.sh 40 SHARD_ID sensitivity
bash workflow/aggregate_full_results.sh 40 sensitivity
```

## Verify the reported table

```bash
Rscript workflow/verify_pair_path_reported_results.R \
  path/to/pair_path_theory_validation_summary.csv
```

Frozen displayed values are in `expected/`. The generalized-eigenvector cases are diagnostic stress tests, not clinically representative outcome models. They show the direction of the sampling-versus-randomization variance separation. The largest reported SIGA-R--RT discrepancy was 0.68 percentage points in a finite-sample synthetic stress direction, so the asymptotic correction is not a uniform finite-sample guarantee.

The study supplies numerical evidence for the pair-path structure but does not prove product-chain regularity for the weighted absolute-imbalance allocation rule.
