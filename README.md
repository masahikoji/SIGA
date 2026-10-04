# SIGA: sampling and conditional randomisation inference

This working revision accompanies Masahiro Kojima's manuscript:

> Squared assignment correlations and randomisation tests of weak null hypotheses under covariate-adaptive allocation.

The historical `v1.0.0` tag identifies the original paper version,
*Fast Power Evaluation under Biased-Coin Minimization: Sampling and Randomization
Calibration*. Preserve that tag and its original source/output snapshots.
The current working version is `1.1.0-dev`; no new public release or DOI is implied.

## Original studies

`production/`, `data/production_results/raw/`, `reference_implementation/` and
`scripts/06_verify_production_results.R` retain the original R workflows and
aggregate checks. The two original full grids used independent outer trials,
with common trials for methods within each grid. The old 100,000-trial results
are not results for the later covariance reconstruction.

## Revision contents

- `validation/code/`: unchanged, observed-data-only paired evaluation code and R adapter.
- `validation/results_10k/`: complete outputs for 280 scenarios, each with 10,000
  trials (2,800,000 trials total), saved plans and extension records.
- `validation/results/`: archived 2,000-trial pilots already included in the 10,000 totals.
- `validation/report/`: checks and deterministic table reproduction.
- `validation/manuscript_20261004/`: the matching manuscript sources.
- `additional_null_study/`: earlier R3 source/results and finite-state programs.
- `theory_checks/`: supplied finite-condition and allocation-direction checks.
- `REPRODUCIBILITY_STATUS.md`: verified scope and remaining provenance limitations.

## Statistical and theoretical scope

The manuscript compares sampling and conditional randomisation variances under an
explicit joint allocation condition. The additional covariance is determined by
squared assignment correlations given the ordered profiles; the variance difference
involves that covariance minus the diagonal matrix of profile probabilities.
Average-effect inference, approximation of the original randomisation test and
variance rescaling of that test are distinct objectives. The corrected test retains
original-scale ties; finite-sample sharp-null exactness is not established.

Joint allocation limits are verified for fixed-size stratified permuted blocks,
stratified Efron biased coins, and stochastic equal-weight absolute-range
Pocock--Simon minimisation including overall balance, for every positive joint law
of two binary factors and more factors satisfying stated finite inequalities.
This is not a result for every minimisation variant or every adaptive design.

The covariance reconstruction retains separately simulated overall/marginal moments
and has the same first-order target for verified balance directions. The 10,000-trial
study supports reduced variance overestimation in the motivating unadjusted settings;
not all variance or tail-probability comparisons improve. Complete comparisons are
included, with paired Monte Carlo standard errors. No new timing claim is inferred.

## Verification from the repository root

Install the dependencies in `validation/code/requirements.txt`, then run:

```bash
python validation/report/reproduce_10k_report.py
python validation/code/run.py test
```

The first command checks/reproduces tables from supplied aggregates, not 2.8 million
new outcome trials. The second checks software/algebra, not performance superiority.
After applying and reviewing this overlay, refresh the root checksum manifest:

```bash
python scripts/refresh_manifest.py --write
python scripts/refresh_manifest.py --check
```

Original results can still be checked with:

```bash
Rscript scripts/06_verify_production_results.R
```

See `validation/README.md` for a fresh rerun using the exact saved 10,000-trial plan.
Keep large active run directories outside cloud-synchronised folders.

## Provenance and publication

The extension retains trial indices 0--1999 and adds 2000--9999 without changing
methods, scenarios, master seed, inner draws or allocation estimates. This is an
increase in Monte Carlo precision, not an independent confirmatory experiment.
Actual plans and source signatures now reconcile the earlier pilot identifiers.
Raw checkpoint files, actual calibration arrays and runtime logs were not supplied
in this update; aggregate checks do not authenticate every originating execution.

The older R3 strong-setting CSV is still explicitly reconstructed from rounded
aggregates; the new 10,000-trial data do not restore that separate original output.
The licence is unchanged. Finalise a new immutable release, archive it with a DOI,
and update the manuscript's code citation before acceptance. Do not invent a
release URL or move the existing v1.0.0 tag.
