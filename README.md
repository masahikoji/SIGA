# SIGA: sampling and randomization calibration under biased-coin minimization

This repository accompanies:

> Kojima, M. *Fast Power Evaluation under Biased-Coin Minimization: Sampling and Randomization Calibration*.

It implements the historical equal-weight absolute-range Pocock-Simon
allocation rule used in the study. It does not substitute a quadratic-
potential allocation criterion.

## Repository structure

- `production/` contains the source snapshots used for the reported numerical
  studies.
- `data/production_results/raw/` contains the corresponding aggregate outputs.
- `scripts/06_verify_production_results.R` recomputes the headline numerical
  summaries from those aggregate outputs.
- `reference_implementation/` contains a consolidated base-R implementation
  for inspection, smoke tests, and small independent reruns.

The production workflows and portable reference implementation are kept
separate so that the provenance of the reported results is explicit. The
binary large-n production workflow includes a documented aggregation
correction for a base-R name-propagation issue; the archived original source
and a unified diff are retained with the active script.

## Included studies

- SIGA-S sampling-targeted full grid: continuous and binary outcomes, 56
  scenario roles and two analyses, using one-path calibration.
- Independent SIGA-R randomization-targeted full grid: the same 56 roles and
  two analyses, using three-path calibration and the `1/n` safeguard.
- Twenty-scenario pair-path stress study.
- SWIFT DIRECT-inspired prospective design simulation based on published
  aggregate planning characteristics.
- Standardized timing and one-path calibration benchmarks.

The two full-grid experiments are independent outer Monte Carlo runs. Methods
within each experiment use common outer trials; RT percentages may therefore
differ slightly between the SIGA-S and SIGA-R tables without contradiction.

## Requirements

- Base R 4.2 or later is recommended.
- No contributed R packages are required by the included simulation code.
- Manuscript-scale reruns are computationally intensive and require sharding.

## Verify the included results

From the repository root:

```bash
bash scripts/00_verify_manifest.sh
Rscript scripts/06_verify_production_results.R
```

The R verifier checks directly that:

- the SIGA-S full grid has 112 comparisons and maximum absolute differences
  0.422, 0.185, and 0.422 percentage points overall, at null boundaries, and
  for power;
- the SIGA-R full grid has 112 comparisons and corresponding maxima 0.322,
  0.205, and 0.322 percentage points;
- the pair-path stress maxima are 0.825 and 0.685 percentage points for SIGA-S
  and SIGA-R;
- the SWIFT DIRECT-inspired powers and 100,000-trial audit match the reported
  results;
- the included standardized CSVs are reproducible from the raw production
  aggregates.

## Regenerate standardized result files

```bash
Rscript scripts/08_make_standardized_results.R
```

## Smoke checks

```bash
Rscript production/common/check_pair_path_engine.R
Rscript production/common/00_verify_rt_kernel_equivalence.R
Rscript reference_implementation/scripts/99_smoke_test.R
```

## Production reruns

See `production/README.md`. Do not use the same cloud-synchronized shard
output directory from multiple machines. Use disjoint shard identifiers and
merge only after completion.

## Method conventions

- Equal-weight sum of absolute overall and active marginal imbalances.
- Probability 0.5 for exact allocation ties.
- Inclusive plus-one Monte Carlo randomization-test p-values.
- Lattice-normal mixture only for unadjusted binary superiority at boundary
  zero in SIGA-S; Gaussian tails otherwise and for all reported SIGA-R tests.
- Two-sided absolute-value superiority tests, upper-tail non-inferiority tests,
  and intersection-union equivalence tests.
- SIGA-R safeguard `max(V_R_raw, V_S/n)`.

## Theoretical scope

The numerical implementation exactly follows the absolute-range rule. The
analysis retains the required one-copy and three-copy allocation limits as
explicit assumptions; simulation is not claimed to prove the unresolved
global stability and additive-functional CLT steps.

## Citation and license

Citation metadata are in `CITATION.cff`. The included license is restrictive;
replace it before public release if permissive reuse is intended.
