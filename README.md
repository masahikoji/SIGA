# SIGA current analysis: computational reproducibility package (10 October 2026)

This directory contains **code, fixed input plans, aggregate results and validation materials only** for the current evaluation of 32 primary plus 20 additional scenarios. It does **not** contain `main.tex`, `supplement.tex`, manuscript bibliography, PDF, or journal page-layout files. The manuscript is maintained separately for journal submission.

The current package is distinct from the historical `v1.0.0` release and older 280-scenario study. The previously published tag must remain unchanged. The current package is a dated snapshot, not yet an immutable release or a DOI deposit.

## Contents

- `code/`: Python/Numba simulation implementation, software requirements and code tests.
- `plans/`: fixed scenario definitions, identifiers, seeds and numerical budgets.
- `data/primary/`, `data/additional/`: aggregate outcomes, rejection counts and variance diagnostics (no participant data).
- `timing/`: original method-specific benchmark repetitions, calibration records and computed summaries.
- `audit/`: verifiers and table/plot recreation from archived aggregates, with frozen SHA-256 digests for the 11 archived table outputs.
- `scripts/`: rerun-plan preparation, deterministic numerical checks and exact small-sample verification.
- `review/`: additional pair-exporter utility used in the checks.

## Fast checks (no outcome resimulation)

From this directory, with the documented dependencies installed:

```sh
python scripts/verify_submission.py
python audit/rebuild_power_summary.py --write
python code/run.py test
python audit/plot_primary_null.py
```

The checks regenerate derived table fragments under `audit/regenerated_tables/` and a figure under `figures/` as **outputs**, not included manuscript sources. They do not rerun the approximately 3.04 million original outcome trials. Do not interpret them as independent verification of every result.

## Rerun the recorded simulations

Install `code/requirements.txt` in a fresh environment (recorded runtime: Python 3.13.16 on macOS arm64, NumPy 2.3.5, SciPy 1.17.0, Numba 0.65.1, llvmlite 0.47.0). Then, for example:

```sh
python scripts/prepare_saved_plan.py --study primary --out "$HOME/SIGA_runs/reproduction_primary"
python code/run.py calibrate --out "$HOME/SIGA_runs/reproduction_primary" --workers 2
python code/run.py run --out "$HOME/SIGA_runs/reproduction_primary" --workers 8
python code/run.py aggregate --out "$HOME/SIGA_runs/reproduction_primary"
```

Use `--study additional` and a separate empty output directory for the other 20 scenarios. The saved plan is checked against its frozen hash and source signature. On another numerical environment use the explicit override documented by `prepare_saved_plan.py`; cross-platform bitwise identity is not guaranteed. Full trial-level checkpoints and calibration arrays are not distributed and must be regenerated.

## Benchmarks and scope

The timing results are projections from measured single-thread benchmark repetitions, not direct runs of 100,000 repeated analyses. Preprocessing and reuse costs are documented in `timing/`. Gaussian approximations avoid per-trial reference regeneration; CRT does not. The associated clinical illustration and historical studies in the root SIGA repository are separate from the current 52-scenario evaluation.

## Archive and DOI

After reviewing this package, publish a new, immutable GitHub release without changing `v1.0.0`, and archive it with Zenodo or another DOI-issuing archive. Insert the actual new release URL and DOI in the journal manuscript before acceptance. No new release or DOI is created by this ZIP.
