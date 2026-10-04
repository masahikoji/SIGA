# Paired covariance-reconstruction evaluation

## Current manuscript record

- `code/`: unchanged simulation code used by the supplied plans; Python 3.11--3.13.
- `code/extend_siga_run.py`: helper used to retain 2,000 trials and add 8,000.
- `results_10k/`: actual saved plans, extension records and complete 10,000-trial aggregates.
- `results/`: archived 2,000-trial aggregates (already included in results_10k).
- `report/reproduce_10k_report.py`: validate all inputs, generate manuscript blocks,
  and compare them with the included self-contained manuscript sources.
- `manuscript_20261004/`: sources corresponding to the generated 10,000-trial tables.
- `tools/prepare_10k_rerun.py`: prepare a fresh rerun from a saved plan without
  regenerating or reselecting scenarios. It does not start the simulation.

The file names `pilot` in older tools and LaTeX labels are retained for stable
references. The current manuscript reports the extended 10,000-trial totals.
Old README/status files inside `code/` describe the development snapshot; this
README and the root REPRODUCIBILITY_STATUS.md describe the supplied results.

## Reproduce tables and checks (not the outcome simulations)

From the repository root, after installing code/requirements.txt:

```bash
python validation/report/reproduce_10k_report.py
```

This reads all 7,280 rejection, 7,840 paired-comparison and 1,504 variance rows.
It retains unfavorable comparisons and tests count/rate/SE/plan/source arithmetic.
Reports are written to `validation/report/generated/10k/` by default.
No input, simulation engine or manuscript is overwritten.

## Optional full rerun

```bash
python validation/tools/prepare_10k_rerun.py --suite factorial --out "$HOME/SIGA_runs/factorial_reproduction_10k"
python validation/code/run.py calibrate --out "$HOME/SIGA_runs/factorial_reproduction_10k" --workers 4
python validation/code/run.py run --out "$HOME/SIGA_runs/factorial_reproduction_10k" --workers 4
python validation/code/run.py aggregate --out "$HOME/SIGA_runs/factorial_reproduction_10k"
```

The helper copies the saved scenarios, seed, allocations and repetition settings,
verifies the analysis signature, and creates a NEW plan with explicit rerun
provenance. It does not claim to have reused trial checkpoints. Use `controls`
or `confirm` in place of `factorial` and a separate output directory for each.
Full reruns are costly and are NOT required just to regenerate the paper tables.
The saved binary probabilities/IDs are reused instead of being regenerated with
platform-dependent floating-point rounding. Bit-for-bit agreement across software
versions and BLAS implementations is not guaranteed.

The source code's `full` profile means 100,000 trials and different inner settings;
it must not be used to reconstruct this 10,000/999/20,000 study.
