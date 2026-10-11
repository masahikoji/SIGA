# Completed follow-up validation (11 October 2026)

This is an additive code-and-results snapshot for the follow-up in Supplement H.4 and Tables 15--17. It contains no manuscript TeX sources, bibliography, manuscript PDFs or font files. Table fragments are generated outputs, not manuscript sources.

## Completed design

- Two-factor strong-heterogeneity models, n=200 and n=1000: 100,000 null trials and 10,000 alternative trials each (220,000 total).
- Fresh A13 evaluation: 100,000 trials.
- Principal reference budget: 4,999 per trial.
- A13 calibration and nested 9,999-reference sensitivity: the same first 2,000 trials, not additional independent trials.
- Total: 320,000 independent trials. All planned trials completed.

## Contents and checks

`program/SIGA_focused_validation_320k_20261011/` is the exact frozen execution source, with upstream allocation/inference files kept unchanged.
`results/` contains the frozen plan, resolved directions, seven calibration arrays, CSVs and execution metadata.
`results_for_review.zip` preserves the original uploaded archive byte for byte.

From this folder, in the recorded SIGA environment:

```sh
python verify_followup_results.py
python rebuild_followup_tables.py
python program/SIGA_focused_validation_320k_20261011/run_focus.py test
```

The first command verifies supplied aggregates/calibrations, not missing trial checkpoints. The second recreates the three paper table fragments. Neither reruns the simulations.

## Important report-text correction

The original `results/report/README.md` contains the obsolete sentence that A13 "is not a high-precision size study". The sealed plan and completion.csv establish 100,000 A13 trials. Only its paired sensitivity subset contains 2,000 trials. The frozen source and raw archive are retained unchanged for provenance; this note corrects the report wording, not the statistical calculations.

## Reproduce the full simulation

Use the dependencies recorded in results/plan.json and in the frozen program. Start a fresh output directory; do not overwrite the recorded results. For example, on the original Mac environment:

```sh
cd program/SIGA_focused_validation_320k_20261011
SIGA_FOCUS_OUT="$HOME/SIGA_runs/followup_reproduction_20261011" bash m3_focus.sh full
```

The original program checks its environment and hashes. Across machines numerical bitwise identity is not guaranteed. Trial-level checkpoints and reference draws are not distributed; rerunning reconstructs them. No new simulation is needed to use the completed results.

## Interpretation

The gap models deliberately examine a sizeable variance difference. Their alternatives use h=1.5, not an 80% power target. Principal procedures estimate effects from observations; true effects enter only explicit diagnostics. The A13 follow-up is a fresh sample and does not replace the original A13. Calibration sensitivity on 2,000 shared trials is not an unconditional precision guarantee.

Upload this entire `followup_20261011` folder into the existing SIGA repository without altering historical tags. Add a link to it in the repository README. No release or DOI is created by this package. The root historical reproducibility-status document should be labelled as an earlier study rather than the current follow-up status.
