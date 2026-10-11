# SIGA: focused late-stage validation (11 October 2026)

This is a **new follow-up study**, not a replacement for the existing 52-scenario results. The package narrows the supplied `SIGA_JRSSB_revision_20261011` program. It uses the original allocation and inference modules without modification; their hashes are checked against `UPSTREAM.json`. It does not contact or modify GitHub, old result folders, or manuscript files.

## Minimal Mac commands

Unzip this directory, open Terminal and enter the directory. The existing Python environment is used; nothing is automatically installed or upgraded.

```sh
cd "$HOME/Downloads/SIGA_focused_validation_20261011"
bash m3_focus.sh smoke
bash m3_focus.sh full
```

Run `full` only after `smoke` ends normally. The full run defaults to 12 worker processes, each with numerical library threads restricted to one. Override with `SIGA_WORKERS=8 bash m3_focus.sh full`. Default Python is `$HOME/.venvs/siga_m3_d012/bin/python`; override only when needed with `SIGA_PYTHON=/absolute/path/to/python`.

The script works from any source location. Results go outside the source tree, to `~/SIGA_runs/jrssb_focused_smoke_320k_20261011` or `~/SIGA_runs/jrssb_focused_full_320k_20261011`. Completed checkpoints are reused on the same command. Existing nonempty directories with a different plan/source/environment are rejected, not overwritten. On macOS, `caffeinate -i` prevents idle sleep during calibration and simulation. Screen sleep is not prohibited.

At completion, upload **`results_for_review.zip` from the full result directory**. This compact archive includes reports, fixed plans, resolved generating parameters and covariance calibration arrays. Trial checkpoints remain locally available.

## Fixed full-run scope

| Scenario | F | n | p_bc | Outcome | Null boundary | Mean effect | Independent outcome trials |
|---|---:|---:|---:|---|---:|---|---:|
| G_C_n200_null | 2 | 200 | 0.95 | continuous | -0.20 | -0.20 | 100,000 |
| G_C_n200_power | 2 | 200 | 0.95 | continuous | -0.20 | -0.20 + 1.5/sqrt(200) | 10,000 |
| G_C_n1000_null | 2 | 1,000 | 0.95 | continuous | -0.20 | -0.20 | 100,000 |
| G_C_n1000_power | 2 | 1,000 | 0.95 | continuous | -0.20 | -0.20 + 1.5/sqrt(1000) | 10,000 |
| A13_diagnostic | 5 | 400 | 0.95 | continuous | 0 | 0 | 100,000 |

Total: **320,000 unique outcome trials**. All original-method comparisons use 4,999 reference assignments per trial. Both unadjusted and main-effect-adjusted scores and both base/reconstructed covariance constructions are retained. These are paired analyses of the same observations, NOT separate independent trials.

### Large-gap settings

These use the structure of Supplementary H.1 S09/S11: independent binary-factor prevalences (0.50, 0.40), balanced target allocation and

    Y(0) = epsilon - (Delta + d_s)/2
    Y(1) = epsilon + (Delta + d_s)/2, epsilon ~ N(0, 0.05^2).

The probability-centred deviation direction is selected from **200,000 independent allocation triples**, separately for n=200 and n=1,000. Its maximum absolute value is one. Another **100,000 allocation triples** estimate the covariances for analysis. Direction selection, analysis calibration, observed assignments, outcome generation and reference assignments have separate random-number streams. The chosen direction is shared by the null/power pair at the same n. No outcome rejection rate is used to select it.

The power shift h=1.5 was inherited from the supplied review program, not fitted after the present results. It is not designed to give 80% power. The strong heterogeneity/low residual variation is a mechanism check, not a claim of typical clinical-trial conditions. Fresh random streams mean these are not the archived S09/S11 trials or exact archived directions.

### A13 diagnostic

A fresh 100,000-trial sample has the A13 generating structure: F=5, n=400, p_bc=0.95, zero stratum mean-effect deviations, residual SD 0.25 and individual-effect noise SD 0.25. Marginal factor probabilities are (0.50, 0.40, 0.30, 0.20, 0.10); the prognostic mean includes main effects and the first-by-second interaction as in the supplied generator.

This is **not a reanalysis of the original A13 checkpoints**, which are not supplied here. A pre-specified subset of 2,000 of the same outcome trials is analysed using three independent 100,000-triple covariance calibrations and reference budgets 4,999 and 9,999. The 4,999 references are the prefix of the 9,999 references; the higher-budget comparison holds calibration fixed. Calibration comparisons hold outcomes and reference paths fixed. The extra analyses must not be pooled as independent outcome trials.

All 100,000 independent A13 trials are used for the main rejection-rate evaluation. The 2,000-trial subset is used only for paired calibration/reference-budget sensitivity; it must not be treated as another independent sample. At true d=0, oracle-effect CRT is checked to equal RT trial by trial.

## Methods and outputs

`rates.csv` contains only RT, SIGA-S, SIGA-R and the **original plug-in CRT**, with base/reconstructed variants. **No external comparator and no CRT_MC are implemented or run.** `diagnostic_rates.csv` contains explicitly labelled oracle-d analyses and the already existing exploratory noise subtraction; neither becomes a new principal method.

Other outputs include paired CRT-minus-RT and SIGA-S-minus-RT comparisons, variance ratios, raw correction decompositions, missing-arm/singleton diagnostics, and A13 calibration/reference-budget sensitivity. `first_order_diagnostics.csv` evaluates null variance formulas with DGP truth and independent finite-n covariance estimates; those values are diagnostic approximations, not known limiting covariances.

## Reproducibility and resume

`plan.json` is frozen at initialization and sealed with a SHA-256 hash. `resolved_design.json` stores the selected direction and complete generating means/variances. Source signatures, calibration hashes and software versions are checked when resuming. Atomic checkpoint writes and an exclusive run-directory lock prevent incompatible concurrent writes. A force-killed process can leave a lock; inspect it and active processes before removing a genuinely stale lock. Do not change code or the plan during a run.

The source package is self-contained with respect to the selected algorithms. For an optional read-only comparison to a local SIGA checkout:

```sh
"$HOME/.venvs/siga_m3_d012/bin/python" run_focus.py audit-upstream --repo /path/to/SIGA
```

## Validation boundary

See `validation/VALIDATION.md`. Unit tests, full pipelines at smoke budgets, worker-count invariance and the production reference budgets have been checked on Linux. The **320,000-trial production study has not been executed here**, and Mac hardware has not been directly tested. No favorable result is assumed. Do not insert numerical conclusions into the manuscript before completing and reviewing the full run.
