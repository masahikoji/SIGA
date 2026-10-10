# SIGA: fixed binary heterogeneity 0.12, single M3 Ultra workflow

See **README_ja.md** for the full Japanese instructions and exact project paths.

This update starts from the previously supplied SIGA_resimulation_20261009 package,
whose vendored upstream revision is recorded in UPSTREAM.json. It does not fetch
moving GitHub code or push changes to a repository. The original license is retained.

## Fixed primary design

32 scenarios: two (F,total n) pairs (2,200)/(5,400), continuous/binary outcomes,
first-factor/highest-order-interaction deviations, and superiority/NI null/power.
Allocation probability is 0.8. Both scores and both covariance constructions are
computed on the same data and reference assignments. Null evaluations use 100,000
trials; power uses 10,000; RT and CRT use 4,999 reference assignments. Each of two
allocation settings uses 100,000 independent calibration triples, once and reused.

For every primary binary DGP, max_s |(p1_s-p0_s)-Delta| is exactly 0.12 (up to floating
point arithmetic), and the profile-weighted deviations have mean zero. Both null
and alternative settings use this amplitude. The actual treatment risks in the
32-scenario plan range from 0.2442566482 to 0.9759945781. Invalid probabilities cause
an error before simulation. The active runner NEVER shrinks effects or clips risks.
The historical factory in design.py retains an explicitly named legacy policy only
for upstream regression tests; run.py / resimulation.py always select policy='error'.
Sharp-null controls, if explicitly requested through the optional suite, use zero
heterogeneity, not 0.12. Optional suites are not automatically added to the M3 run.

## Single-Mac workflow

Use standard native-arm64 CPython 3.13. Dependencies are pinned, including llvmlite.
The setup installs into ~/.venvs/siga_m3_d012, not into Conda base or a Dropbox venv.

```bash
bash m3.sh setup
bash m3.sh smoke
bash m3.sh prepare
bash m3.sh benchmark
bash m3.sh run
bash m3.sh archive
```

`prepare` is restartable and freezes the primary plan, checks every binary risk,
and computes calibration. `benchmark` times all 32 scenarios at 8/12/16/20/24
workers (subject to CPU count) and selects the fastest projected setting. This is
only a local timing estimate. `run` uses that setting; `run --workers 16` overrides
it. Without a benchmark, the fallback is min(12,CPU count). `run` also aggregates
when all trials are present. Completed chunks are verified and reused on restart.
Use `bash m3.sh status` from a second terminal to inspect progress.

Active output: ~/SIGA_runs/main_d012_20261009. Logs: ~/SIGA_runs/logs/.
Numba/Python caches: ~/.cache/siga_m3_d012/. A process lock prevents two mutating
M3 workflow commands from using the same output at once. `caffeinate -is` runs
with long commands. Do not intentionally sleep, shut down, or close the terminal
while running. Ctrl+C stops the run; restart the same command to resume.

`archive` requires a completed production summary. It makes a local ZIP containing
all checkpoints, calibration, summaries, logs, and code, then copies the archive
to SIGA_results beside this program folder. This is the requested Dropbox project
when installed using README_ja.md. Local files are NOT deleted.

## Outputs and interpretation

- plan.json / plan.csv: frozen DGP, budgets, seeds, identities, and source signature.
- binary_probability_check.csv: all per-stratum binary risks and true deviations.
- binary_probability_summary.csv / preflight.json: deterministic design checks.
- calibration/: saved allocation-only covariance estimates, shared across trials.
- trials/: checkpointed p-values, decisions and variance diagnostics.
- summary/type1_error.csv, summary/power.csv: seven primary methods, separated by goal.
- summary/paired_comparisons.csv: paired differences/MCSEs.
- summary/variance_diagnostics.csv: reference variance, sampling variance at null,
  missing-arm strata/mass and numerical fallbacks.
- summary/STATUS.json: requires complete=true AND production_budget=true for final use.

CSV rates and MCSEs are probabilities: multiply by 100 for percent/percentage points.
MCSEs condition on the shared calibration, rather than averaging over calibration
uncertainty. SIGA-R approximates RT; agreement with RT is not a size-validity result.
Thirteen inherited computational methods are saved, but the primary files contain
only RT and the base/reconstructed S-normal, R-normal and data-based CRT methods.
Lattice and noise-adjusted variants remain supplementary/exploratory.

## Optional checks and suites (not automatic)

`bash m3.sh precision` audits scenarios 0,8,16,24 with 1,000 common outer trials,
4,999 vs 9,999 reference draws and three calibration batches. It does not alter
the 4,999-draw main result. There is no automatic statistical acceptance threshold.
The audit is a deterministic prefix, not a new independent outcome study.

Direct `run.py` also retains controls, reconstruction, sample_size and pbc95 suites.
See docs/suite_workloads.csv. Use new output directories for different suites.
All heterogeneous binary settings in these suites use 0.12; controls use zero.
The sample_size suite retains source size-specific alternative effects, so it is
not a fixed-effect power curve. The pbc95 suite changes only pbc relative to matched
continuous-outcome reconstruction settings; it is not the old combined stress DGP.

## Verification and limitations

Read docs/TEST_REPORT.md. Tests and orchestration were run on Linux x86_64 with the
pinned numerical environment, not on a physical M3 Ultra. The macOS shell was syntax
checked; the native Python/architecture checks and sleep assertions require the Mac.
The full 1,760,000-trial primary simulation was NOT performed here. Do not replace
manuscript results with the included smoke/benchmark checks. No old 0.15-amplitude
results are imported. Code and numerical-environment changes invalidate existing
frozen plans rather than silently mixing results.
