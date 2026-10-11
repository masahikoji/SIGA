# Focused follow-up results

Profile: full; completed 320,000 / 320,000 unique trials.
Pointwise conditional Monte Carlo uncertainty; inspect independent-calibration and high-reference audits separately.

- rates.csv: the four principal procedures only; RT plus base/reconstructed SIGA-S, SIGA-R and the original plug-in CRT.
- diagnostic_rates.csv: oracle-effect analyses and the existing exploratory noise subtraction, NOT new principal procedures.
- paired_comparisons.csv: matched trial comparisons, including CRT-RT and SIGA-S-RT. Differences are probabilities; multiply by 100 for percentage points.
- variance_diagnostics.csv: ratios of means and quantiles of scale factors. VR_MC is a Monte Carlo diagnostic of the reference variance, not a new CRT implementation.
- effect_decomposition.csv: signed decomposition of the RAW variance correction, not a decomposition of rejection rates.
- audit_rates.csv and audit_paired.csv: A13 on the SAME 2,000 observations across nested reference budgets and independent calibration batches. No pooling as independent trials.
- calibration_sensitivity.csv: descriptive spread across three calibration batches. It is not a confidence interval for calibration uncertainty.
- first_order_diagnostics.csv: null-only variance formulas using DGP quantities and independent finite-n calibration estimates. They are not exact limiting variances or a proof of asymptotic sign.
- completion.csv: all planned scenarios, including unrun ones.

A13_diagnostic is a fresh diagnostic sample, not a reanalysis of the original A13 trial checkpoints, and is not a high-precision size study. The four large-gap scenarios use independent new simulations and retain h=1.5 for the power alternatives; no 80% power target was fitted. Calibration and direction selection are independent of the outcomes. No external comparator and no CRT_MC are run.
