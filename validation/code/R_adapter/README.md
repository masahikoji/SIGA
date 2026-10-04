# Integration into the existing R engine

Source `covariance_refinement.R` after the existing `siga_pair_path_engine.R`
and `r3_engine_additions.R`. Do not replace those files wholesale. The returned
sampling-variance list is compatible with `corrected_variance()` in the supplied
R3 add-on.

1. Augment each cached allocation calibration with
   `cal <- add_balancing_covariance(cal)`.
2. For the same observed score, calculate `siga_sampling_variance()` and
   `siga_sampling_variance_refined()`.
3. Calculate `estimate_stratum_effect_deviations()` from observed outcomes,
   assignments and strata. Pass exactly the same estimated vector to both
   reference-variance corrections.
4. Reuse the same generated reference assignment sequences for every method.
5. Record all scenarios, rank-fallback events, variance floors and paired
   rejection counts. Do not replace the original output files.

The R adapter uses the already stored first-copy covariance of D/sqrt(n) to
recover Xi = n L Cov(D/sqrt(n)) L'. That covariance has the original engine's
small numerical PSD adjustment. The Python runner accumulates Xi directly from
X=L D before any count standardisation. Both are first-order equivalent; small
floating-point differences are not Monte Carlo performance differences.

This environment did not have an R interpreter and had no network route for
installing one. The R adapter was source-reviewed and its algebra checked in
Python, but the R test script was not executed here. The full Python runner was
executed. Run `Rscript test_adapter.R` before integrating with an R simulation.

`paired_analysis.R` also supplies `compare_observed_trial()`, which evaluates
both covariance constructions on one observed trial using the same reference
sequences and data-estimated effects. It handles one boundary or the two
boundaries of equivalence, and applies the corrected lattice region when
appropriate. The input contains no generating effect or outcome model.

```r
source("original_engine/siga_pair_path_engine.R")
source("original_engine/r3_engine_additions.R")
source("covariance_refinement.R")
source("paired_analysis.R")
# Optional C acceleration:
load_rt_kernel("original_engine")
cal <- add_balancing_covariance(cal)
answer <- compare_observed_trial(y, A, X, boundaries=c(-.45,.45),
    tails=c("greater","less"), calibration=cal, alpha=.05,
    B=4999L, reference_seed=reference_seed)
```

The `original_engine` files are unchanged snapshots supplied by the author.
They contain optional model-based analysis utilities; the paired wrapper does
not call them. Do not infer that the original R3 driver's model-based primary
columns have become data-based merely by sourcing this adapter.
