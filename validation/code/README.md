# SIGA covariance-refinement validation

This package compares the original count-rescaled covariance with the
reconstruction that retains simulated overall/marginal imbalance moments.
It does not change the allocation rule, observed score or reference allocation
law. **No full-scale validation result is supplied or implied.**

## What has and has not been tested

The complete Python/Numba runner was executed in Linux with Python 3.13.5,
NumPy 2.3.5, SciPy 1.17.0 and Numba 0.65.1. Algebra tests and short execution
checks cover all 312 planned scenarios. Those checks use only 8 outer trials
and 49 reference draws per scenario and are NOT performance studies.

`R_adapter/` provides drop-in functions and a paired, observed-data analysis
wrapper for the author's supplied R3 engine. The R functions were reviewed
but not executed here because the environment has no R interpreter. Run its
R test before integration. Original R source snapshots are included unchanged.
The Python runner is an independent implementation, not a byte-for-byte
reproduction of the earlier R random streams or reported tables.

## Install and run a software check

Use a native Python 3.11, 3.12 or 3.13 installation. Do not copy another Mac's
virtual environment. Installation needs an internet connection for packages.
The pinned versions are the versions tested here; this is not a claim that
they are the newest versions. macOS itself was not available for testing.

```sh
bash setup.sh
source .venv/bin/activate
python run.py init --out results_smoke --suite smoke --profile smoke
python run.py calibrate --out results_smoke --workers 2
python run.py run --out results_smoke --workers 2
python run.py aggregate --out results_smoke
```

The first invocation compiles the allocation functions through Numba. No R,
C compiler or third-party R package is required for the Python runner.
`python run.py test` reruns the independent allocation/algebra/interface tests.

## Prespecified scenario families

The scenario generator does not examine rejection rates or covariance
eigenvectors to choose effect directions, seeds, sample sizes or thresholds.
The design was assembled after the earlier exploratory diagnostics; it is
**not a retrospective claim of preregistration before method development**.
The validation uses new, independently keyed random streams.

| Suite | Scenarios | Purpose |
| --- | ---: | --- |
| confirm | 224 | Four original F/n combinations, two outcome types, seven null/power roles, two configurations, first-factor and highest-order-interaction directions |
| factorial | 40 | Separate coin probability, residual SD and heterogeneity in F=5,n=400 continuous outcomes; include zero stratum heterogeneity |
| controls | 16 | Constant-effect/zero sharp-null controls across four designs and both coins/outcomes |
| correlated | 32 | The supplement's fixed correlated-profile perturbation; no covariance-optimised direction |
| all | 312 | All four suites, without deleting unfavorable scenarios |

Both adjusted and unadjusted scores are evaluated in every scenario. All
methods use the same outer trial, the same covariance simulations and the
same reference assignments. Equivalence uses the maximum of BOTH one-sided
p-values on the same trial. It is not approximated by a single boundary test.

In `confirm`, the moderate continuous configuration uses SD=1, individual-
effect SD=.25 and maximum stratum-effect deviation=.5 with p=.80. The strong
configuration uses .25, .10 and 1 respectively with p=.95. This reproduces
the combined change in the earlier stress configuration, not a causal
comparison of coin probabilities. `factorial` varies p in {.80,.95}, SD in
{.25,1}, and stratum-deviation magnitude in {0,.5,1}, with individual-effect
SD held at .25. In particular it does not secretly turn off individual
heterogeneity when stratum heterogeneity is zero.

Control risks use the same .60 marginal-risk logistic model and fixed
coefficients as the supplied R3 code. Binary effect deviations use the
original, deterministic feasibility convention: multiply the entire centred
direction by .9 until all generated treatment risks lie in (.02,.98).
This occurs while freezing the plan, before simulation. The applied factor
and actual heterogeneity are recorded for EVERY scenario in plan.csv/json.
No clipping depends on a method's performance. Existing alternative effects
and equivalence margins are retained; effects can vary between sample sizes,
so these rows must not be joined into a fixed-effect sample-size search.

## Main comparisons do not use outcome truth

The analysis API accepts observed y, treatment, stratum, baseline factors,
allocation-only covariance estimates and the TESTED boundary b. It does not
accept the generating treatment effect, true stratum effect vector, true
sampling variance, reference variance or potential outcomes.

Both original and refined SIGA-R/CRT use the SAME within-stratum observed arm-
mean differences minus b. If an arm is absent in a stratum, the manuscript's
zero-fill convention is retained and the number of estimable strata is
reported. It is not hidden or replaced by a true effect.

The DGP necessarily specifies outcome means/effects; these values are used
only to generate outcomes, declare the null/power role and audit the DGP.
No oracle-effect testing method is included in the primary output. The
reference-draw variance is a DIAGNOSTIC, computed after the proposed
estimates; it is never fed into either method's critical value.

Default `--profile-source known_design` uses the specified design-stage
profile probabilities in allocation simulations for BOTH methods. These are
not unknown treatment effects. To examine profile-law estimation separately,
initialise another plan with `--profile-source external_estimated`. It uses
an independent planning sample of max(5000,10*n) profiles and Laplace-smoothed
frequencies, not outcomes. The true and planning frequencies are recorded.
This option is a sensitivity analysis, not a claim that the profile law is
known in every completed trial.

## Profiles, workload, and staged use

| Profile | Outer trials/scenario | Reference draws/trial | Allocation replicates/design |
| --- | ---: | ---: | ---: |
| smoke | 8 | 49 | 300 |
| pilot | 2,000 | 999 | 20,000 |
| full | 100,000 | 4,999 | 100,000 |

`full` matches the high-precision counts used in the manuscript and can take
substantial time. It is NOT the default. `all/full` implies about 156 billion
complete reference assignment sequences and should not be launched blindly.
The program prints the exact workload before any simulations. Start with
software checks and a pilot, inspect paired uncertainty, and choose the
precision of a new frozen run based on that uncertainty, not on which method
wins. Full-scale output can require many GB of disk space.

A diagnostic pilot:

```sh
python run.py init --out results_factorial --suite factorial --profile pilot
python run.py calibrate --out results_factorial --workers 4
python run.py run --out results_factorial --workers 4
python run.py aggregate --out results_factorial
```

The confirmatory scenario family is requested with `--suite confirm`.
`--suite all` adds the controls and correlated laws. The plan is written and
hashed BEFORE execution. Changing code or the plan after initialisation is
rejected; use a new directory for a revised plan. Existing valid checkpoints
are reused, not overwritten. A failed trial stops its chunk; it is not removed
from the denominator.

`SIGA_BLAS_THREADS=1` is the default to prevent nested BLAS oversubscription.
Choose a moderate `--workers` value and leave memory/cores for other work.
Seeds depend on the scenario, trial and stage, not the worker or shard, so
resuming/reordering jobs does not change trial values.

For disjoint sharding across machines, first prepare the identical plan and
calibration, copy the entire output directory to each machine, then run e.g.
`--shard 1/4`, `2/4`, `3/4`, `4/4`. Use local, nonsynchronised working folders;
do not have multiple machines write the same cloud-synchronised directory.
Merge only completed nonoverlapping trial subdirectories after all machines
stop. Aggregation checks hashes, identities, missing and duplicate trials.
The software does not select favorable shards or scenarios.

Optional `--scenarios 0,1,...` facilitates debugging a planned subset. Default
aggregation refuses an incomplete plan. `--allow-partial` writes a separately
labelled partial summary and cannot be mistaken for full validation.

## Methods and interpretation of outputs

`summary/rejection_rates.csv`: rates, integer counts, binomial SEs, Wilson
intervals and nominal-level differences for every scenario and method.
`summary/paired_comparisons.csv`: paired rejection differences/SEs and both
directional discordance counts, plus differences in p-values.
`summary/variance_diagnostics.csv`: sampling and conditional-reference moment
diagnostics, within/between decomposition, rank fallbacks and variance floors.
`summary/coverage.csv`: completeness of every planned scenario.

Primary implementable comparisons are `original/refined_S_normal`,
`original/refined_R_lattice` (lattice only where applicable, normal otherwise),
and `original/refined_CRT_data`. The R-normal columns are matched normal-only
ablations. The CRT-lattice columns approximate the actual corrected rejection
region using the REFERENCE variance and rescaled threshold with raw ties;
elsewhere their first-order Gaussian approximation uses S. `CRT_noise` is an
exploratory data-based noise subtraction, not a selected primary procedure.

Closeness of SIGA-R to nominal alpha and closeness to the reference RT are
DIFFERENT targets. A conservative RT can be well reproduced by a conservative
approximation. Similarly, an accurate reference variance does not alone prove
normality of its tails. Do not report favorable nominal rates as evidence of
RT reproduction if the RT rates differ. Save both targets and all scenarios.

Outer Monte Carlo SEs condition on the simulated covariance matrices shared
across trials. They do not include all covariance-calibration uncertainty.
A second complete plan with another declared seed can assess that source;
do not choose which completed seed to publish based on its result.

## Covariance refinement and scope

The full formula and proof are in the accompanying manuscript's Appendix A.4.
It requires known asymptotic balance directions with Sigma_D L'=0. It is not
imposed on unrestricted independent allocation or arbitrary covariates.
The candidate is PSD and matches Xi on the full-rank event; if rank fails,
the whole original covariance is used. No penalty multiplier, eigen-direction
chosen by outcome performance, extra jitter, mid-p rule or significance-level
tuning is introduced. Xi is estimated from the same outcome-free allocations.
The construction does not assert that the true finite-sample cross block is
zero, or that finite-sample performance must improve in every scenario.
