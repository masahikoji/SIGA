# Final delivery audit: methodology update and validation code

## What is delivered

Full main.tex and supplement.tex add the allocation-only covariance
reconstruction, algebraic moment-matching and recoding invariance, rank fallback,
and first-order equivalence proof. The optional stratum-effect noise correction
and the corrected binary lattice rejection region are also documented with
their conditions. No performance ranking or full-scale result is manufactured.

The base is the last integrated manuscript, not the old text snapshots. Every
existing numerical table body and label sequence was compared with that base.
The numerical sections explicitly retain the ORIGINAL covariance construction.
No numerical entries or timing values were repurposed for the new method.

## Fairness and role of the earlier pilot

The earlier independent diagnostic pilot used model-known stratum deviations
to isolate variance error. It is not an implementable-method comparison.
The current primary Python analysis accepts only observed outcomes, assignments,
profiles, tested boundaries and allocation-only covariance estimates. Both
methods use the same data-based deviation estimator. True DGP means/effects
are used only by the generator and DGP audit, not by inferential methods.

The fixed validation suites include the previous design combinations, both
scores, two effect directions, factorial separation, sharp-null controls and
correlated profile laws. This is a post-development validation plan, not a
claim of preregistration before the exploratory pilot. There is no outcome-
based tuning, favorable-seed search, covariance-eigenvector effect selection
or silent omission of failed trials. Calibration/outer/reference streams are
separately keyed; every comparison uses common simulated trials.

Known design-stage profile probabilities are common inputs to both methods.
The external-estimated-profile option uses independent planning profiles.
No outcome-truth testing method is included in primary output. Optional
noise subtraction is retained separately as an exploratory data-based column.

## Execution checks

Python unit/algebra tests passed. All 312 scenarios completed software checks
with 8 outer trials, 49 reference assignments and 300 calibration replicates.
These 2496 trials are SOFTWARE checks, not evidence of method superiority.
A separate four-scenario external-profile check completed using two disjoint
shards. The end-to-end aggregation, null roles, true effect declarations and
paired comparisons were exercised.

The independent absolute-range allocator was compared with the production
sign-rule implementation. Tests cover PSD/moment matching, recoding invariance,
empty strata/full-rank and rank-fallback, trace formula, data-only API,
shift invariance, deterministic streams, lattice ties and all DGP feasibility
checks. The full-precision validation still has to be run by the user.

The Python environment was Linux/Python3.13.5/NumPy2.3.5/SciPy1.17.0/Numba0.65.1.
The package has Mac-compatible paths and spawn-based parallelism, but was not
executed on macOS. R was unavailable. The R adapter was source-reviewed, with
its matrix algebra checked in Python, but its R self-test was NOT executed.
The original R source snapshots retain their old optional model-based utilities;
the new observed-data wrapper does not call them.

