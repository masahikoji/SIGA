# Reproducibility status: 2026-10-04

## Extended evaluation checked

The supplied package has actual plan.json and EXTENSION_PROVENANCE.json files
for factorial (40), controls (16) and broad/confirm (224) suites. Every scenario
has 10,000 complete trials. These are 2,800,000 trials, including the 560,000 pilot
trials. Two score analyses do not double the independent trial count.

All source signatures equal:
`c70e90a5f7f6060a1daa5fe3cf62ac219a29be6a3e7d5f0a795031a69c867182`.

Checks verify saved plan hashes, exact linkage to the retained simulation source,
reverse extension to the pilot hash, common seed 2026100319, 999 reference draws,
20,000 allocation replicates, 13 analysis columns per score, and all exported
scenario fields. Pilot IDs match the saved plans, including the binary scenarios.
Differences between 10,000- and 2,000-trial rejection/discordance counts are feasible
for 8,000 added trials. Source/plan checks do not authenticate absent raw checkpoints.

Rejection counts, binomial MCSEs, Wilson intervals, paired differences/discordances,
and the reference-variance decomposition have been checked. No failures, zero
reference variances, variance-floor or rank-fallback activations occur in the
provided summaries. Table reproduction does not discard unfavorable rows.

## What is not provided or not newly executed

The actual runtime environment.json, calibration NPZ arrays and trial checkpoints
were not supplied in the current package. The code is runnable from saved plans,
but the original 2.8-million-trial execution has not been rerun here. The R adapter
has not been execution-tested here. Older theoretical verification files retain
their original verification status; document compilation is not a new proof audit.

The earlier R3 strong summary is still the explicitly labelled reconstructed CSV
under additional_null_study/. It must not be relabelled as the original aggregator
output. Its original execution provenance is separate from this complete 10,000-
trial aggregate record. Existing limitations are not repaired by new validation.

## Statistical interpretation

The reconstruction and fixed scenarios preceded the 2,000-trial pilot. Increasing
the number of trials did not change formulas or settings. Both methods use the
same observed outcomes and estimated stratum effects; no true outcome parameter
enters the primary analysis. Known design-stage profile probabilities are shared.
No independent confirmation claim is made for totals that retain the pilot.

At a true rejection rate of 5%, 10,000 trials give about 0.22 percentage points
of MCSE. Paired method comparisons use their discordance-based MCSE. These are
scenario-specific and conditional on common estimated allocation covariances,
not simultaneous accuracy guarantees. Inner Monte Carlo error also remains.
No correlated-profile run or new reconstruction-time benchmark is reported.

## Public repository

Read-only check on 2026-10-04 found public main at
8c3560802f67b72aaf74fa7140c8b453a8dfdd24. This overlay updates descriptions and adds
local files; no push, release or DOI creation was performed. Retain v1.0.0 and
finalise a new versioned archive only after review. Before acceptance JRSSB asks
for public reproducibility code with URL and DOI in the manuscript:
https://academic.oup.com/jrsssb/pages/general-instructions
