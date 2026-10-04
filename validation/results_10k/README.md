# Extended 10,000-trial results

The three directories contain the supplied plans, extension records and complete
aggregate files. None of these input files has been edited. There are 280 scenarios
and 2,800,000 outcome trials, with two score analyses on the same trial.
The first 2,000 trials per scenario are included, not an independent comparison.
Archived pilot summaries are in ../results/ and must not be pooled again with
these extended summaries. Reference draws remain 999 and allocation replicates
remain 20,000 per design. Use ../report/reproduce_10k_report.py for checks/tables.

plan.json and EXTENSION_PROVENANCE.json were supplied; environment.json,
actual calibration NPZ arrays and trial-level checkpoints were not. This package
supports verification of the reported arithmetic and linkage to the unchanged
source, and full reruns from the saved scenario definitions. It does not provide
a raw-checkpoint audit of the original execution. Some extension records contain
local machine paths; they are retained exactly for provenance. Review those paths
before public release. No credentials or participant data are included.
