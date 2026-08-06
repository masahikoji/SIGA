# Production provenance

This release distinguishes three kinds of material:

1. `production/`: frozen production drivers used to reproduce the reported simulation configurations and aggregate outputs;
2. `data/production_results/raw/`: machine-generated aggregate outputs from
   those runs; and
3. `reference_implementation/`: a consolidated portable implementation for
   inspection and small independent checks.

The production aggregate files, not LaTeX transcriptions, are the source of
all checks in `scripts/06_verify_production_results.R`. `PROVENANCE.csv`
records SHA-256 hashes for every included production source and aggregate
file. The raw replicate-level shard files and calibration caches are omitted
because they are very large; the source scripts retain checkpointing and
aggregation logic.

The SWIFT DIRECT section is a trial-inspired simulation using published
aggregate planning characteristics. It is not a participant-level reanalysis
and does not reconstruct the completed trial's deterministic minimization
path.

The targeted pair-path sensitivity release uses the 12-scenario supplemental configuration (base seed 20260805, 120 shards). The included production driver is a semantically exact deparse of the 105,077-byte source used for that run; the original source SHA-256 was `35e0adb0da0b7b584fb093cee5f4869c2ad7e79242cd7b02c0c0ea2e19988ca6`. The only post-run change in the driver is the `unname()` correction in the Wilson rejection-summary helper, which affects aggregate column names but not simulated trials. The publication-facing table excludes Oracle-N and method-minus-reference columns; those quantities remain only in the raw aggregate for diagnostic provenance.
