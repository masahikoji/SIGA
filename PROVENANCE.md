# Production provenance

This release distinguishes three kinds of material:

1. `production/`: source snapshots used for the reported simulations;
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
