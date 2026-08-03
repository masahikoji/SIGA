# Portable reference implementation

This directory contains a consolidated, base-R implementation of the methods
for reading, smoke testing, and small independent reruns. It is useful for
understanding the algorithms, but it is **not** the authoritative source of
the numerical values reported in the manuscript.

The exact scripts used to generate the reported results are preserved under
`../production/`, and their machine-generated aggregate outputs are under
`../data/production_results/raw/`.

Run from the repository root:

```bash
Rscript reference_implementation/scripts/99_smoke_test.R
```

The manuscript-scale production workflows should be run using the scripts and
environment variables documented in `../production/README.md`.
