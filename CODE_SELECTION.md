# Code selection from the development archive

The submitted archive contained many superseded programs. This supplement
retains only the final fixed-score continuous programs, the latest binary
large-design and small-design programs, the latest SWIFT DIRECT-inspired power
program, the latest corrected timing program, and the RT-kernel equivalence
check.

The older continuous `*_vs_randomization_revised.R` files were excluded because
they do not implement the final fixed treatment-score comparator. Earlier
binary, case-study and timing versions were superseded. The old final-table v6
script was also excluded because it expects 36 timing rows, whereas the current
manuscript contains the complete 48-row timing grid. Its role is replaced by
`workflow/build_manuscript_outputs.R` and
`workflow/create_manuscript_tables.R`, which fail if the current output grid is
incomplete.

Production R files are byte-identical to the selected source-archive files.
Checksums are listed in `MANIFEST.csv` and `SHA256SUMS`.
