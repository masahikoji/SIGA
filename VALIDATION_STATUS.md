# Validation status

- The manuscript and Supplementary Material TeX sources were compiled successfully with resolved citations and cross-references.
- Shell workflow scripts passed `bash -n` static syntax checks during package assembly.
- File inventories and SHA-256 checksums were generated after assembly.
- The current assembly environment did not provide `Rscript`; therefore, the R smoke tests and manuscript-scale simulations were not re-executed during packaging.
- Before public release and journal upload, run `bash workflow/run_smoke_tests.sh` and `cd theory_extension && bash workflow/run_smoke_test.sh` in the intended R environment.
- Full manuscript results require the production profiles and substantial computation; reduced smoke tests do not reproduce manuscript estimates.
