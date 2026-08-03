# Release checklist

Before creating the public GitHub release or uploading Supplementary Software:

```bash
bash scripts/00_verify_manifest.sh
Rscript scripts/06_verify_production_results.R
Rscript production/common/check_pair_path_engine.R
Rscript production/common/00_verify_rt_kernel_equivalence.R
Rscript reference_implementation/scripts/99_smoke_test.R
Rscript scripts/07_session_info.R > session_info.txt
```

The package assembly environment did not contain R, so the R commands above
must be run in an author environment. The production workflows themselves had
already generated the included 100,000-replicate aggregate outputs; this
checklist validates the assembled release tree and portable smoke tests.

After adding `session_info.txt`, regenerate `MANIFEST_SHA256.txt`, commit the
snapshot, and create tag/release `v1.0.0`.

The included `LICENSE` is restrictive. Replace it before public release if a
permissive open-source license is intended.
