# Release checklist

Before creating the public GitHub release:

```bash
bash scripts/00_verify_manifest.sh
Rscript scripts/06_verify_production_results.R
Rscript production/common/check_pair_path_engine.R
Rscript production/common/00_verify_rt_kernel_equivalence.R
Rscript reference_implementation/scripts/99_smoke_test.R
Rscript scripts/07_session_info.R > session_info.txt
```

All checks should pass before the release is tagged. After adding
`session_info.txt`, regenerate `MANIFEST_SHA256.txt`, commit the final snapshot,
and create tag/release `v1.0.0`.

The included `LICENSE` is restrictive. Replace it before public release if a
permissive open-source license is intended.
