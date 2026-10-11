# Difference from the supplied review program

1. The default `core` suite (18 scenarios / 1.35 million trials) is removed. The only full default is four large-gap scenarios (220,000 trials) plus a 100,000-trial A13 diagnostic.
2. The Ye-Yi-Shao comparator is not included. No availability-based filtering is performed.
3. The proposed CRT using a Monte Carlo reference-variance denominator is not included. Reference sample variance remains a diagnostic only.
4. The original plug-in CRT and SIGA implementations are unchanged. Oracle-d and existing estimated-noise-subtraction outputs are separate diagnostic reports.
5. Extra calibration and higher-reference sensitivity is limited to A13: three calibration batches and nested 4,999/9,999 references on the same observations. Large-gap scenarios use one independent analysis calibration per n and a separate direction calibration.
6. All results use a new master seed (2026101107). The smoke seed is different (2026101191). This does not reproduce the older seed/direction files.
7. No software installation, GitHub change, manuscript rewrite or original-study overwrite is performed by the simulation program.
8. The larger-effect validation is intended to connect an established theoretical formula to an implementable estimated-effect procedure. It is not an external superiority benchmark and not a search for favorable parameter settings.
