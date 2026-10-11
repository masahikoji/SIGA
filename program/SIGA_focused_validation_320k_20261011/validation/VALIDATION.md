# Validation actually executed

- 20 focused unit tests passed (see the final smoke logs, which start with the unit-test run).
- Original allocation.py and core.py match both the donor manifest and the earlier code-only archive byte for byte.
- Full init -> calibrate -> run -> report -> export pipeline completed at smoke budgets: 36 unique trials covering all five settings, both scores, both covariance constructions and all A13 audit branches.
- Repeating with one versus two worker processes produced identical values in 90 non-timing arrays across 10 checkpoints.
- Running the same completed command again reused the checkpoints and retained completed status.
- A separate 15-trial smoke run used the production reference counts (4,999 and 9,999) and exercised all five settings. It used reduced calibration counts and is not a performance study.
- The production plan count is tested: 222,000 unique outcome trials, seven calibration jobs and 900,000 total allocation triples, including independent direction selection and A13 sensitivity batches.

Not performed: the 222,000-trial production study, the original 3.04-million-trial study, or direct execution on an Apple M3 Ultra. The smoke results are not evidence of type-I error or power. The prior donor program's reported 12,000-trial pilot is not claimed as newly reproduced here.
