# Change to full simulation counts (11 October 2026)

- A13: independent outcome trials changed from 2,000 to **100,000**.
- Large-gap null settings remain 100,000 each; alternatives remain 10,000 each. Total: **320,000** independent outcome trials.
- All primary analyses use **4,999** reference assignments per outcome trial.
- The A13 paired sensitivity audit remains restricted to the first **2,000** trials: three independent covariance calibrations and 4,999 versus 9,999 references. This subset is NOT an independent replication.
- Direction and covariance calibration counts remain 200,000 and 100,000, respectively.
- The default output paths were changed to distinct `*_320k_20261011` folders to prevent accidental reuse of the old frozen 222k plan.
- No inference code or upstream vendor modules were modified. The 320,000-trial run has not been performed.
