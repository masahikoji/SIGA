#!/usr/bin/env python3
from pathlib import Path
import subprocess, sys
root=Path(__file__).resolve().parents[1]
for script in ['audit/rebuild_tables.py','audit/verify_complete_benchmark.py','scripts/check_n3.py']:
    subprocess.run([sys.executable,str(root/script)],cwd=root,check=True)
print('Submission tables, timing calculations and exact three-participant identities verified.')
