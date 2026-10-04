#!/usr/bin/env python3
"""Prepare (but do not execute) a fresh rerun from the supplied 10,000-trial plan.

This keeps saved binary risks/scenario IDs, the master seed, inner draws and
calibration counts. It creates a new output directory and explicit rerun provenance.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
from pathlib import Path
import shlex
import sys


def stable_id(value: dict) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True,
                                    separators=(',', ':')).encode()).hexdigest()


def signature(code: Path) -> str:
    h = hashlib.sha256()
    for path in sorted((code/'siga_validation').glob('*.py')) + [code/'run.py']:
        h.update(path.relative_to(code).as_posix().encode())
        h.update(path.read_bytes())
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite', choices=['factorial', 'controls', 'confirm'], required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    validation = Path(__file__).resolve().parents[1]
    code = validation/'code'
    source = validation/'results_10k'/('results_'+args.suite+'_10k')/'plan.json'
    original = json.loads(source.read_text())
    check = dict(original)
    original_hash = check.pop('plan_hash')
    if stable_id(check) != original_hash:
        raise ValueError('Supplied plan hash does not match its contents.')
    if original['source_signature'] != signature(code):
        raise ValueError('Simulation code does not match the original plan signature.')
    if original['configuration'] != {'outer':10000,'reference':999,'calibration':20000,'chunk':100}:
        raise ValueError('Unexpected repetition configuration.')
    out = args.out.expanduser().resolve()
    if out.exists():
        raise ValueError('Destination already exists; use a new directory: '+str(out))
    plan = copy.deepcopy(original)
    plan.pop('extension', None)
    plan.pop('plan_hash')
    plan['profile'] = 'frozen_10000_rerun'
    plan['rerun_of'] = original_hash
    plan['plan_hash'] = stable_id(plan)
    provenance = {'source_plan_hash': original_hash, 'rerun_plan_hash': plan['plan_hash'],
                  'reused_checkpoint_trials': 0,
                  'settings_changed': False, 'scenario_definitions_preserved': True,
                  'source_signature': original['source_signature'],
                  'note': 'Fresh rerun with original random streams; not independent evidence.'}
    out.mkdir(parents=True)
    (out/'plan.json').write_text(json.dumps(plan, indent=2, sort_keys=True, allow_nan=False)+'\n')
    (out/'RERUN_PROVENANCE.json').write_text(json.dumps(provenance, indent=2)+'\n')
    print('Prepared a new plan only. No simulations or network writes have started.')
    print('Retained 10,000 outcome trials / 999 reference draws / 20,000 allocation replicates.')
    for command in ['calibrate', 'run', 'aggregate']:
        args2 = [sys.executable, str(code/'run.py'), command, '--out', str(out)]
        if command != 'aggregate':
            args2 += ['--workers', '4']
        print(shlex.join(args2))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError) as exc:
        sys.exit('ERROR: '+str(exc))
