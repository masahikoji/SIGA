#!/usr/bin/env python3
"""Prepare a NEW run from the exact reported scenario arrays; never alter an old run.

An explicit --allow-environment-change is required on a different platform or
software version. This preserves scenario IDs, arrays and seeds but does not
promise bitwise-identical floating-point output across platforms.
"""
from pathlib import Path
import argparse, copy, json, sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from siga_validation.design import stable_id
from siga_validation.experiment import numerical_environment, source_signature
from siga_validation.resimulation import freeze

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',choices=['primary','additional'],required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--allow-environment-change',action='store_true')
    args=p.parse_args()
    original=json.loads((ROOT/'plans'/f'{args.study}.json').read_text())
    plan=copy.deepcopy(original);given=plan.pop('plan_hash')
    if stable_id(plan)!=given: raise RuntimeError('The archived plan checksum does not match.')
    if source_signature()!=plan['source_signature']: raise RuntimeError('The archived simulation code has changed.')
    expected=plan['numerical_environment'];actual=numerical_environment()
    if actual!=expected and not args.allow_environment_change:
        raise RuntimeError('The numerical environment differs from the reported run. Restore it or explicitly use --allow-environment-change for an independent rerun. Existing checkpoints must not be reused.')
    plan['replication_of_plan_hash']=given
    plan['replication_environment_changed']=actual!=expected
    frozen=freeze(args.out.expanduser(),plan)
    print('Prepared new run:',args.out.expanduser())
    print('Original scenario arrays and IDs retained; original plan hash:',given)
    print('New run hash:',frozen['plan_hash'])
    print('Trials planned:',sum(x['outer'] for x in frozen['scenarios']))
if __name__=='__main__': main()
