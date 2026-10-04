#!/usr/bin/env python3
"""Write/check SHA-256 manifest after a locally reviewed repository update."""
from pathlib import Path
import argparse, hashlib
SKIP={'.git','.venv','__pycache__','.pytest_cache','.mypy_cache','.DS_Store'}
SUFFIX={'.pyc','.nbc','.nbi'}
def selected(root:Path):
    for q in sorted(root.rglob('*')):
        r=q.relative_to(root)
        if not q.is_file() or q.is_symlink() or any(x in SKIP for x in r.parts):continue
        if q.suffix in SUFFIX or r.as_posix()=='MANIFEST_SHA256.txt':continue
        if r.parts[:3]==('validation','report','generated'):continue
        yield q

def main():
    p=argparse.ArgumentParser(description=__doc__)
    g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent.parent)
    a=p.parse_args();root=a.root.resolve();manifest=root/'MANIFEST_SHA256.txt'
    if a.write:
        content=''.join(f'{hashlib.sha256(q.read_bytes()).hexdigest()}  {q.relative_to(root).as_posix()}\n' for q in selected(root))
        manifest.write_text(content);print('WROTE',manifest);return
    if not manifest.is_file():raise SystemExit('No manifest. Review changes, then run --write once.')
    failures=[];count=0
    for line in manifest.read_text().splitlines():
        if not line.strip() or line.startswith('#'):continue
        digest,name=line.split(maxsplit=1);name=name.lstrip('* ');q=root/name
        if not q.is_file() or hashlib.sha256(q.read_bytes()).hexdigest()!=digest:failures.append(name)
        count+=1
    if failures:raise SystemExit('Checksum failure(s):\n'+'\n'.join(failures))
    print('PASS',count,'checksums')
if __name__=='__main__':main()
