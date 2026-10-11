from __future__ import annotations
import contextlib
import csv
import hashlib
import json
import os
import platform
import sys
import time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
# A snapshot is used, not whichever similarly named module happens to be installed.
sys.path.insert(0, str(ROOT / 'vendor'))


def stable_hash(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def file_hash(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def rng_for(master_seed, *keys):
    raw = json.dumps([int(master_seed), *keys], separators=(',', ':')).encode()
    entropy = np.frombuffer(hashlib.sha256(raw).digest(), dtype='<u4').tolist()
    return np.random.Generator(np.random.PCG64DXSM(np.random.SeedSequence(entropy)))


def atomic_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f'.tmp.{os.getpid()}')
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + '\n', encoding='utf-8')
    os.replace(tmp, path)


def atomic_npz(path, **arrays):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f'.tmp.{os.getpid()}')
    with open(tmp, 'wb') as f:
        np.savez_compressed(f, **arrays)
    os.replace(tmp, path)


def write_csv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if not rows:
        path.write_text('', encoding='utf-8')
        return
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def environment():
    import scipy
    import numba
    return dict(python=platform.python_version(), numpy=np.__version__,
                scipy=scipy.__version__, numba=numba.__version__,
                platform=platform.platform(), machine=platform.machine())


def verify_vendor():
    record = json.loads((ROOT / 'UPSTREAM.json').read_text())
    for rel, entry in record['files'].items():
        if file_hash(ROOT / rel) != entry['sha256']:
            raise RuntimeError(f'Vendored upstream file was changed: {rel}')
    return record


def source_signature():
    paths = [ROOT / 'run_focus.py', ROOT / 'UPSTREAM.json']
    paths += sorted((ROOT / 'revision').glob('*.py'))
    paths += sorted((ROOT / 'vendor/siga_validation').glob('*.py'))
    return stable_hash({str(p.relative_to(ROOT)): file_hash(p) for p in paths})


def load_plan(root, check_environment=True):
    plan = json.loads((Path(root) / 'plan.json').read_text())
    expected = plan['plan_hash']
    actual = stable_hash({k: v for k, v in plan.items() if k != 'plan_hash'})
    if actual != expected:
        raise RuntimeError('plan.json was changed after freezing. Start a new output directory.')
    verify_vendor()
    if source_signature() != plan['source_signature']:
        raise RuntimeError('Executable source changed. Do not mix versions in one run directory.')
    if check_environment:
        now = environment()
        for key in ('python', 'numpy', 'scipy', 'numba', 'machine'):
            if now[key] != plan['environment'][key]:
                raise RuntimeError(f'Environment mismatch for {key}: {now[key]} != {plan["environment"][key]}. Start a new run.')
    return plan


@contextlib.contextmanager
def run_lock(root, operation):
    """Prevent two command invocations from writing the same run directory."""
    path = Path(root) / '.revision.lock'
    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as e:
        raise RuntimeError(f'{path} already exists. Check for an active process; remove ONLY a confirmed stale lock.') from e
    try:
        with os.fdopen(fd, 'w') as f:
            json.dump(dict(pid=os.getpid(), operation=operation, time=time.time(), host=platform.node()), f)
        yield
    finally:
        path.unlink(missing_ok=True)


def se(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.std(values, ddof=1) / np.sqrt(len(values))) if len(values) > 1 else float('nan')
