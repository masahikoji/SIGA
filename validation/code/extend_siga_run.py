#!/usr/bin/env python3
"""Safely extend a completed SIGA validation run to more outer trials.

The analysis code, scenarios, seed, reference draws per trial and calibration
replicates are unchanged. Existing completed trial chunks are copied into a new
output directory with a new frozen plan hash; only configuration['outer'] is
increased. The source run is never modified.

Place this file in the SIGA_validation_code directory and run it with the same
Python environment used for the validation package.
"""
from __future__ import annotations

import argparse
import copy
import json
import os
import platform
import shutil
import sys
from pathlib import Path

import numpy as np

from siga_validation.design import stable_id
from siga_validation.experiment import load_plan, source_signature, atomic_json, atomic_npz


def expected_chunks(outer: int, chunk: int):
    for start in range(0, outer, chunk):
        stop = min(start + chunk, outer)
        yield start, stop


def verify_source_complete(root: Path, plan: dict) -> None:
    cfg = plan["configuration"]
    outer = int(cfg["outer"])
    chunk = int(cfg["chunk"])

    # Calibration files must exist for every calibration key.
    for key in sorted({s["calibration_key"] for s in plan["scenarios"]}):
        p = root / "calibration" / f"{key}.npz"
        if not p.exists():
            raise RuntimeError(f"Missing calibration file: {p}")
        with np.load(p, allow_pickle=False) as z:
            if str(z["plan_hash"]) != plan["plan_hash"]:
                raise RuntimeError(f"Calibration plan hash mismatch: {p}")
            if int(z["B0"]) != int(cfg["calibration"]):
                raise RuntimeError(f"Calibration replicate count mismatch: {p}")

    # Every originally planned trial chunk must be present and exact.
    for sc in plan["scenarios"]:
        d = root / "trials" / sc["id"]
        for start, stop in expected_chunks(outer, chunk):
            p = d / f"{start:07d}_{stop:07d}.npz"
            if not p.exists():
                raise RuntimeError(f"Missing completed source chunk: {p}")
            with np.load(p, allow_pickle=False) as z:
                if str(z["plan_hash"]) != plan["plan_hash"]:
                    raise RuntimeError(f"Trial plan hash mismatch: {p}")
                if str(z["scenario_id"]) != sc["id"]:
                    raise RuntimeError(f"Scenario id mismatch: {p}")
                idx = z["trial_index"]
                if not np.array_equal(idx, np.arange(start, stop)):
                    raise RuntimeError(f"Trial indices mismatch: {p}")


def rewrite_npz(src: Path, dst: Path, new_plan_hash: str) -> None:
    with np.load(src, allow_pickle=False) as z:
        arrays = {k: z[k] for k in z.files}
    arrays["plan_hash"] = np.array(new_plan_hash)
    atomic_npz(dst, **arrays)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True, help="Completed source run directory, e.g. results_factorial_pilot")
    ap.add_argument("--out", required=True, help="New output directory for the extended run")
    ap.add_argument("--target-outer", type=int, default=10000)
    args = ap.parse_args()

    source = Path(args.source).expanduser().resolve()
    out = Path(args.out).expanduser().resolve()
    target = int(args.target_outer)

    if not source.exists():
        raise SystemExit(f"Source directory does not exist: {source}")
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"Output directory is not empty: {out}")

    # load_plan verifies that the CURRENT analysis source is exactly the source
    # that created the frozen run. Adding this helper file does not change the
    # package source signature because source_signature hashes run.py and
    # siga_validation/*.py only.
    plan = load_plan(source)
    old_outer = int(plan["configuration"]["outer"])
    if target <= old_outer:
        raise SystemExit(f"target-outer must exceed the source outer count ({old_outer}).")

    verify_source_complete(source, plan)

    # Only the number of outer trials is changed. Reference draws, calibration
    # size, scenario definitions, seed and analysis source remain identical.
    new_plan = copy.deepcopy(plan)
    old_hash = new_plan.pop("plan_hash")
    new_plan["profile"] = f"{plan.get('profile', 'run')}_extended_{target}"
    new_plan["configuration"]["outer"] = target
    new_plan["source_signature"] = source_signature()
    new_plan["extension"] = {
        "source_plan_hash": old_hash,
        "source_outer": old_outer,
        "target_outer": target,
        "reused_trial_indices": [0, old_outer - 1],
        "changed_quantities": ["configuration.outer"],
        "unchanged_reference_draws": int(plan["configuration"]["reference"]),
        "unchanged_calibration_replicates": int(plan["configuration"]["calibration"]),
        "note": "Existing completed trials are reused; no method, scenario, seed, reference-draw or calibration setting is changed."
    }
    new_plan["plan_hash"] = stable_id(new_plan)
    new_hash = new_plan["plan_hash"]

    out.mkdir(parents=True, exist_ok=True)
    atomic_json(out / "plan.json", new_plan)
    if (source / "plan.csv").exists():
        shutil.copy2(source / "plan.csv", out / "plan.csv")
    if (source / "environment.json").exists():
        shutil.copy2(source / "environment.json", out / "source_environment.json")
    atomic_json(out / "extension_environment.json", {
        "python": sys.version,
        "platform": platform.platform(),
        "source_signature": source_signature(),
        "source_directory": str(source),
        "source_plan_hash": old_hash,
        "new_plan_hash": new_hash,
    })

    # Reuse the exact allocation-only calibration values, rewriting only the
    # plan hash so the normal runner can verify the derived frozen plan.
    for key in sorted({s["calibration_key"] for s in plan["scenarios"]}):
        rewrite_npz(source / "calibration" / f"{key}.npz",
                    out / "calibration" / f"{key}.npz", new_hash)

    # Reuse all completed source trial chunks 0,...,old_outer-1 unchanged apart
    # from the plan-hash metadata. The source directory is never altered.
    chunk = int(plan["configuration"]["chunk"])
    copied_chunks = 0
    copied_trials = 0
    for sc in plan["scenarios"]:
        for start, stop in expected_chunks(old_outer, chunk):
            src = source / "trials" / sc["id"] / f"{start:07d}_{stop:07d}.npz"
            dst = out / "trials" / sc["id"] / src.name
            rewrite_npz(src, dst, new_hash)
            copied_chunks += 1
            copied_trials += stop - start

    atomic_json(out / "EXTENSION_PROVENANCE.json", {
        "source_directory": str(source),
        "source_plan_hash": old_hash,
        "new_plan_hash": new_hash,
        "source_outer_per_scenario": old_outer,
        "target_outer_per_scenario": target,
        "scenarios": len(plan["scenarios"]),
        "reused_trial_chunks": copied_chunks,
        "reused_trial_records": copied_trials,
        "reference_draws_per_trial": int(plan["configuration"]["reference"]),
        "calibration_replicates_per_design": int(plan["configuration"]["calibration"]),
        "analysis_source_signature": source_signature(),
    })

    print(f"Verified complete source run: {source}")
    print(f"Created extended frozen run: {out}")
    print(f"Scenarios: {len(plan['scenarios'])}")
    print(f"Reused outer trials per scenario: {old_outer}")
    print(f"Target outer trials per scenario: {target}")
    print(f"Remaining outer trials per scenario: {target-old_outer}")
    print(f"Reference draws/trial unchanged: {plan['configuration']['reference']}")
    print(f"Calibration replicates/design unchanged: {plan['configuration']['calibration']}")
    print("New plan hash:", new_hash)
    print("Next: python run.py calibrate --out <NEW_OUT> --workers 4")
    print("Then: python run.py run --out <NEW_OUT> --workers 4")
    print("Then: python run.py aggregate --out <NEW_OUT>")


if __name__ == "__main__":
    main()
