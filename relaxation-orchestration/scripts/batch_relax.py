#!/usr/bin/env python3
"""Batch geometry optimization wrapper around mlip_platform.

Scans a directory tree for input structures and runs MLIP relaxation on
each one in turn, producing the standard mlip_platform output set
(opt.traj, opt_final.vasp, opt_convergence.csv) per structure. Designed
for catalysis workflows: a clean slab plus N slab+adsorbate inputs from
the adsorbate-placement skill, all relaxed with consistent settings.

Requires:
    pip install mlip-platform ase

Usage:
    # Single structure
    python batch_relax.py --structure slab.vasp \\
        --mlip uma-s-1p1 --uma-task oc20 --fmax 0.03

    # All input.vasp files under a tree (e.g. adsorbate-placement output)
    python batch_relax.py --tree placements/ --pattern 'input.vasp' \\
        --mlip uma-s-1p1 --uma-task oc20 --fmax 0.03

    # Resume — skip directories that already have a converged opt run
    python batch_relax.py --tree placements/ --resume

    # Dry run — print what would be relaxed, without invoking the MLIP
    python batch_relax.py --tree placements/ --dry-run
"""

import argparse
import csv
import sys
import time
from pathlib import Path


def find_structures(tree: Path, pattern: str) -> list[Path]:
    """Return sorted list of structure files matching `pattern` under `tree`."""
    if not tree.exists():
        print(f"ERROR: tree not found: {tree}", file=sys.stderr)
        sys.exit(1)
    return sorted(tree.rglob(pattern))


def is_already_converged(structure: Path) -> bool:
    """Heuristic check: a previous run produced opt_final.vasp + a CSV with
    a final fmax row below the threshold recorded in opt_params.txt."""
    out_dir = structure.parent
    final_vasp = out_dir / "opt_final.vasp"
    conv_csv = out_dir / "opt_convergence.csv"
    params = out_dir / "opt_params.txt"
    if not (final_vasp.exists() and conv_csv.exists() and params.exists()):
        return False
    try:
        text = params.read_text()
        return "Converged:         True" in text
    except OSError:
        return False


def relax_one(structure: Path, *, mlip: str, uma_task: str, optimizer: str,
              fmax: float, max_steps: int, verbose: bool) -> dict:
    """Run a single optimization. Returns a result row.

    Imports happen inside the function so unit tests can monkeypatch the
    module without the heavy mlip_platform import on collection.
    """
    from ase.io import read
    from mlip_platform.cli.utils import setup_calculator
    from mlip_platform.core.optimize import run_optimization

    atoms = read(str(structure))
    atoms = setup_calculator(atoms, mlip=mlip, uma_task=uma_task)
    t0 = time.time()
    converged = run_optimization(
        atoms=atoms, optimizer=optimizer, fmax=fmax, max_steps=max_steps,
        output_dir=structure.parent, model_name=mlip, verbose=verbose,
    )
    return {
        "structure": str(structure),
        "converged": bool(converged),
        "duration_s": round(time.time() - t0, 2),
        "final_vasp": str(structure.parent / "opt_final.vasp"),
        "trajectory": str(structure.parent / "opt.traj"),
    }


def write_summary(rows: list[dict], path: Path) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Batch geometry optimization with mlip_platform.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument("--structure", type=Path, help="Single input structure (.vasp)")
    grp.add_argument("--tree", type=Path, help="Directory tree to search for inputs")
    parser.add_argument("--pattern", default="input.vasp",
                        help="Filename pattern under --tree (default: input.vasp)")
    parser.add_argument("--mlip", default="uma-s-1p1",
                        help="MLIP model (default: uma-s-1p1; or 'auto')")
    parser.add_argument("--uma-task", default="oc20",
                        choices=["omat", "oc20", "omol", "odac"],
                        help="UMA task head (default: oc20 for catalysis surfaces)")
    parser.add_argument("--optimizer", default="fire",
                        help="Optimizer: fire, bfgs, lbfgs, ... (default: fire)")
    parser.add_argument("--fmax", type=float, default=0.03,
                        help="Force convergence threshold eV/Å (default: 0.03)")
    parser.add_argument("--max-steps", type=int, default=300,
                        help="Maximum optimization steps (default: 300)")
    parser.add_argument("--resume", action="store_true",
                        help="Skip structures whose parent dir already has a converged run")
    parser.add_argument("--dry-run", action="store_true",
                        help="List the structures that would be relaxed and exit")
    parser.add_argument("--summary", default="batch_relax.csv",
                        help="Summary CSV path (default: batch_relax.csv)")
    parser.add_argument("--verbose", action="store_true", default=True)
    args = parser.parse_args()

    if args.structure:
        structures = [args.structure]
    else:
        structures = find_structures(args.tree, args.pattern)
        if not structures:
            print(f"ERROR: no files matched '{args.pattern}' under {args.tree}",
                  file=sys.stderr)
            sys.exit(1)

    pending, skipped = [], []
    for s in structures:
        if args.resume and is_already_converged(s):
            skipped.append(s)
        else:
            pending.append(s)

    print(f"Structures total: {len(structures)} | pending: {len(pending)} | "
          f"skipped (resume): {len(skipped)}")
    for s in pending:
        print(f"  → {s}")

    if args.dry_run:
        return

    rows = []
    for s in pending:
        print(f"\n=== {s} ===")
        try:
            row = relax_one(
                s, mlip=args.mlip, uma_task=args.uma_task,
                optimizer=args.optimizer, fmax=args.fmax,
                max_steps=args.max_steps, verbose=args.verbose,
            )
        except Exception as e:
            row = {"structure": str(s), "converged": False,
                   "duration_s": 0.0, "final_vasp": "", "trajectory": "",
                   "error": str(e)}
            print(f"  ERROR: {e}", file=sys.stderr)
        rows.append(row)

    summary_path = (args.tree if args.tree else args.structure.parent) / args.summary \
        if not Path(args.summary).is_absolute() else Path(args.summary)
    write_summary(rows, summary_path)
    print(f"\nSummary: {summary_path}")
    failed = [r for r in rows if not r.get("converged")]
    if failed:
        print(f"⚠️  {len(failed)}/{len(rows)} did not converge")
        sys.exit(2)


if __name__ == "__main__":
    main()
