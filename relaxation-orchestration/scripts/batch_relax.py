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
import json
import sys
import time
from pathlib import Path

VALID_UMA_TASKS = ("omat", "oc20", "omol", "odac")
RUN_META_FILENAME = "compcat_run.json"


def _import_mlip_platform():
    """Lazy-import mlip_platform with a friendly error if it's missing.

    Returns (setup_calculator, run_optimization, detect_mlip).
    """
    try:
        from mlip_platform.cli.utils import setup_calculator, detect_mlip
        from mlip_platform.core.optimize import run_optimization
        return setup_calculator, run_optimization, detect_mlip
    except ImportError as e:
        print(
            "ERROR: mlip-platform is not available in this environment.\n"
            "  Install it (e.g. `pip install mlip-platform`) or activate the env\n"
            "  that ships it. The --dry-run flag does not require mlip-platform.\n"
            f"  Original ImportError: {e}",
            file=sys.stderr,
        )
        sys.exit(2)


def resolve_mlip(name: str) -> str:
    """Resolve an --mlip argument to a concrete model name.

    Pass-through for explicit names. For 'auto', delegate to mlip_platform's
    detect_mlip() so the choice tracks the platform's logic instead of being
    pinned in this skill.
    """
    if name != "auto":
        return name
    _, _, detect_mlip = _import_mlip_platform()
    return detect_mlip()


def find_structures(tree: Path, pattern: str) -> list[Path]:
    """Return sorted list of structure files matching `pattern` under `tree`."""
    if not tree.exists():
        print(f"ERROR: tree not found: {tree}", file=sys.stderr)
        sys.exit(1)
    return sorted(tree.rglob(pattern))


def is_already_converged(structure: Path, fmax: float | None = None) -> bool:
    """Decide if a prior run satisfied the requested convergence.

    Behavioral check (not string-grep): the run is considered converged if
      (a) opt_final.vasp + opt_convergence.csv both exist, AND
      (b) the last row of opt_convergence.csv has fmax <= the requested
          threshold (or compcat_run.json's recorded fmax if no `fmax` is
          passed in here).
    """
    out_dir = structure.parent
    final_vasp = out_dir / "opt_final.vasp"
    conv_csv = out_dir / "opt_convergence.csv"
    if not (final_vasp.exists() and conv_csv.exists()):
        return False

    threshold = fmax
    meta = out_dir / RUN_META_FILENAME
    if threshold is None and meta.exists():
        try:
            threshold = float(json.loads(meta.read_text()).get("fmax", 0))
        except (OSError, ValueError, json.JSONDecodeError):
            threshold = None
    if threshold is None:
        return False

    try:
        with open(conv_csv) as f:
            rows = list(csv.DictReader(f))
        if not rows:
            return False
        # mlip_platform writes fmax(eV/Å) — pick whatever fmax-ish column exists
        last = rows[-1]
        for key in ("fmax(eV/Å)", "fmax(eV/A)", "fmax", "Fmax"):
            if key in last:
                return float(last[key]) <= threshold
        return False
    except (OSError, ValueError):
        return False


def has_vacuum(atoms, axis_threshold: float = 5.0) -> list[bool]:
    """Return a bool per cell axis: True if there's >= axis_threshold Å of
    empty space in that direction. A bulk has all False; a slab has exactly
    one True (typically c); a molecule in a vacuum box has all True."""
    import numpy as np

    pos = atoms.get_positions()
    cell = atoms.get_cell()
    flags = []
    for i in range(3):
        axis = cell[i]
        norm = np.linalg.norm(axis)
        if norm < 1e-6:
            flags.append(True)  # degenerate cell axis = effectively unbounded
            continue
        unit = axis / norm
        proj = pos @ unit
        flags.append((norm - (proj.max() - proj.min())) >= axis_threshold)
    return flags


def classify_structure(atoms) -> str:
    """Return one of: 'bulk', 'slab', 'molecule', 'ambiguous'."""
    flags = has_vacuum(atoms)
    n_vac = sum(flags)
    if n_vac == 0:
        return "bulk"
    if n_vac == 1:
        return "slab"
    if n_vac >= 2 and len(atoms) <= 30:
        return "molecule"
    return "ambiguous"


def infer_uma_task(structures: list[Path]) -> str | None:
    """Return a UMA task head only when every input is unambiguous and the
    task is one we can pick safely (omat for bulks, omol for small isolated
    molecules). Otherwise return None — caller refuses with a help message."""
    from ase.io import read

    classes = set()
    for s in structures:
        try:
            atoms = read(str(s))
        except Exception:
            return None
        classes.add(classify_structure(atoms))

    if classes == {"bulk"}:
        return "omat"
    if classes == {"molecule"}:
        return "omol"
    return None


def write_run_metadata(out_dir: Path, *, mlip: str, uma_task: str,
                        fmax: float, optimizer: str) -> None:
    """Persist the run settings next to opt_final.vasp so downstream
    skills (adsorption-energy) can verify task-head consistency."""
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / RUN_META_FILENAME).write_text(json.dumps({
        "mlip": mlip,
        "uma_task": uma_task,
        "fmax": fmax,
        "optimizer": optimizer,
    }, indent=2))


def relax_one(structure: Path, *, mlip: str, uma_task: str, optimizer: str,
              fmax: float, max_steps: int, verbose: bool) -> dict:
    """Run a single optimization. Returns a result row.

    Imports happen inside the function so unit tests can monkeypatch the
    module without the heavy mlip_platform import on collection, and so
    `--dry-run` works in environments without mlip-platform installed.
    """
    from ase.io import read

    setup_calculator, run_optimization, _ = _import_mlip_platform()

    atoms = read(str(structure))
    atoms = setup_calculator(atoms, mlip=mlip, uma_task=uma_task)
    t0 = time.time()
    converged = run_optimization(
        atoms=atoms, optimizer=optimizer, fmax=fmax, max_steps=max_steps,
        output_dir=structure.parent, model_name=mlip, verbose=verbose,
    )
    write_run_metadata(structure.parent, mlip=mlip, uma_task=uma_task,
                        fmax=fmax, optimizer=optimizer)
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
    parser.add_argument("--mlip", default="auto",
                        help="MLIP model: 'auto' (let mlip_platform.detect_mlip pick), "
                             "or an explicit name like 'uma-s-1p1', 'uma-s-1p2', "
                             "'mace', '7net-mf-ompa'. Default: auto.")
    parser.add_argument("--uma-task", default="auto",
                        choices=["auto", *VALID_UMA_TASKS],
                        help="UMA task head. 'auto' picks omat for bulks "
                             "and omol for isolated molecules; for slabs / "
                             "adsorbate systems you must pass an explicit "
                             "task (oc20 metal surfaces, oc22 oxides, oc25 "
                             "solid–liquid interfaces — see SKILL.md).")
    parser.add_argument("--optimizer", default="bfgs",
                        help="Optimizer: bfgs, lbfgs, fire, ... Default 'bfgs' "
                             "matches mlip_platform's choice. Switch to 'fire' "
                             "if BFGS fails to converge on noisy MLIP forces.")
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

    if not args.dry_run:
        # Resolve --mlip auto early so the choice is logged and persisted into
        # compcat_run.json. Skip in dry-run to keep the dry path mlip_platform-free.
        try:
            args.mlip = resolve_mlip(args.mlip)
            if args.mlip != "auto":
                print(f"MLIP backend: {args.mlip}")
        except SystemExit:
            raise

    if args.uma_task == "auto":
        inferred = infer_uma_task(structures)
        if inferred is None:
            print(
                "ERROR: could not safely infer --uma-task for the given "
                "structures. The 'auto' setting only handles all-bulks "
                "(omat) and all-molecules (omol). For slabs / adsorbate "
                "systems pass an explicit --uma-task. See the SKILL.md "
                "task-head decision table:\n"
                "  oc20 = metal surface chemistry (Pt, Cu, Ni, ...)\n"
                "  oc22 = oxide surfaces (CoOOH, Co3O4, TiO2, ...) "
                "[requires mlip_platform support]\n"
                "  oc25 = solid–liquid interfaces with explicit solvent "
                "[requires mlip_platform support]\n"
                "  omat = bulk crystals\n"
                "  omol = isolated molecules in vacuum",
                file=sys.stderr)
            sys.exit(2)
        print(f"Auto-inferred --uma-task = {inferred}")
        args.uma_task = inferred

    pending, skipped = [], []
    for s in structures:
        if args.resume and is_already_converged(s, fmax=args.fmax):
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
