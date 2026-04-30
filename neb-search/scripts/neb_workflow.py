#!/usr/bin/env python3
"""NEB workflow with pre- and post-run validation.

Three subcommands:
    validate-setup  — check IS/FS pair before running NEB
    validate-path   — check trajectory after NEB has run
    run             — orchestrate the whole thing (validate IS/FS, run via
                      mlip_platform CustomNEB, validate the path)

The validators are ports of utilities developed during prior catalysis runs;
they catch atom-identity swaps, IS/FS geometric mismatches, excessive
diffusion, inter-image jumps, and non-monotonic bond-distance evolution.

Requires:
    pip install ase numpy scipy
    pip install mlip-platform     # for the `run` subcommand only

Usage:
    python neb_workflow.py validate-setup \\
        --initial IS.vasp --final FS.vasp --adsorbate-indices 36 37

    python neb_workflow.py validate-path \\
        --trajectory A2B.traj --adsorbate-indices 36 37 \\
        --monotonic-pairs 36-37

    python neb_workflow.py run \\
        --initial IS.vasp --final FS.vasp \\
        --num-images 7 --adsorbate-indices 36 37 \\
        --mlip uma-s-1p1 --uma-task oc20 --fmax 0.05 \\
        --output-dir neb_v1/
"""

import argparse
import sys
from pathlib import Path

# Make the local utils/ importable when the script is run directly.
sys.path.insert(0, str(Path(__file__).resolve().parent))


def parse_pair(s: str) -> tuple[int, int]:
    """Parse '36-37' → (36, 37)."""
    a, b = s.split("-")
    return (int(a), int(b))


def cmd_validate_setup(args):
    from ase.io import read
    from utils.neb_validation import validate_neb_setup

    initial = read(str(args.initial))
    final = read(str(args.final))
    result = validate_neb_setup(
        initial, final,
        adsorbate_indices=args.adsorbate_indices,
        max_disp=args.max_disp,
        check_crossing=not args.no_crossing_check,
        verbose=not args.quiet,
    )
    return 0 if result["valid"] else 2


def cmd_validate_path(args):
    from utils.neb_validation import validate_neb_path

    monotonic = [parse_pair(p) for p in (args.monotonic_pairs or [])]
    result = validate_neb_path(
        str(args.trajectory),
        adsorbate_indices=args.adsorbate_indices,
        expected_monotonic_pairs=monotonic or None,
        max_jump=args.max_jump,
        verbose=not args.quiet,
    )
    return 0 if result["valid"] else 2


def cmd_run(args):
    from ase.io import read
    from mlip_platform.core.neb import CustomNEB
    from utils.neb_validation import validate_neb_setup, validate_neb_path

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    initial = read(str(args.initial))
    final = read(str(args.final))

    # Pre-flight validation
    if not args.skip_pre_validate:
        print("\n=== Pre-flight: validate-setup ===")
        pre = validate_neb_setup(
            initial, final,
            adsorbate_indices=args.adsorbate_indices,
            max_disp=args.max_disp,
            check_crossing=True,
            verbose=True,
        )
        if not pre["valid"] and not args.force:
            print("\n❌ Pre-flight validation FAILED. "
                  "Re-run with --force to NEB anyway, or fix the IS/FS pair.",
                  file=sys.stderr)
            return 2

    # Image-count heuristic note
    if args.barrier_estimate is not None:
        suggested = max(5, int(round(args.barrier_estimate / 0.1)))
        if args.num_images < suggested:
            print(f"⚠️  num_images={args.num_images} may be coarse for an "
                  f"estimated {args.barrier_estimate} eV barrier "
                  f"(rule of thumb: ≥{suggested}).")

    # Run the NEB
    neb = CustomNEB(
        initial=initial, final=final,
        num_images=args.num_images,
        mlip=args.mlip, uma_task=args.uma_task,
        fmax=args.fmax,
        output_dir=str(output_dir),
    )
    neb.interpolate_idpp()
    neb.run_neb(climb=not args.no_climb, max_steps=args.max_steps)

    # Post-flight validation
    traj_path = output_dir / "A2B.traj"
    monotonic = [parse_pair(p) for p in (args.monotonic_pairs or [])]
    if traj_path.exists() and not args.skip_post_validate:
        print("\n=== Post-flight: validate-path ===")
        post = validate_neb_path(
            str(traj_path),
            adsorbate_indices=args.adsorbate_indices,
            expected_monotonic_pairs=monotonic or None,
            max_jump=args.max_jump,
            verbose=True,
        )
        if not post["valid"]:
            print("\n⚠️  Path validation flagged issues — see above. "
                  "The NEB ran, but the path may be unphysical.",
                  file=sys.stderr)
            return 2
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)

    # ---- validate-setup ----
    s1 = sub.add_parser("validate-setup",
                         help="Pre-flight check on IS/FS pair (no NEB run)")
    s1.add_argument("--initial", type=Path, required=True)
    s1.add_argument("--final", type=Path, required=True)
    s1.add_argument("--adsorbate-indices", type=int, nargs="+", required=True)
    s1.add_argument("--max-disp", type=float, default=3.0,
                    help="Max allowed per-atom displacement (Å, default: 3.0)")
    s1.add_argument("--no-crossing-check", action="store_true")
    s1.add_argument("--quiet", action="store_true")
    s1.set_defaults(func=cmd_validate_setup)

    # ---- validate-path ----
    s2 = sub.add_parser("validate-path",
                         help="Post-flight check on a NEB trajectory")
    s2.add_argument("--trajectory", type=Path, required=True)
    s2.add_argument("--adsorbate-indices", type=int, nargs="+", required=True)
    s2.add_argument("--monotonic-pairs", nargs="*",
                    help="Pairs like '36-37' that should change monotonically along the path")
    s2.add_argument("--max-jump", type=float, default=1.0,
                    help="Max allowed inter-image displacement per atom (Å, default: 1.0)")
    s2.add_argument("--quiet", action="store_true")
    s2.set_defaults(func=cmd_validate_path)

    # ---- run ----
    s3 = sub.add_parser("run", help="Validate, run NEB, and validate path")
    s3.add_argument("--initial", type=Path, required=True)
    s3.add_argument("--final", type=Path, required=True)
    s3.add_argument("--adsorbate-indices", type=int, nargs="+", required=True)
    s3.add_argument("--num-images", type=int, default=7,
                    help="Number of images including endpoints (default: 7)")
    s3.add_argument("--mlip", default="uma-s-1p1")
    s3.add_argument("--uma-task", default="oc20",
                    choices=["omat", "oc20", "omol", "odac"])
    s3.add_argument("--fmax", type=float, default=0.05)
    s3.add_argument("--max-steps", type=int, default=600)
    s3.add_argument("--no-climb", action="store_true")
    s3.add_argument("--max-disp", type=float, default=3.0)
    s3.add_argument("--max-jump", type=float, default=1.0)
    s3.add_argument("--monotonic-pairs", nargs="*")
    s3.add_argument("--barrier-estimate", type=float,
                    help="Approximate barrier in eV (used for an image-count sanity check)")
    s3.add_argument("--output-dir", type=Path, default=Path("."))
    s3.add_argument("--skip-pre-validate", action="store_true")
    s3.add_argument("--skip-post-validate", action="store_true")
    s3.add_argument("--force", action="store_true",
                    help="Run NEB even if pre-flight validation fails")
    s3.set_defaults(func=cmd_run)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main() or 0)
