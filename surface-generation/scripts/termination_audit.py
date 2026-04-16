#!/usr/bin/env python3
"""Audit terminations for a given surface — verify that pymatgen finds all possible cuts.

The core idea: the number of unique terminations equals the number of distinct
atomic planes in the oriented unit cell (OUC) along the surface normal,
modulo symmetry. If pymatgen's clustering (ftol) merges planes that are
actually distinct, terminations are silently lost.

This script:
1. Builds the OUC for a given Miller index
2. Extracts z-coordinates of all atoms, identifies planes at fine resolution
3. Shows the gap structure — where cuts happen and which gaps are close to ftol
4. Sweeps ftol to show how termination count varies
5. Flags any ftol-sensitive terminations

Usage:
    python termination_audit.py --input bulk.vasp --miller 1 1 0
    python termination_audit.py --input bulk.vasp --miller 1 1 0 --ftol-range 0.01 0.5 20
"""

import argparse
import json
import sys
from pathlib import Path


def load_structure(filepath):
    """Load a structure from file."""
    from pymatgen.core import Structure

    filepath = Path(filepath)
    if not filepath.exists():
        print(f"ERROR: File not found: {filepath}", file=sys.stderr)
        sys.exit(1)

    if filepath.suffix == ".json":
        with open(filepath) as f:
            return Structure.from_dict(json.load(f))
    return Structure.from_file(str(filepath))


def get_ouc_and_planes(structure, miller_index, fine_tol=0.01):
    """Get oriented unit cell and identify atomic planes along surface normal.

    Returns the OUC, z-coordinates, and plane assignments at fine tolerance.
    """
    import numpy as np
    from pymatgen.core.surface import SlabGenerator

    gen = SlabGenerator(
        structure,
        miller_index=miller_index,
        min_slab_size=1.0,  # doesn't matter, we only want the OUC
        min_vacuum_size=1.0,
        center_slab=False,
        lll_reduce=True,
        primitive=True,
    )

    ouc = gen.oriented_unit_cell

    # Get fractional z-coordinates and convert to Cartesian height
    frac_z = np.array([site.frac_coords[2] for site in ouc])
    cart_z = frac_z * ouc.lattice.c  # project onto c-axis

    # Sort atoms by z
    sort_idx = np.argsort(cart_z)
    sorted_z = cart_z[sort_idx]
    sorted_species = [ouc[i].species_string for i in sort_idx]

    # Cluster into planes using fine tolerance
    planes = []
    current_plane = [0]
    for i in range(1, len(sorted_z)):
        if sorted_z[i] - sorted_z[current_plane[-1]] < fine_tol:
            current_plane.append(i)
        else:
            planes.append(current_plane)
            current_plane = [i]
    planes.append(current_plane)

    # Compute plane z-positions (mean of atoms in each plane)
    plane_z = [np.mean(sorted_z[p]) for p in planes]
    plane_species = [
        ", ".join(sorted([sorted_species[j] for j in p]))
        for p in planes
    ]

    # Compute gaps between adjacent planes (including periodic wrap)
    gaps = []
    for i in range(len(plane_z) - 1):
        gaps.append(plane_z[i + 1] - plane_z[i])
    # Periodic gap (top to bottom of next repeat)
    periodic_gap = (ouc.lattice.c - plane_z[-1]) + plane_z[0]
    gaps.append(periodic_gap)

    return ouc, sorted_z, sorted_species, planes, plane_z, plane_species, gaps


def sweep_ftol(structure, miller_index, ftol_min=0.01, ftol_max=0.5, n_steps=20,
               min_slab_size=10.0, min_vacuum_size=15.0):
    """Sweep ftol and count terminations at each value."""
    import numpy as np
    from pymatgen.core.surface import SlabGenerator

    ftols = np.linspace(ftol_min, ftol_max, n_steps)
    results = []

    gen = SlabGenerator(
        structure,
        miller_index=miller_index,
        min_slab_size=min_slab_size,
        min_vacuum_size=min_vacuum_size,
        center_slab=True,
        lll_reduce=True,
        primitive=True,
    )

    for ftol in ftols:
        slabs = gen.get_slabs(ftol=ftol)
        shifts = sorted([s.shift for s in slabs])
        results.append({
            'ftol': ftol,
            'n_terminations': len(slabs),
            'shifts': shifts,
        })

    return results


def audit(structure, miller_index, fine_tol=0.01, ftol_min=0.01, ftol_max=0.5,
          n_ftol_steps=20, min_slab_size=10.0, min_vacuum_size=15.0):
    """Full termination audit for a given Miller index."""
    import numpy as np

    miller_str = f"({''.join(str(i) for i in miller_index)})"

    print(f"\n{'=' * 70}")
    print(f"TERMINATION AUDIT: {miller_str}")
    print(f"{'=' * 70}")

    # Step 1: Analyze the oriented unit cell
    ouc, sorted_z, sorted_species, planes, plane_z, plane_species, gaps = \
        get_ouc_and_planes(structure, miller_index, fine_tol)

    print(f"\nOriented Unit Cell:")
    print(f"  Atoms: {ouc.num_sites}")
    print(f"  c-height: {ouc.lattice.c:.4f} Å")
    print(f"  Atomic planes (tol={fine_tol} Å): {len(planes)}")

    print(f"\n{'Plane':<8} {'z (Å)':<12} {'Species':<30} {'Gap to next (Å)':<18}")
    print("-" * 70)
    for i, (z, sp, gap) in enumerate(zip(plane_z, plane_species, gaps)):
        flag = " ← SMALL" if gap < 0.2 else ""
        print(f"  {i:<6} {z:<12.4f} {sp:<30} {gap:<18.4f}{flag}")

    # Step 2: Identify which gaps are vulnerable to ftol merging
    print(f"\n--- Gap Analysis ---")
    min_gap = min(gaps)
    max_gap = max(gaps)
    print(f"  Smallest gap: {min_gap:.4f} Å")
    print(f"  Largest gap:  {max_gap:.4f} Å")
    print(f"  Default ftol: 0.1 Å")

    vulnerable = [(i, g) for i, g in enumerate(gaps) if g < 0.3]
    if vulnerable:
        print(f"\n  WARNING: {len(vulnerable)} gap(s) < 0.3 Å — terminations may be lost at default ftol!")
        for i, g in vulnerable:
            print(f"    Gap {i}→{(i+1) % len(planes)}: {g:.4f} Å "
                  f"(between {plane_species[i]} and {plane_species[(i+1) % len(planes)]})")
        print(f"\n  Recommendation: use --ftol {min_gap * 0.4:.3f} to resolve all planes")
    else:
        print(f"\n  All gaps > 0.3 Å — default ftol=0.1 should find all terminations.")

    # Step 3: Sweep ftol
    print(f"\n--- ftol Sweep ---")
    results = sweep_ftol(
        structure, miller_index,
        ftol_min=ftol_min, ftol_max=ftol_max, n_steps=n_ftol_steps,
        min_slab_size=min_slab_size, min_vacuum_size=min_vacuum_size,
    )

    print(f"{'ftol (Å)':<12} {'Terminations':<15} {'Shifts'}")
    print("-" * 60)
    prev_n = None
    for r in results:
        marker = ""
        if prev_n is not None and r['n_terminations'] != prev_n:
            marker = " ← CHANGE"
        shifts_str = ", ".join(f"{s:.4f}" for s in r['shifts'])
        print(f"  {r['ftol']:<10.4f} {r['n_terminations']:<15} [{shifts_str}]{marker}")
        prev_n = r['n_terminations']

    # Step 4: Summary
    max_terminations = max(r['n_terminations'] for r in results)
    default_result = min(results, key=lambda r: abs(r['ftol'] - 0.1))

    print(f"\n--- Summary ---")
    print(f"  Max terminations found (across all ftol): {max_terminations}")
    print(f"  Terminations at ftol=0.1: {default_result['n_terminations']}")

    if default_result['n_terminations'] < max_terminations:
        # Find the ftol that gives max terminations
        best_ftol = [r['ftol'] for r in results if r['n_terminations'] == max_terminations]
        print(f"  *** DEFAULT ftol=0.1 MISSES TERMINATIONS ***")
        print(f"  Use ftol ≤ {max(best_ftol):.4f} to find all {max_terminations}")
    else:
        print(f"  Default ftol=0.1 finds all terminations. ✓")

    print(f"{'=' * 70}\n")

    return {
        'miller_index': miller_index,
        'n_planes': len(planes),
        'gaps': gaps,
        'min_gap': min_gap,
        'max_terminations': max_terminations,
        'default_terminations': default_result['n_terminations'],
        'sweep': results,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Audit terminations — verify pymatgen finds all possible surface cuts"
    )
    parser.add_argument("--input", "-i", required=True,
                        help="Input bulk structure file")
    parser.add_argument("--miller", "-m", nargs=3, type=int, required=True,
                        help="Miller indices h k l")
    parser.add_argument("--fine-tol", type=float, default=0.01,
                        help="Fine tolerance for plane detection (default: 0.01 Å)")
    parser.add_argument("--ftol-range", nargs=3, type=float, default=[0.01, 0.5, 20],
                        help="ftol sweep: min max n_steps (default: 0.01 0.5 20)")
    parser.add_argument("--min-slab-size", type=float, default=10.0,
                        help="Min slab size for ftol sweep (default: 10.0)")
    parser.add_argument("--min-vacuum-size", type=float, default=15.0,
                        help="Min vacuum size for ftol sweep (default: 15.0)")

    args = parser.parse_args()

    structure = load_structure(args.input)
    miller_index = tuple(args.miller)

    audit(
        structure, miller_index,
        fine_tol=args.fine_tol,
        ftol_min=args.ftol_range[0],
        ftol_max=args.ftol_range[1],
        n_ftol_steps=int(args.ftol_range[2]),
        min_slab_size=args.min_slab_size,
        min_vacuum_size=args.min_vacuum_size,
    )


if __name__ == "__main__":
    main()
