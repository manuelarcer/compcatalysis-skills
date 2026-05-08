#!/usr/bin/env python3
"""Compute surface energies from optimized slab and bulk structures.

E_surf = (E_slab - (N_slab / N_bulk) * E_bulk) / (2 * A)

Requires:
    pip install ase pymatgen

Usage:
    python surface_energy.py --bulk-traj bulk/opt.traj --slab-dirs clean/t0 clean/t1 --output clean/surfenergy.log
    python surface_energy.py --bulk-traj bulk/opt.traj --slab-dirs clean/t0 clean/t1 --slab-traj opt.traj
"""

import argparse
import sys
from pathlib import Path


def read_energy_from_traj(traj_path):
    """Read the final energy from an ASE trajectory file."""
    from ase.io import read

    traj_path = Path(traj_path)
    if not traj_path.exists():
        print(f"ERROR: Trajectory not found: {traj_path}", file=sys.stderr)
        sys.exit(1)

    atoms = read(str(traj_path), index=-1)
    energy = atoms.get_potential_energy()
    return energy, atoms


def read_surface_area(slab_atoms):
    """Compute surface area from slab cell vectors (|a x b|)."""
    import numpy as np

    cell = slab_atoms.get_cell()
    cross = np.cross(cell[0], cell[1])
    return np.linalg.norm(cross)


def compute_slab_thickness(slab_atoms, layer_tol=0.5):
    """Layer-mean slab thickness in Å: <z_top_layer> - <z_bottom_layer>.

    Atoms are sorted by z and split into atomic layers wherever the gap
    between consecutive z values exceeds `layer_tol` Å. The thickness is
    the difference between the mean z of the topmost and bottommost
    layers. This is more robust than a raw `z_max - z_min` span, which is
    skewed by single protruding atoms (e.g., dangling H on hydroxylated
    surfaces, an apical O on a buckled metal-oxide layer).

    Falls back to the raw span if no clear interlayer gap is detected.
    """
    import numpy as np

    z = np.sort(slab_atoms.get_positions()[:, 2])
    if len(z) < 2:
        return 0.0
    gaps = np.diff(z)
    breaks = np.where(gaps > layer_tol)[0]
    if len(breaks) == 0:
        return float(z[-1] - z[0])
    bottom = z[: breaks[0] + 1]
    top = z[breaks[-1] + 1 :]
    return float(top.mean() - bottom.mean())


def compute_surface_energy(e_slab, n_slab, e_bulk, n_bulk, area):
    """Compute surface energy in eV/Å² and J/m².

    E_surf = (E_slab - (N_slab / N_bulk) * E_bulk) / (2 * A)
    """
    e_surf_ev_a2 = (e_slab - (n_slab / n_bulk) * e_bulk) / (2 * area)
    e_surf_j_m2 = e_surf_ev_a2 * 16.0217663  # eV/Å² to J/m²
    return e_surf_ev_a2, e_surf_j_m2


def find_traj_in_dir(slab_dir, traj_name="opt.traj"):
    """Find trajectory file in a slab directory."""
    path = Path(slab_dir) / traj_name
    if path.exists():
        return path
    # Try common alternatives
    for alt in ["opt.traj", "slab_opt.traj"]:
        alt_path = Path(slab_dir) / alt
        if alt_path.exists():
            return alt_path
    # Check for any .traj file
    trajs = list(Path(slab_dir).glob("*.traj"))
    if len(trajs) == 1:
        return trajs[0]
    elif len(trajs) > 1:
        print(f"WARNING: Multiple .traj files in {slab_dir}, using {trajs[0]}", file=sys.stderr)
        return trajs[0]
    return None


def find_vasp_in_dir(slab_dir):
    """Find the unrelaxed slab VASP file to get termination info."""
    parent = Path(slab_dir).parent
    term_name = Path(slab_dir).name  # e.g. "t0"
    # Look for original slab file in parent directory
    for f in parent.glob("*.vasp"):
        if term_name in f.stem:
            return f
    return None


def write_surfenergy_log(output_path, bulk_info, slab_results):
    """Write the surface energy comparison log."""
    output_path = Path(output_path)

    lines = []
    lines.append("=" * 70)
    lines.append("SURFACE ENERGY COMPARISON")
    lines.append("=" * 70)
    lines.append("")
    lines.append(f"Bulk reference: {bulk_info['path']}")
    lines.append(f"  Energy: {bulk_info['energy']:.6f} eV")
    lines.append(f"  Atoms:  {bulk_info['n_atoms']}")
    lines.append(f"  E/atom: {bulk_info['energy'] / bulk_info['n_atoms']:.6f} eV")
    lines.append("")

    # Sort by surface energy
    slab_results.sort(key=lambda x: x['e_surf_j_m2'])

    # Table header
    has_thickness = any('thickness' in r for r in slab_results)
    if has_thickness:
        lines.append(f"{'Term':<8} {'E_slab (eV)':<16} {'Atoms':<8} {'Area (Å²)':<12} "
                     f"{'Thickness (Å)':<14} {'E_surf (eV/Å²)':<16} "
                     f"{'E_surf (J/m²)':<14} {'Stable?':<8}")
        lines.append("-" * 100)
    else:
        lines.append(f"{'Term':<8} {'E_slab (eV)':<16} {'Atoms':<8} {'Area (Å²)':<12} "
                     f"{'E_surf (eV/Å²)':<16} {'E_surf (J/m²)':<14} {'Stable?':<8}")
        lines.append("-" * 86)

    for i, r in enumerate(slab_results):
        stable = "<<" if i == 0 else ""
        if has_thickness:
            t = r.get('thickness', float('nan'))
            lines.append(f"{r['name']:<8} {r['e_slab']:<16.6f} {r['n_atoms']:<8} "
                         f"{r['area']:<12.4f} {t:<14.4f} "
                         f"{r['e_surf_ev_a2']:<16.6f} "
                         f"{r['e_surf_j_m2']:<14.4f} {stable:<8}")
        else:
            lines.append(f"{r['name']:<8} {r['e_slab']:<16.6f} {r['n_atoms']:<8} "
                         f"{r['area']:<12.4f} {r['e_surf_ev_a2']:<16.6f} "
                         f"{r['e_surf_j_m2']:<14.4f} {stable:<8}")

    lines.append("")
    best = slab_results[0]
    lines.append(f"Most stable termination: {best['name']} "
                 f"(E_surf = {best['e_surf_j_m2']:.4f} J/m²)")

    if len(slab_results) > 1:
        diff = slab_results[1]['e_surf_j_m2'] - slab_results[0]['e_surf_j_m2']
        lines.append(f"Energy difference to next: {diff:.4f} J/m²")

    if has_thickness:
        lines.append("")
        lines.append("Thickness = <z_top_layer> - <z_bottom_layer> "
                     "(layer-mean, robust to dangling atoms; layer_tol = 0.5 Å).")

    lines.append("=" * 70)

    content = "\n".join(lines) + "\n"

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(content)

    # Also print to stdout
    print(content)

    return best


def main():
    parser = argparse.ArgumentParser(
        description="Compute and compare surface energies from optimized slabs"
    )
    parser.add_argument("--bulk-traj", required=True,
                        help="Path to bulk optimization trajectory")
    parser.add_argument("--slab-dirs", nargs="+", required=True,
                        help="Directories containing optimized slab trajectories (e.g., clean/t0 clean/t1)")
    parser.add_argument("--slab-traj", default="opt.traj",
                        help="Trajectory filename within each slab dir (default: opt.traj)")
    parser.add_argument("--output", "-o", default="surfenergy.log",
                        help="Output log file (default: surfenergy.log)")

    args = parser.parse_args()

    # Read bulk
    e_bulk, bulk_atoms = read_energy_from_traj(args.bulk_traj)
    n_bulk = len(bulk_atoms)
    bulk_info = {
        'path': args.bulk_traj,
        'energy': e_bulk,
        'n_atoms': n_bulk,
    }

    # Process each slab directory
    slab_results = []
    for slab_dir in args.slab_dirs:
        slab_dir = Path(slab_dir)
        traj_path = find_traj_in_dir(slab_dir, args.slab_traj)

        if traj_path is None:
            print(f"WARNING: No trajectory found in {slab_dir}, skipping", file=sys.stderr)
            continue

        e_slab, slab_atoms = read_energy_from_traj(traj_path)
        n_slab = len(slab_atoms)
        area = read_surface_area(slab_atoms)
        thickness = compute_slab_thickness(slab_atoms)

        e_surf_ev_a2, e_surf_j_m2 = compute_surface_energy(
            e_slab, n_slab, e_bulk, n_bulk, area
        )

        slab_results.append({
            'name': slab_dir.name,
            'path': str(slab_dir),
            'traj': str(traj_path),
            'e_slab': e_slab,
            'n_atoms': n_slab,
            'area': area,
            'thickness': thickness,
            'e_surf_ev_a2': e_surf_ev_a2,
            'e_surf_j_m2': e_surf_j_m2,
        })

    if not slab_results:
        print("ERROR: No valid slab trajectories found", file=sys.stderr)
        sys.exit(1)

    write_surfenergy_log(args.output, bulk_info, slab_results)


if __name__ == "__main__":
    main()
