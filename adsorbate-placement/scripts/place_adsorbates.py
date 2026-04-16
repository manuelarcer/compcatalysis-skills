#!/usr/bin/env python3
"""Place adsorbates on surface ontop sites of a slab.

Supports built-in OER adsorbates (O, OH, OOH) and arbitrary adsorbates
loaded from structure files. Detects surface sites by species, computes
bonding height from covalent radii, and writes one structure per
(adsorbate, site) combination.

Requires:
    pip install ase numpy

Usage:
    # Built-in adsorbates
    python place_adsorbates.py --slab slab.vasp --adsorbates O OH OOH --sites Co O

    # Arbitrary adsorbate from file (binding atom index 0 = first atom)
    python place_adsorbates.py --slab slab.vasp --adsorbate-file CCH.xyz --binding-atom 0 --sites Co

    # Override auto-computed height
    python place_adsorbates.py --slab slab.vasp --adsorbates OH --sites Co --height 2.0
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from ase.io import read, write
from ase import Atom, Atoms
from ase.constraints import FixAtoms
from ase.data import covalent_radii, atomic_numbers


# -- Height estimation --

def estimate_height(surface_symbol, binding_symbol):
    """Estimate bonding height from covalent radii of surface and binding atoms.

    Returns the sum of covalent radii, which is a reasonable initial guess
    for the bond length. The MLIP relaxation will refine this.

    Parameters
    ----------
    surface_symbol : str
        Chemical symbol of the surface atom (e.g. 'Co', 'Pt').
    binding_symbol : str
        Chemical symbol of the adsorbate atom that binds to the surface.

    Returns
    -------
    float
        Estimated height in Å.
    """
    r_surf = covalent_radii[atomic_numbers[surface_symbol]]
    r_bind = covalent_radii[atomic_numbers[binding_symbol]]
    return r_surf + r_bind


# -- Adsorbate geometry builders (built-in presets) --

def build_O():
    """Single oxygen atom at the origin."""
    return [("O", np.array([0.0, 0.0, 0.0]))], "O"


def build_OH(d_OH=0.98):
    """OH: O at origin, H straight up."""
    atoms = [
        ("O", np.array([0.0, 0.0, 0.0])),
        ("H", np.array([0.0, 0.0, d_OH])),
    ]
    return atoms, "O"


def build_OOH(tilt_deg=60, d_OO=1.21, d_OH=0.98):
    """OOH: O1 at origin, O2 tilted from vertical, H extends along O1->O2.

    The tilt angle is measured from the surface normal (z-axis).
    O2 is placed at d_OO from O1, tilted toward the surface plane.
    H is placed at d_OH from O2, extending along the O1->O2 direction.
    """
    tilt_rad = np.radians(tilt_deg)
    o2_offset = np.array([
        d_OO * np.sin(tilt_rad),
        0.0,
        d_OO * np.cos(tilt_rad),
    ])
    direction = o2_offset / np.linalg.norm(o2_offset)
    h_offset = o2_offset + d_OH * direction
    atoms = [
        ("O", np.array([0.0, 0.0, 0.0])),
        ("O", o2_offset),
        ("H", h_offset),
    ]
    return atoms, "O"


# Maps name -> (builder_func, returns (atoms_list, binding_symbol))
ADSORBATE_BUILDERS = {
    "O": build_O,
    "OH": build_OH,
    "OOH": build_OOH,
}


def load_adsorbate_from_file(filepath, binding_atom_index=0):
    """Load an arbitrary adsorbate from a structure file.

    The adsorbate is recentered so that the binding atom is at the origin.
    Atoms are oriented with the binding atom closest to z=0 (surface side).

    Parameters
    ----------
    filepath : str or Path
        Path to adsorbate structure file (XYZ, VASP, CIF, etc.).
    binding_atom_index : int
        Index of the atom that binds to the surface (0-based).

    Returns
    -------
    tuple of (list of (symbol, offset), str)
        Adsorbate atoms as (symbol, offset_from_binding_atom) pairs,
        and the chemical symbol of the binding atom.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        print(f"ERROR: Adsorbate file not found: {filepath}", file=sys.stderr)
        sys.exit(1)

    mol = read(str(filepath))
    if binding_atom_index >= len(mol):
        print(f"ERROR: Binding atom index {binding_atom_index} out of range "
              f"(adsorbate has {len(mol)} atoms, indices 0-{len(mol)-1})",
              file=sys.stderr)
        sys.exit(1)

    symbols = mol.get_chemical_symbols()
    positions = mol.get_positions()
    binding_symbol = symbols[binding_atom_index]
    binding_pos = positions[binding_atom_index].copy()

    # Recenter: binding atom at origin
    offsets = positions - binding_pos

    # Reorder: binding atom first, then sorted by distance from binding atom
    distances = np.linalg.norm(offsets, axis=1)
    order = [binding_atom_index] + sorted(
        [i for i in range(len(mol)) if i != binding_atom_index],
        key=lambda i: distances[i]
    )

    atoms_list = []
    for i in order:
        atoms_list.append((symbols[i], offsets[i].copy()))

    # Orient so that the adsorbate extends upward (positive z) from binding atom
    # Find center of mass of non-binding atoms relative to binding atom
    if len(mol) > 1:
        non_binding_offsets = np.array([offsets[i] for i in order[1:]])
        com_direction = non_binding_offsets.mean(axis=0)
        if com_direction[2] < 0:
            # Flip z: adsorbate was pointing downward
            for j in range(len(atoms_list)):
                sym, off = atoms_list[j]
                atoms_list[j] = (sym, off * np.array([1, 1, -1]))

    return atoms_list, binding_symbol


# -- Surface site detection --

def find_surface_sites(atoms, depth=2.0, site_species=None):
    """Find ontop surface sites grouped by species.

    Parameters
    ----------
    atoms : ase.Atoms
        The slab structure.
    depth : float
        Depth from topmost atom to consider as surface (Å).
    site_species : list of str, optional
        Which species to include. If None, all species found at surface.

    Returns
    -------
    dict[str, list[dict]]
        Maps species symbol to list of site dicts with keys:
        'index' (atom index), 'position' (xyz array), 'label' (e.g. 'Co_top_0').
    """
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    z_max = positions[:, 2].max()
    z_cutoff = z_max - depth

    # Get fractional coordinates for sorting
    try:
        frac = atoms.get_scaled_positions()
    except Exception:
        frac = positions

    # Collect surface atoms
    surface_atoms = []
    for i, (pos, sym) in enumerate(zip(positions, symbols)):
        if pos[2] > z_cutoff:
            if site_species is None or sym in site_species:
                surface_atoms.append({
                    "index": i,
                    "symbol": sym,
                    "position": pos.copy(),
                    "frac_xy": (frac[i, 0], frac[i, 1]),
                })

    # Group by species
    from collections import defaultdict
    groups = defaultdict(list)
    for sa in surface_atoms:
        groups[sa["symbol"]].append(sa)

    # Sort each group by fractional (x, y) for reproducible labeling
    result = {}
    for sym in sorted(groups.keys()):
        group = sorted(groups[sym], key=lambda a: (round(a["frac_xy"][0], 4),
                                                    round(a["frac_xy"][1], 4)))
        for i, site in enumerate(group):
            site["label"] = f"{sym}_top_{i}"
        result[sym] = group

    return result


# -- Adsorbate placement --

def place_adsorbate(slab, site_index, adsorbate_atoms, height):
    """Place an adsorbate on a slab at a given site.

    Parameters
    ----------
    slab : ase.Atoms
        Clean slab (will be copied, not modified).
    site_index : int
        Index of the surface atom to adsorb on.
    adsorbate_atoms : list of (symbol, offset)
        Adsorbate geometry: list of (symbol, xyz_offset) where offsets are
        relative to the binding atom position.
    height : float
        Height above the site atom (Å).

    Returns
    -------
    ase.Atoms
        New slab with adsorbate appended.
    """
    new_slab = slab.copy()

    # Preserve existing FixAtoms constraint
    old_fixed = set()
    for c in slab.constraints:
        if isinstance(c, FixAtoms):
            old_fixed.update(c.index)

    site_pos = slab.positions[site_index]
    base = np.array([site_pos[0], site_pos[1], site_pos[2] + height])

    for symbol, offset in adsorbate_atoms:
        new_slab.append(Atom(symbol, position=base + offset))

    # Re-apply constraints (only to original fixed atoms)
    if old_fixed:
        new_slab.set_constraint(FixAtoms(indices=list(old_fixed)))

    return new_slab


# -- Main generation logic --

def generate_all_placements(slab, adsorbate_specs, site_species, output_dir,
                            height=None, depth=2.0):
    """Generate all (adsorbate, site) combinations and write structures.

    Parameters
    ----------
    slab : ase.Atoms
        Clean slab structure.
    adsorbate_specs : list of dict
        Each dict has keys: 'name' (str), 'atoms' (list of (symbol, offset)),
        'binding_symbol' (str).
    site_species : list of str
        Which surface species to place on.
    output_dir : Path
        Root output directory.
    height : float or None
        Height above site atom (Å). If None, auto-computed from covalent radii.
    depth : float
        Depth for surface atom detection (Å).

    Returns
    -------
    list of dict
        Summary of each placement.
    """
    output_dir = Path(output_dir)
    sites = find_surface_sites(slab, depth=depth, site_species=site_species)

    if not sites:
        print("ERROR: No surface sites found for the specified species.", file=sys.stderr)
        print(f"  Species requested: {site_species}", file=sys.stderr)
        surface_syms = set(slab.get_chemical_symbols()[i]
                          for i, z in enumerate(slab.positions[:, 2])
                          if z > slab.positions[:, 2].max() - depth)
        print(f"  Species at surface (depth={depth} Å): {sorted(surface_syms)}", file=sys.stderr)
        sys.exit(1)

    results = []

    for spec in adsorbate_specs:
        ads_name = spec["name"]
        ads_atoms = spec["atoms"]
        binding_sym = spec["binding_symbol"]
        ads_dir = output_dir / ads_name

        for sym, site_list in sites.items():
            # Compute height for this (surface_species, binding_atom) pair
            if height is not None:
                h = height
            else:
                h = estimate_height(sym, binding_sym)

            for site in site_list:
                label = site["label"]
                site_dir = ads_dir / label
                site_dir.mkdir(parents=True, exist_ok=True)

                placed = place_adsorbate(slab, site["index"], ads_atoms, h)
                out_path = site_dir / "input.vasp"
                write(str(out_path), placed, format="vasp")

                results.append({
                    "adsorbate": ads_name,
                    "site_label": label,
                    "site_index": site["index"],
                    "site_species": sym,
                    "binding_symbol": binding_sym,
                    "height": h,
                    "n_atoms": len(placed),
                    "path": str(out_path),
                })

    return results


def print_summary(sites, results):
    """Print a summary of detected sites and generated structures."""
    print("\n=== Surface Sites ===")
    total_sites = sum(len(v) for v in sites.values())
    print(f"Found {total_sites} ontop site(s):")
    for sym in sorted(sites.keys()):
        print(f"  {sym}: {len(sites[sym])} site(s)")
        for s in sites[sym]:
            pos = s["position"]
            print(f"    {s['label']:>15s}  atom {s['index']:3d}  "
                  f"({pos[0]:7.3f}, {pos[1]:7.3f}, {pos[2]:7.3f})")

    print(f"\n=== Generated Structures ===")
    print(f"{'Adsorbate':<12} {'Site':<15} {'Height':<8} {'Atoms':<8} {'Path'}")
    print("-" * 80)
    for r in results:
        print(f"{r['adsorbate']:<12} {r['site_label']:<15} {r['height']:<8.2f} "
              f"{r['n_atoms']:<8} {r['path']}")
    print(f"\nTotal: {len(results)} structure(s) written.")


def main():
    parser = argparse.ArgumentParser(
        description="Place adsorbates on surface ontop sites.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  # Built-in OER adsorbates
  %(prog)s --slab slab.vasp --adsorbates O OH OOH --sites Co O

  # Arbitrary adsorbate from file (binds through atom 0)
  %(prog)s --slab slab.vasp --adsorbate-file CCH.xyz --binding-atom 0 --sites Co

  # Override auto-height
  %(prog)s --slab slab.vasp --adsorbates OH --sites Co --height 2.0

  # Multiple file-based adsorbates
  %(prog)s --slab slab.vasp --adsorbate-file CO.xyz HCOO.xyz --binding-atom 0 0 --sites Pt
""",
    )
    parser.add_argument("--slab", required=True,
                        help="Input clean slab structure file (VASP, CIF, XYZ)")
    parser.add_argument("--adsorbates", "-a", nargs="+",
                        help="Built-in adsorbate(s): O, OH, OOH")
    parser.add_argument("--adsorbate-file", "-f", nargs="+",
                        help="Adsorbate structure file(s) (XYZ, VASP, CIF)")
    parser.add_argument("--binding-atom", "-b", nargs="+", type=int, default=[0],
                        help="Index of binding atom in each adsorbate file (default: 0)")
    parser.add_argument("--sites", "-s", nargs="+", required=True,
                        help="Surface species to adsorb onto (e.g. Co O Pt)")
    parser.add_argument("--output-dir", "-o", default=".",
                        help="Root output directory (default: current dir)")
    parser.add_argument("--height", type=float, default=None,
                        help="Height above surface site (Å). If omitted, "
                             "auto-computed from covalent radii of surface + binding atoms")
    parser.add_argument("--depth", type=float, default=2.0,
                        help="Depth from top for surface detection (Å, default: 2.0)")
    args = parser.parse_args()

    if not args.adsorbates and not args.adsorbate_file:
        parser.error("Provide --adsorbates and/or --adsorbate-file")

    # Validate input file
    slab_path = Path(args.slab)
    if not slab_path.exists():
        print(f"ERROR: Slab file not found: {slab_path}", file=sys.stderr)
        sys.exit(1)

    # Load slab
    slab = read(str(slab_path))
    print(f"Loaded slab: {slab.get_chemical_formula()} ({len(slab)} atoms)")

    # Build adsorbate specs
    adsorbate_specs = []

    # Built-in adsorbates
    if args.adsorbates:
        for name in args.adsorbates:
            if name not in ADSORBATE_BUILDERS:
                print(f"ERROR: Unknown built-in adsorbate '{name}'. "
                      f"Available: {', '.join(ADSORBATE_BUILDERS.keys())}",
                      file=sys.stderr)
                sys.exit(1)
            atoms, binding_sym = ADSORBATE_BUILDERS[name]()
            adsorbate_specs.append({
                "name": name,
                "atoms": atoms,
                "binding_symbol": binding_sym,
            })

    # File-based adsorbates
    if args.adsorbate_file:
        binding_indices = args.binding_atom
        # Extend binding indices if fewer than files
        if len(binding_indices) < len(args.adsorbate_file):
            binding_indices = binding_indices + [0] * (len(args.adsorbate_file) - len(binding_indices))

        for fpath, bidx in zip(args.adsorbate_file, binding_indices):
            atoms, binding_sym = load_adsorbate_from_file(fpath, bidx)
            name = Path(fpath).stem
            adsorbate_specs.append({
                "name": name,
                "atoms": atoms,
                "binding_symbol": binding_sym,
            })

    # Print parameters
    height_str = f"{args.height} Å" if args.height is not None else "auto (covalent radii)"
    print(f"\nParameters:")
    print(f"  Adsorbates:   {', '.join(s['name'] for s in adsorbate_specs)}")
    print(f"  Site species: {', '.join(args.sites)}")
    print(f"  Height:       {height_str}")
    print(f"  Depth:        {args.depth} Å")
    print(f"  Output dir:   {args.output_dir}")

    if args.height is None:
        print(f"\n  Auto-height estimates (covalent radii):")
        seen = set()
        for spec in adsorbate_specs:
            for sym in args.sites:
                pair = (sym, spec["binding_symbol"])
                if pair not in seen:
                    seen.add(pair)
                    h = estimate_height(sym, spec["binding_symbol"])
                    print(f"    {sym}—{spec['binding_symbol']} "
                          f"({spec['name']}): {h:.2f} Å")

    # Detect sites (for summary)
    sites = find_surface_sites(slab, depth=args.depth, site_species=args.sites)

    # Generate placements
    results = generate_all_placements(
        slab,
        adsorbate_specs=adsorbate_specs,
        site_species=args.sites,
        output_dir=args.output_dir,
        height=args.height,
        depth=args.depth,
    )

    print_summary(sites, results)


if __name__ == "__main__":
    main()
