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


# -- Orientation search --

def rotate_adsorbate_offsets(adsorbate_atoms, angle_deg):
    """Rotate adsorbate atom offsets around the surface normal (z-axis).

    Parameters
    ----------
    adsorbate_atoms : list of (symbol, offset)
        Offsets are relative to the binding atom at the origin.
    angle_deg : float
        Rotation angle in degrees (positive = counter-clockwise viewed from +z).

    Returns
    -------
    list of (symbol, offset)
        New list with rotated offsets. Symbols and the binding atom (which
        sits at the origin and is therefore invariant under z-rotation)
        are unchanged.
    """
    angle_rad = np.deg2rad(angle_deg)
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    R = np.array([[c, -s, 0.0],
                  [s,  c, 0.0],
                  [0.0, 0.0, 1.0]])
    return [(sym, R @ np.asarray(offset, dtype=float))
            for sym, offset in adsorbate_atoms]


def min_nonsite_distance(slab, site_index, placed_atoms):
    """Minimum distance between any appended adsorbate atom and any non-site slab atom.

    Excludes the site atom (intentional bonding partner). Used to score
    candidate orientations — larger = better clearance.
    """
    n_slab = len(slab)
    n_total = len(placed_atoms)
    if n_total <= n_slab:
        return float("inf")
    min_d = float("inf")
    for ai in range(n_slab, n_total):
        for si in range(n_slab):
            if si == site_index:
                continue
            d = placed_atoms.get_distance(ai, si, mic=True)
            if d < min_d:
                min_d = d
    return min_d


def find_best_orientation(slab, site_index, adsorbate_atoms, height,
                          n_rotations=12, clash_factor=0.7):
    """Try N rotations of the adsorbate around the surface normal and return
    the orientation that maximizes the minimum distance to non-site slab
    atoms. For adsorbates that are symmetric under z-rotation (single atom,
    or any geometry with all offsets on the z-axis) every angle gives the
    same result and the search collapses to a single placement.

    The strict ``clash_factor`` threshold is still applied at the end on
    the best-rotation placement; if even that orientation clashes, the
    caller may choose to skip.

    Parameters
    ----------
    slab : ase.Atoms
        Clean slab.
    site_index : int
        Site atom index in slab.
    adsorbate_atoms : list of (symbol, offset)
        Adsorbate geometry (binding atom at origin).
    height : float
        Bonding height above the site atom (Å).
    n_rotations : int
        Number of equally-spaced rotations to try in [0°, 360°).
        Set to 1 to disable orientation search (use the input geometry as-is).
    clash_factor : float
        Hard clash threshold for the final placement.

    Returns
    -------
    dict
        Keys: ``angle_deg`` (best angle), ``atoms`` (rotated offset list),
        ``placed`` (full slab+adsorbate with best rotation), ``min_distance``
        (min distance to any non-site slab atom), ``clashes`` (list of
        hard clashes at clash_factor for the best orientation).
    """
    if n_rotations < 1:
        n_rotations = 1
    angles = [i * 360.0 / n_rotations for i in range(n_rotations)]

    best = None
    for angle in angles:
        rotated = rotate_adsorbate_offsets(adsorbate_atoms, angle)
        placed = place_adsorbate(slab, site_index, rotated, height)
        min_d = min_nonsite_distance(slab, site_index, placed)
        if best is None or min_d > best["min_distance"]:
            best = {
                "angle_deg": float(angle),
                "atoms": rotated,
                "placed": placed,
                "min_distance": min_d,
            }
    best["clashes"] = detect_clashes(slab, site_index, best["placed"],
                                      clash_factor=clash_factor)
    return best


# -- Clash detection --

def detect_clashes(slab, site_index, placed_atoms, clash_factor=0.7):
    """Detect close contacts between newly placed adsorbate atoms and existing slab atoms.

    For each appended adsorbate atom (positions beyond ``len(slab)``), scan all
    slab atoms except the site atom (which is the intentional bonding partner)
    and flag any that sit closer than ``clash_factor × (r_cov(slab) + r_cov(ads))``.

    The default factor of 0.7 catches gross overlaps (e.g. a placed O landing
    < 0.7 Å from a lattice H — i.e. *inside* what would be the new O–H bond)
    without flagging ordinary close-packing in a corrugated surface.

    Parameters
    ----------
    slab : ase.Atoms
        Clean slab (used for length / site lookup).
    site_index : int
        Index of the site atom in the slab — excluded from the clash check
        because it sits at the intended bond distance from the binding atom.
    placed_atoms : ase.Atoms
        Full slab+adsorbate structure with adsorbate atoms appended at the end.
    clash_factor : float
        Multiplier on the sum of covalent radii. Lower = stricter (more
        clashes flagged); higher = more permissive.

    Returns
    -------
    list of dict
        One entry per clash with keys ``ads_idx, ads_symbol, slab_idx,
        slab_symbol, distance, threshold``.
    """
    n_slab = len(slab)
    n_total = len(placed_atoms)
    clashes = []
    for ai in range(n_slab, n_total):
        ads_sym = placed_atoms[ai].symbol
        r_ads = covalent_radii[atomic_numbers[ads_sym]]
        for si in range(n_slab):
            if si == site_index:
                continue
            slab_sym = placed_atoms[si].symbol
            r_slab = covalent_radii[atomic_numbers[slab_sym]]
            threshold = clash_factor * (r_slab + r_ads)
            d = placed_atoms.get_distance(ai, si, mic=True)
            if d < threshold:
                clashes.append({
                    "ads_idx": ai,
                    "ads_symbol": ads_sym,
                    "slab_idx": si,
                    "slab_symbol": slab_sym,
                    "distance": float(d),
                    "threshold": float(threshold),
                })
    return clashes


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
                            height=None, depth=2.0,
                            clash_factor=0.7, allow_clashes=False,
                            n_rotations=12):
    """Generate all (adsorbate, site) combinations and write structures.

    For each (adsorbate, site) pair the script tries ``n_rotations`` equally
    spaced rotations of the adsorbate around the surface normal and keeps
    the orientation that maximizes the minimum distance from any adsorbate
    atom to any non-site slab atom. This avoids tilted adsorbates (e.g. the
    standard 60°-tilted OOH) accidentally pointing into a neighboring
    lattice atom — without changing the adsorbate's intrinsic geometry.

    By default, placements that still produce a hard clash (an adsorbate
    atom within ``clash_factor × (r_cov_slab + r_cov_ads)`` of any non-site
    slab atom) after the orientation search are skipped — these almost
    always come from picking a site that is already saturated (e.g. a
    lattice O that already has an H above it, where adding a new O at
    sum-of-radii height lands the new O ~0.4 Å from the H). Set
    ``allow_clashes=True`` to write them anyway (e.g. if the user is
    intentionally studying replacement / proton-transfer chemistry).

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
    clash_factor : float
        Multiplier on covalent-radii sum used as the clash threshold.
    allow_clashes : bool
        If True, write clashing placements with a warning. Default False (skip).

    Returns
    -------
    tuple (list, list)
        ``(written, skipped)`` — list of placements actually written and list
        of placements skipped due to clashes. Each entry is a dict with the
        same keys; skipped entries also carry ``"clashes"``.
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

    written = []
    skipped = []

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
                best = find_best_orientation(
                    slab, site["index"], ads_atoms, h,
                    n_rotations=n_rotations, clash_factor=clash_factor,
                )
                placed = best["placed"]
                clashes = best["clashes"]

                entry = {
                    "adsorbate": ads_name,
                    "site_label": label,
                    "site_index": site["index"],
                    "site_species": sym,
                    "binding_symbol": binding_sym,
                    "height": h,
                    "rotation_deg": best["angle_deg"],
                    "min_nonsite_distance": best["min_distance"],
                    "n_atoms": len(placed),
                    "clashes": clashes,
                }

                if clashes and not allow_clashes:
                    entry["path"] = None
                    skipped.append(entry)
                    continue

                site_dir = ads_dir / label
                site_dir.mkdir(parents=True, exist_ok=True)
                out_path = site_dir / "input.vasp"
                write(str(out_path), placed, format="vasp")
                entry["path"] = str(out_path)
                written.append(entry)

    return written, skipped


def print_summary(sites, written, skipped):
    """Print a summary of detected sites, generated structures, and skipped clashes."""
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
    print(f"{'Adsorbate':<12} {'Site':<14} {'Height':<7} {'Rot°':<6} "
          f"{'min-d (Å)':<10} {'Path'}")
    print("-" * 92)
    for r in written:
        warn = "  ⚠ clash kept" if r["clashes"] else ""
        rot = r.get("rotation_deg", 0.0)
        mind = r.get("min_nonsite_distance", float("nan"))
        print(f"{r['adsorbate']:<12} {r['site_label']:<14} {r['height']:<7.2f} "
              f"{rot:<6.0f} {mind:<10.3f} {r['path']}{warn}")
    print(f"\nTotal: {len(written)} structure(s) written.")

    if skipped:
        print(f"\n=== Skipped (clash detected — pass --allow-clashes to keep) ===")
        print(f"{'Adsorbate':<12} {'Site':<15} {'Closest non-site contact':<55}")
        print("-" * 80)
        for r in skipped:
            c = min(r["clashes"], key=lambda x: x["distance"])
            contact = (f"{c['ads_symbol']}{c['ads_idx']}↔{c['slab_symbol']}{c['slab_idx']}"
                       f"  d={c['distance']:.3f} Å  (threshold {c['threshold']:.3f})")
            print(f"{r['adsorbate']:<12} {r['site_label']:<15} {contact}")
        print(f"\nSkipped: {len(skipped)} placement(s) due to clashes.")


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
    parser.add_argument("--clash-factor", type=float, default=0.7,
                        help="Clash threshold as a fraction of "
                             "(r_cov_slab + r_cov_ads). Lower = stricter. "
                             "Default 0.7 catches gross overlaps "
                             "(e.g. new O on a lattice OH where the existing H "
                             "would sit ~0.4 Å from the new O).")
    parser.add_argument("--allow-clashes", action="store_true",
                        help="Write clashing placements anyway (with a "
                             "warning). Use when intentionally studying "
                             "replacement / proton-transfer chemistry.")
    parser.add_argument("--rotations", type=int, default=12,
                        help="Number of rotations of the adsorbate around "
                             "the surface normal to try, picking the one "
                             "that maximizes clearance from non-site slab "
                             "atoms (default: 12). Set to 1 to disable "
                             "orientation search.")
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
    written, skipped = generate_all_placements(
        slab,
        adsorbate_specs=adsorbate_specs,
        site_species=args.sites,
        output_dir=args.output_dir,
        height=args.height,
        depth=args.depth,
        clash_factor=args.clash_factor,
        allow_clashes=args.allow_clashes,
        n_rotations=args.rotations,
    )

    print_summary(sites, written, skipped)


if __name__ == "__main__":
    main()
