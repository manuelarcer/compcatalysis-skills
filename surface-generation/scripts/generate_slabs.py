#!/usr/bin/env python3
"""Generate surface slabs from bulk crystal structures using pymatgen.

Requires:
    pip install pymatgen ase

Usage:
    python generate_slabs.py --input Pt_bulk.vasp --miller 1 1 1
    python generate_slabs.py --input TiO2_bulk.vasp --miller 1 1 0 --all-terminations
    python generate_slabs.py --input Pt_bulk.vasp --max-index 2 --output-dir Pt_surfaces/
    python generate_slabs.py --input Pt_bulk.vasp --max-index 2 --list-surfaces
"""

import argparse
import json
import sys
from pathlib import Path


def load_structure(filepath):
    """Load a structure from file (POSCAR, CIF, XYZ, JSON)."""
    from pymatgen.core import Structure

    filepath = Path(filepath)
    if not filepath.exists():
        print(f"ERROR: Input file not found: {filepath}", file=sys.stderr)
        print("Provide a bulk structure file (POSCAR, CIF, XYZ, or JSON).", file=sys.stderr)
        sys.exit(1)

    if filepath.suffix == ".json":
        with open(filepath) as f:
            return Structure.from_dict(json.load(f))
    return Structure.from_file(str(filepath))


def parse_bonds(bond_strings):
    """Parse bond strings like 'Ti-O:2.0' into dict for SlabGenerator.

    Returns dict mapping (species1, species2) -> max_bond_length.
    """
    if not bond_strings:
        return None

    from pymatgen.core import Element

    bonds = {}
    for bs in bond_strings:
        try:
            pair, dist = bs.split(":")
            sp1, sp2 = pair.split("-")
            bonds[(Element(sp1), Element(sp2))] = float(dist)
        except (ValueError, KeyError) as e:
            print(f"ERROR: Invalid bond format '{bs}'. Use 'Element1-Element2:distance'",
                  file=sys.stderr)
            print(f"  Example: Ti-O:2.0", file=sys.stderr)
            sys.exit(1)
    return bonds


def get_top_surface_signature(slab, depth=1.5, z_bin=0.1):
    """Get a structural signature for the top surface of a slab.

    Considers both composition AND z-structure of the top atoms, so that
    surfaces with the same species counts but different corrugation or
    stacking are distinguished.

    Atoms within `depth` Å of the topmost atom are grouped into sub-layers
    (binned by z with tolerance `z_bin`), and the signature encodes
    species counts per sub-layer from top to bottom.
    """
    import numpy as np

    z = np.array([site.coords[2] for site in slab])
    z_max = z.max()

    # Collect top atoms
    top_atoms = []
    for site in slab:
        if site.coords[2] > z_max - depth:
            top_atoms.append((site.species_string, site.coords[2]))

    if not top_atoms:
        return "empty"

    # Sort by z descending (topmost first)
    top_atoms.sort(key=lambda x: -x[1])

    # Bin into sub-layers
    layers = []
    current_layer = [top_atoms[0]]
    for sp, zz in top_atoms[1:]:
        if current_layer[-1][1] - zz < z_bin:
            current_layer.append((sp, zz))
        else:
            layers.append(current_layer)
            current_layer = [(sp, zz)]
    layers.append(current_layer)

    # Build signature: each sub-layer as sorted species counts
    parts = []
    for layer in layers:
        species_count = {}
        for sp, _ in layer:
            species_count[sp] = species_count.get(sp, 0) + 1
        parts.append("_".join(f"{k}{v}" for k, v in sorted(species_count.items())))

    return "|".join(parts)


def filter_unique_top_surfaces(slabs, depth=1.5):
    """Filter slabs to keep only those with unique top-surface compositions.

    pymatgen's StructureMatcher considers a slab equivalent to its
    top/bottom flip, but for catalysis only the exposed (top) surface
    matters. This function deduplicates by top-surface composition.
    """
    seen = {}
    unique = []
    for slab in slabs:
        sig = get_top_surface_signature(slab, depth)
        if sig not in seen:
            seen[sig] = slab.shift
            unique.append(slab)
    return unique


def is_stoichiometric(slab, bulk_structure):
    """Check whether a slab has the same reduced composition as its bulk parent."""
    from pymatgen.core import Composition
    return (Composition(slab.composition).reduced_formula
            == Composition(bulk_structure.composition).reduced_formula)


def filter_sym_stoich(slabs, bulk_structure, require_symmetric=False,
                      require_stoichiometric=False):
    """Filter slabs by symmetry and/or stoichiometry constraints.

    Returns
    -------
    kept : list of Slab
        Slabs satisfying every active constraint.
    dropped : list of (slab, reason)
        Rejected slabs with a short reason string for reporting.
    """
    if not (require_symmetric or require_stoichiometric):
        return list(slabs), []

    kept, dropped = [], []
    for slab in slabs:
        reasons = []
        if require_symmetric and not slab.is_symmetric():
            reasons.append("asymmetric")
        if require_stoichiometric and not is_stoichiometric(slab, bulk_structure):
            reasons.append("non-stoichiometric")
        if reasons:
            dropped.append((slab, "+".join(reasons)))
        else:
            kept.append(slab)
    return kept, dropped


def generate_slabs_for_miller(structure, miller_index, min_slab_size=10.0,
                               min_vacuum_size=15.0, center_slab=True,
                               symmetrize=False, ftol=0.1,
                               max_broken_bonds=0, bonds=None):
    """Generate all unique slabs for a given Miller index.

    Disables pymatgen's symmetry filtering and instead deduplicates by
    top-surface composition — because for catalysis, a slab with surface A
    on top is NOT equivalent to the same slab flipped with surface B on top.
    """
    from pymatgen.core.surface import SlabGenerator

    gen = SlabGenerator(
        structure,
        miller_index=miller_index,
        min_slab_size=min_slab_size,
        min_vacuum_size=min_vacuum_size,
        center_slab=center_slab,
        lll_reduce=True,
        primitive=True,
    )

    # Get ALL slabs without pymatgen's symmetry filter
    slabs = gen.get_slabs(
        bonds=bonds,
        ftol=ftol,
        max_broken_bonds=max_broken_bonds,
        symmetrize=symmetrize,
        filter_out_sym_slabs=False,
    )

    # Deduplicate by top-surface composition (catalysis-aware)
    slabs = filter_unique_top_surfaces(slabs)

    return slabs


def get_distinct_miller_indices(structure, max_index):
    """Get all symmetrically distinct Miller indices up to max_index."""
    from pymatgen.core.surface import get_symmetrically_distinct_miller_indices

    return get_symmetrically_distinct_miller_indices(structure, max_index)


def apply_selective_dynamics(slab, fix_bottom=0.5):
    """Apply selective dynamics to a slab — fix atoms in the bottom fraction.

    Args:
        slab: pymatgen Slab or Structure object
        fix_bottom: fraction of atoms to fix (by z-coordinate), default 0.5.
                    Set to 0 to disable constraints.

    Returns:
        Structure with selective_dynamics site property set.
    """
    if fix_bottom <= 0:
        return slab

    import numpy as np

    z_coords = np.array([site.coords[2] for site in slab])
    z_min, z_max = z_coords.min(), z_coords.max()
    z_threshold = z_min + fix_bottom * (z_max - z_min)

    selective_dynamics = []
    n_fixed = 0
    for site in slab:
        if site.coords[2] <= z_threshold:
            selective_dynamics.append([False, False, False])
            n_fixed += 1
        else:
            selective_dynamics.append([True, True, True])

    slab.add_site_property("selective_dynamics", selective_dynamics)
    return slab, n_fixed


def write_slab(slab, output_path, fmt="vasp"):
    """Write a slab structure to file."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if fmt == "vasp":
        slab.to(filename=str(output_path), fmt="poscar")
    elif fmt == "cif":
        slab.to(filename=str(output_path), fmt="cif")
    elif fmt == "xyz":
        slab.to(filename=str(output_path), fmt="xyz")
    elif fmt == "json":
        with open(output_path, "w") as f:
            json.dump(slab.as_dict(), f, indent=2)
    else:
        raise ValueError(f"Unknown format: {fmt}. Use: vasp, cif, xyz, json")


def miller_str(miller_index):
    """Format Miller index as string, e.g. (1,1,1) -> '111'."""
    return "".join(str(i) for i in miller_index)


def miller_display(miller_index):
    """Format Miller index for display, e.g. (1,1,1) -> '(111)'."""
    return f"({''.join(str(i) for i in miller_index)})"


def slab_filename(formula, miller_index, term_index, fmt="vasp"):
    """Generate standard filename for a slab."""
    ext = {"vasp": ".vasp", "cif": ".cif", "xyz": ".xyz", "json": ".json"}
    return f"{formula}_{miller_str(miller_index)}_t{term_index}{ext.get(fmt, '.vasp')}"


def print_slab_summary(slab, miller_index, term_index, n_fixed=0):
    """Print summary of a generated slab."""
    sym = "yes" if slab.is_symmetric() else "no"
    polar = "yes" if slab.is_polar() else "no"
    lattice = slab.lattice
    fixed_str = f", fixed={n_fixed}/{slab.num_sites}" if n_fixed > 0 else ""
    print(f"  Termination {term_index}: shift={slab.shift:.4f}, "
          f"atoms={slab.num_sites}{fixed_str}, "
          f"area={slab.surface_area:.2f} Å², "
          f"symmetric={sym}, polar={polar}")
    print(f"    Lattice: a={lattice.a:.3f} b={lattice.b:.3f} c={lattice.c:.3f} Å")


def list_surfaces(structure, max_index, min_slab_size=10.0, min_vacuum_size=15.0,
                  ftol=0.1, bonds=None, max_broken_bonds=0):
    """List all symmetrically distinct surfaces up to max_index."""
    indices = get_distinct_miller_indices(structure, max_index)

    print(f"\nSymmetrically distinct surfaces (max index={max_index}):")
    print(f"{'Miller':<10} {'Terms':<8} {'Atoms (1st)':<12} {'Area (Å²)':<12} {'Symmetric':<12}")
    print("-" * 55)

    total_slabs = 0
    for mi in sorted(indices, key=lambda x: sum(abs(i) for i in x)):
        slabs = generate_slabs_for_miller(
            structure, mi,
            min_slab_size=min_slab_size,
            min_vacuum_size=min_vacuum_size,
            ftol=ftol,
            bonds=bonds,
            max_broken_bonds=max_broken_bonds,
        )
        if not slabs:
            print(f"{miller_display(mi):<10} {'0':<8} {'-':<12} {'-':<12} {'-':<12}")
            continue

        sym = "yes" if slabs[0].is_symmetric() else "no"
        print(f"{miller_display(mi):<10} {len(slabs):<8} {slabs[0].num_sites:<12} "
              f"{slabs[0].surface_area:<12.2f} {sym:<12}")
        total_slabs += len(slabs)

    print(f"\nTotal: {len(indices)} distinct Miller indices, {total_slabs} unique terminations")
    return indices


def process_single_miller(structure, miller_index, all_terminations=False,
                           output=None, output_dir=".", fmt="vasp",
                           min_slab_size=10.0, min_vacuum_size=15.0,
                           center_slab=True, symmetrize=False,
                           ftol=0.1, max_broken_bonds=0, bonds=None,
                           fix_bottom=0.5,
                           require_symmetric=False,
                           require_stoichiometric=False):
    """Generate slabs for a single Miller index."""
    from pymatgen.core import Composition

    formula = Composition(structure.composition).reduced_formula
    slabs = generate_slabs_for_miller(
        structure, miller_index,
        min_slab_size=min_slab_size,
        min_vacuum_size=min_vacuum_size,
        center_slab=center_slab,
        symmetrize=symmetrize,
        ftol=ftol,
        max_broken_bonds=max_broken_bonds,
        bonds=bonds,
    )

    if not slabs:
        print(f"No valid slabs found for {miller_display(miller_index)}.", file=sys.stderr)
        print("Try: --max-broken-bonds 1 or adjust --bonds", file=sys.stderr)
        sys.exit(1)

    n_total = len(slabs)
    slabs, dropped = filter_sym_stoich(
        slabs, structure,
        require_symmetric=require_symmetric,
        require_stoichiometric=require_stoichiometric,
    )
    if dropped:
        print(f"  Dropped {len(dropped)}/{n_total} termination(s) by filter:")
        for slab, reason in dropped:
            comp = slab.composition.reduced_formula
            print(f"    shift={slab.shift:.4f}, atoms={slab.num_sites}, "
                  f"comp={comp}, reason={reason}")
    if not slabs:
        print(f"All terminations rejected by --require-symmetric / "
              f"--require-stoichiometric for {miller_display(miller_index)}.",
              file=sys.stderr)
        sys.exit(1)

    print(f"\n{formula} {miller_display(miller_index)}: {len(slabs)} termination(s)")

    if all_terminations:
        output_dir = Path(output_dir)
        written = []
        for i, slab in enumerate(slabs):
            n_fixed = 0
            if fix_bottom > 0:
                slab, n_fixed = apply_selective_dynamics(slab, fix_bottom)
            print_slab_summary(slab, miller_index, i, n_fixed)
            term_dir = output_dir / f"t{i}"
            term_dir.mkdir(parents=True, exist_ok=True)
            fname = slab_filename(formula, miller_index, i, fmt)
            out_path = term_dir / fname
            write_slab(slab, out_path, fmt)
            print(f"    Saved: {out_path}")
            written.append((slab, out_path))
        return written
    else:
        slab = slabs[0]
        n_fixed = 0
        if fix_bottom > 0:
            slab, n_fixed = apply_selective_dynamics(slab, fix_bottom)
        print_slab_summary(slab, miller_index, 0, n_fixed)
        if len(slabs) > 1:
            print(f"  Note: {len(slabs) - 1} other termination(s) available. "
                  f"Use --all-terminations to generate all.")

        if output is None:
            output = slab_filename(formula, miller_index, 0, fmt)
        write_slab(slab, output, fmt)
        print(f"  Saved: {output}")
        return [(slab, output)]


def process_max_index(structure, max_index, output_dir=".", fmt="vasp",
                       min_slab_size=10.0, min_vacuum_size=15.0,
                       center_slab=True, symmetrize=False,
                       ftol=0.1, max_broken_bonds=0, bonds=None,
                       fix_bottom=0.5,
                       require_symmetric=False,
                       require_stoichiometric=False):
    """Generate most stable termination for all distinct surfaces up to max_index."""
    from pymatgen.core import Composition

    formula = Composition(structure.composition).reduced_formula
    indices = get_distinct_miller_indices(structure, max_index)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nGenerating surfaces for {formula} (max Miller index = {max_index})")
    print(f"Found {len(indices)} symmetrically distinct Miller indices\n")

    all_written = []
    for mi in sorted(indices, key=lambda x: sum(abs(i) for i in x)):
        slabs = generate_slabs_for_miller(
            structure, mi,
            min_slab_size=min_slab_size,
            min_vacuum_size=min_vacuum_size,
            center_slab=center_slab,
            symmetrize=symmetrize,
            ftol=ftol,
            max_broken_bonds=max_broken_bonds,
            bonds=bonds,
        )

        if not slabs:
            print(f"{miller_display(mi)}: no valid slabs (all break bonds)")
            continue

        n_total = len(slabs)
        slabs, dropped = filter_sym_stoich(
            slabs, structure,
            require_symmetric=require_symmetric,
            require_stoichiometric=require_stoichiometric,
        )
        if not slabs:
            print(f"{miller_display(mi)}: all {n_total} terminations rejected by filter")
            continue
        if dropped:
            print(f"{miller_display(mi)}: dropped {len(dropped)}/{n_total} "
                  f"terminations by filter")

        slab = slabs[0]
        n_fixed = 0
        if fix_bottom > 0:
            slab, n_fixed = apply_selective_dynamics(slab, fix_bottom)
        fname = slab_filename(formula, mi, 0, fmt)
        out_path = output_dir / fname
        write_slab(slab, out_path, fmt)

        sym = "yes" if slab.is_symmetric() else "no"
        fixed_str = f", fixed={n_fixed}/{slab.num_sites}" if n_fixed > 0 else ""
        print(f"{miller_display(mi)}: {slab.num_sites} atoms{fixed_str}, "
              f"area={slab.surface_area:.2f} Å², "
              f"symmetric={sym}, "
              f"{len(slabs)} termination(s) → {out_path}")
        all_written.append((slab, out_path))

    print(f"\nGenerated {len(all_written)} slab files in {output_dir}/")
    return all_written


def main():
    parser = argparse.ArgumentParser(
        description="Generate surface slabs from bulk structures"
    )
    parser.add_argument("--input", "-i", required=True,
                        help="Input bulk structure file (POSCAR, CIF, XYZ, JSON)")
    parser.add_argument("--miller", "-m", nargs=3, type=int,
                        help="Miller indices h k l (e.g., 1 1 1)")
    parser.add_argument("--max-index", type=int,
                        help="Max Miller index for enumerating all surfaces")
    parser.add_argument("--min-slab-size", type=float, default=10.0,
                        help="Minimum slab thickness in Å (default: 10.0)")
    parser.add_argument("--min-vacuum-size", type=float, default=15.0,
                        help="Minimum vacuum thickness in Å (default: 15.0)")
    parser.add_argument("--all-terminations", action="store_true",
                        help="Generate all unique terminations")
    parser.add_argument("--symmetrize", action="store_true",
                        help="Force symmetric slabs (top = bottom)")
    parser.add_argument("--center-slab", action="store_true", default=True,
                        help="Center slab in vacuum (default: True)")
    parser.add_argument("--output", "-o", help="Output filename (single slab)")
    parser.add_argument("--output-dir", default=".",
                        help="Output directory (multiple slabs)")
    parser.add_argument("--format", "-f", default="vasp",
                        choices=["vasp", "cif", "xyz", "json"],
                        help="Output format (default: vasp)")
    parser.add_argument("--list-surfaces", action="store_true",
                        help="List surfaces and exit (no files written)")
    parser.add_argument("--ftol", type=float, default=0.1,
                        help="Tolerance for termination clustering in Å (default: 0.1)")
    parser.add_argument("--max-broken-bonds", type=int, default=0,
                        help="Max broken bonds when cleaving (default: 0)")
    parser.add_argument("--bonds", nargs="+",
                        help="Bond pairs to preserve, e.g. Ti-O:2.0 Si-O:1.8")
    parser.add_argument("--fix-bottom", type=float, default=0.5,
                        help="Fraction of atoms to fix at bottom (default: 0.5). "
                             "Set to 0 to disable constraints.")
    parser.add_argument("--require-symmetric", action="store_true",
                        help="Drop terminations whose slab is not symmetric "
                             "(top != bottom). Affects --all-terminations and "
                             "--max-index modes.")
    parser.add_argument("--require-stoichiometric", action="store_true",
                        help="Drop terminations whose composition does not match "
                             "bulk stoichiometry. Pairs naturally with "
                             "--require-symmetric for clean Wulff calculations.")

    args = parser.parse_args()

    # Validate arguments
    if not args.miller and not args.max_index:
        print("ERROR: Provide --miller h k l or --max-index N", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    structure = load_structure(args.input)
    bonds = parse_bonds(args.bonds)

    # List surfaces mode
    if args.list_surfaces:
        if not args.max_index:
            print("ERROR: --list-surfaces requires --max-index", file=sys.stderr)
            sys.exit(1)
        list_surfaces(
            structure, args.max_index,
            min_slab_size=args.min_slab_size,
            min_vacuum_size=args.min_vacuum_size,
            ftol=args.ftol,
            bonds=bonds,
            max_broken_bonds=args.max_broken_bonds,
        )
        return

    # Max index mode
    if args.max_index:
        process_max_index(
            structure, args.max_index,
            output_dir=args.output_dir,
            fmt=args.format,
            min_slab_size=args.min_slab_size,
            min_vacuum_size=args.min_vacuum_size,
            center_slab=args.center_slab,
            symmetrize=args.symmetrize,
            ftol=args.ftol,
            max_broken_bonds=args.max_broken_bonds,
            bonds=bonds,
            fix_bottom=args.fix_bottom,
            require_symmetric=args.require_symmetric,
            require_stoichiometric=args.require_stoichiometric,
        )
        return

    # Single Miller index mode
    miller_index = tuple(args.miller)
    process_single_miller(
        structure, miller_index,
        all_terminations=args.all_terminations,
        output=args.output,
        output_dir=args.output_dir,
        fmt=args.format,
        min_slab_size=args.min_slab_size,
        min_vacuum_size=args.min_vacuum_size,
        center_slab=args.center_slab,
        symmetrize=args.symmetrize,
        ftol=args.ftol,
        max_broken_bonds=args.max_broken_bonds,
        bonds=bonds,
        fix_bottom=args.fix_bottom,
        require_symmetric=args.require_symmetric,
        require_stoichiometric=args.require_stoichiometric,
    )


if __name__ == "__main__":
    main()
