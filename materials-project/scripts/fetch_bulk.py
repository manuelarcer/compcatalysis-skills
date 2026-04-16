#!/usr/bin/env python3
"""Fetch bulk crystal structures from the Materials Project database.

Requires:
    pip install mp-api pymatgen

Usage:
    python fetch_bulk.py --formula Pt --output Pt_bulk.vasp
    python fetch_bulk.py --mp-id mp-126 --output Pt_bulk.vasp
    python fetch_bulk.py --formula Pt --conventional --output Pt_bulk_conv.vasp
    python fetch_bulk.py --formula TiO2 --list-polymorphs
    python fetch_bulk.py --formula Pt Cu Ni --output-dir bulk_structures/
"""

import argparse
import json
import os
import sys
from pathlib import Path


def get_api_key(cli_key=None):
    """Get MP API key from CLI arg or environment variable."""
    key = cli_key or os.environ.get("MP_API_KEY")
    if not key:
        print("ERROR: No Materials Project API key found.", file=sys.stderr)
        print("", file=sys.stderr)
        print("Set it as an environment variable:", file=sys.stderr)
        print("  export MP_API_KEY='your_key_here'", file=sys.stderr)
        print("", file=sys.stderr)
        print("Or pass it directly:", file=sys.stderr)
        print("  python fetch_bulk.py --api-key YOUR_KEY ...", file=sys.stderr)
        print("", file=sys.stderr)
        print("Get a key at: https://next-gen.materialsproject.org/api", file=sys.stderr)
        sys.exit(1)
    return key


def fetch_by_mpid(api_key, mp_id):
    """Fetch a single structure by Materials Project ID."""
    from mp_api.client import MPRester

    with MPRester(api_key) as mpr:
        docs = mpr.materials.summary.search(
            material_ids=[mp_id],
            fields=[
                "material_id",
                "formula_pretty",
                "structure",
                "symmetry",
                "energy_above_hull",
                "formation_energy_per_atom",
                "band_gap",
                "nsites",
            ],
        )
        if not docs:
            print(f"No structure found for '{mp_id}'", file=sys.stderr)
            sys.exit(1)
        return docs[0]


def fetch_by_formula(api_key, formula, max_ehull=0.025):
    """Fetch structures matching a chemical formula, sorted by stability."""
    from mp_api.client import MPRester

    with MPRester(api_key) as mpr:
        docs = mpr.materials.summary.search(
            formula=formula,
            energy_above_hull=(0, max_ehull),
            fields=[
                "material_id",
                "formula_pretty",
                "structure",
                "symmetry",
                "energy_above_hull",
                "formation_energy_per_atom",
                "band_gap",
                "nsites",
            ],
        )
    return sorted(docs, key=lambda d: d.energy_above_hull)


def list_polymorphs(api_key, formula, max_ehull=1.0):
    """List all polymorphs for a given formula."""
    from mp_api.client import MPRester

    with MPRester(api_key) as mpr:
        docs = mpr.materials.summary.search(
            formula=formula,
            energy_above_hull=(0, max_ehull),
            fields=[
                "material_id",
                "formula_pretty",
                "symmetry",
                "energy_above_hull",
                "formation_energy_per_atom",
                "nsites",
            ],
        )

    docs = sorted(docs, key=lambda d: d.energy_above_hull)

    if not docs:
        print(f"No structures found for formula '{formula}' (max Ehull = {max_ehull} eV/atom)")
        return

    print(f"\nPolymorphs for {formula} ({len(docs)} found):")
    print(f"{'Material ID':<15} {'Space Group':<20} {'Ehull (eV/at)':<15} {'Sites':<8} {'Form. E (eV/at)':<15}")
    print("-" * 75)
    for doc in docs:
        sg = doc.symmetry.symbol if doc.symmetry else "?"
        print(
            f"{str(doc.material_id):<15} "
            f"{sg:<20} "
            f"{doc.energy_above_hull:<15.4f} "
            f"{doc.nsites:<8} "
            f"{doc.formation_energy_per_atom:<15.4f}"
        )


def to_conventional(structure):
    """Convert to conventional standard cell."""
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    sga = SpacegroupAnalyzer(structure)
    return sga.get_conventional_standard_structure()


def write_structure(structure, output_path, fmt="vasp"):
    """Write structure to file in the specified format."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if fmt == "vasp":
        structure.to(filename=str(output_path), fmt="poscar")
    elif fmt == "cif":
        structure.to(filename=str(output_path), fmt="cif")
    elif fmt == "xyz":
        structure.to(filename=str(output_path), fmt="xyz")
    elif fmt == "json":
        with open(output_path, "w") as f:
            json.dump(structure.as_dict(), f, indent=2)
    else:
        raise ValueError(f"Unknown format: {fmt}. Use: vasp, cif, xyz, json")


def print_summary(doc, structure):
    """Print a summary of the fetched structure."""
    lattice = structure.lattice
    sg = doc.symmetry.symbol if doc.symmetry else "?"
    sg_num = doc.symmetry.number if doc.symmetry else "?"

    print(f"\n--- Structure Summary ---")
    print(f"Material ID:    {doc.material_id}")
    print(f"Formula:        {doc.formula_pretty}")
    print(f"Space group:    {sg} (#{sg_num})")
    print(f"Lattice (Å):    a={lattice.a:.4f}  b={lattice.b:.4f}  c={lattice.c:.4f}")
    print(f"Angles (°):     α={lattice.alpha:.2f}  β={lattice.beta:.2f}  γ={lattice.gamma:.2f}")
    print(f"Volume (ų):   {lattice.volume:.2f}")
    print(f"Num atoms:      {structure.num_sites}")
    print(f"Ehull (eV/at):  {doc.energy_above_hull:.4f}")
    print(f"Form E (eV/at): {doc.formation_energy_per_atom:.4f}")
    if doc.band_gap is not None:
        print(f"Band gap (eV):  {doc.band_gap:.3f}")
    print(f"------------------------")


def process_single(api_key, formula=None, mp_id=None, max_ehull=0.025,
                    conventional=False, output=None, fmt="vasp"):
    """Fetch and write a single structure."""
    if mp_id:
        doc = fetch_by_mpid(api_key, mp_id)
        structure = doc.structure
    elif formula:
        docs = fetch_by_formula(api_key, formula, max_ehull)
        if not docs:
            print(f"No structures found for '{formula}' with Ehull ≤ {max_ehull} eV/atom",
                  file=sys.stderr)
            print("Try: --max-ehull 0.1 or --list-polymorphs", file=sys.stderr)
            sys.exit(1)

        doc = docs[0]
        structure = doc.structure

        if len(docs) > 1:
            print(f"Note: {len(docs)} polymorphs found. Using most stable ({doc.material_id}).")
            print(f"  Run with --list-polymorphs to see all options.")
    else:
        print("ERROR: Provide --formula or --mp-id", file=sys.stderr)
        sys.exit(1)

    if conventional:
        structure = to_conventional(structure)
        print(f"Converted to conventional cell ({structure.num_sites} atoms)")

    # Default output filename
    if output is None:
        ext = {"vasp": ".vasp", "cif": ".cif", "xyz": ".xyz", "json": ".json"}
        output = f"{doc.formula_pretty}_bulk{ext.get(fmt, '.vasp')}"

    write_structure(structure, output, fmt)
    print_summary(doc, structure)
    print(f"\nSaved to: {output}")

    return doc, structure


def main():
    parser = argparse.ArgumentParser(
        description="Fetch bulk structures from Materials Project"
    )
    parser.add_argument("--formula", nargs="+", help="Chemical formula(s)")
    parser.add_argument("--mp-id", help="Materials Project ID (e.g., mp-126)")
    parser.add_argument("--output", "-o", help="Output filename")
    parser.add_argument("--output-dir", default=".", help="Output directory (for multiple formulas)")
    parser.add_argument("--format", "-f", default="vasp",
                        choices=["vasp", "cif", "xyz", "json"],
                        help="Output format (default: vasp)")
    parser.add_argument("--conventional", action="store_true",
                        help="Use conventional standard cell")
    parser.add_argument("--max-ehull", type=float, default=0.025,
                        help="Max energy above hull in eV/atom (default: 0.025)")
    parser.add_argument("--list-polymorphs", action="store_true",
                        help="List polymorphs and exit")
    parser.add_argument("--api-key", help="MP API key (overrides MP_API_KEY env var)")

    args = parser.parse_args()

    api_key = get_api_key(args.api_key)

    # List polymorphs mode
    if args.list_polymorphs:
        if not args.formula:
            print("ERROR: --list-polymorphs requires --formula", file=sys.stderr)
            sys.exit(1)
        for formula in args.formula:
            list_polymorphs(api_key, formula)
        return

    # Fetch by mp-id
    if args.mp_id:
        process_single(api_key, mp_id=args.mp_id, conventional=args.conventional,
                        output=args.output, fmt=args.format)
        return

    # Fetch by formula(s)
    if not args.formula:
        parser.print_help()
        sys.exit(1)

    if len(args.formula) == 1:
        process_single(api_key, formula=args.formula[0], max_ehull=args.max_ehull,
                        conventional=args.conventional, output=args.output, fmt=args.format)
    else:
        # Multiple formulas → output to directory
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        ext = {"vasp": ".vasp", "cif": ".cif", "xyz": ".xyz", "json": ".json"}

        for formula in args.formula:
            out_file = output_dir / f"{formula}_bulk{ext.get(args.format, '.vasp')}"
            print(f"\n{'='*50}")
            print(f"Fetching: {formula}")
            print(f"{'='*50}")
            try:
                process_single(api_key, formula=formula, max_ehull=args.max_ehull,
                                conventional=args.conventional, output=str(out_file),
                                fmt=args.format)
            except SystemExit:
                print(f"Skipping {formula} (not found)", file=sys.stderr)
                continue


if __name__ == "__main__":
    main()
