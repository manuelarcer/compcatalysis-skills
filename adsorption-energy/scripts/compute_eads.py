#!/usr/bin/env python3
"""Compute adsorption energies from relaxed structure trajectories.

E_ads(X) = E(slab+X) - E(slab) - coeff * E(X_ref)

For OER intermediates the gas-phase references are typically:
  E_ads(O)   = E(slab+O)   - E(slab) - 0.5 * E(O2)        (--coeff 0.5)
  E_ads(OH)  = E(slab+OH)  - E(slab) - (E(H2O) - 0.5*E(H2))   # build manually
  E_ads(CO)  = E(slab+CO)  - E(slab) - E(CO)
  E_ads(OOH) = E(slab+OOH) - E(slab) - (2*E(H2O) - 1.5*E(H2))  # build manually

Usage:
    # Single calculation
    python compute_eads.py \\
        --slab clean/opt.traj \\
        --slab-ads CO_top_0/opt.traj \\
        --ads-ref CO_gas/opt.traj \\
        --label CO_top_0 \\
        --output eads.csv

    # Half-O2 reference for O*
    python compute_eads.py \\
        --slab clean/opt.traj \\
        --slab-ads O_fcc_0/opt.traj \\
        --ads-ref O2_gas/opt.traj \\
        --coeff 0.5 --label O_fcc_0 \\
        --output eads.csv --append

    # Batch via JSON config
    python compute_eads.py --config eads.json --output eads.csv
"""

import argparse
import csv
import json
import sys
from pathlib import Path


def read_final_energy(path: Path) -> float:
    """Read the final-frame potential energy from an ASE trajectory."""
    from ase.io import read

    if not path.exists():
        raise FileNotFoundError(f"trajectory not found: {path}")
    atoms = read(str(path), index=-1)
    return float(atoms.get_potential_energy())


def adsorption_energy(e_slab: float, e_slab_ads: float, e_ref: float,
                      coeff: float = 1.0) -> float:
    """E_ads = E(slab+ads) - E(slab) - coeff * E(ref). Negative = exothermic."""
    return e_slab_ads - e_slab - coeff * e_ref


def compute_one(slab: Path, slab_ads: Path, ads_ref: Path,
                label: str, coeff: float = 1.0) -> dict:
    e_slab = read_final_energy(slab)
    e_slab_ads = read_final_energy(slab_ads)
    e_ref = read_final_energy(ads_ref)
    e_ads = adsorption_energy(e_slab, e_slab_ads, e_ref, coeff)
    return {
        "label": label,
        "slab": str(slab),
        "slab_ads": str(slab_ads),
        "ads_ref": str(ads_ref),
        "coeff": coeff,
        "E_slab_eV": round(e_slab, 6),
        "E_slab_ads_eV": round(e_slab_ads, 6),
        "E_ref_eV": round(e_ref, 6),
        "E_ads_eV": round(e_ads, 6),
    }


def compute_from_config(config: dict) -> list[dict]:
    """Run all calculations declared in a JSON config.

    Schema:
        {
          "slab": "<path>",                              # shared clean-slab traj
          "calculations": [
            {"label": str, "slab_ads": "<path>",
             "ref": "<path>", "coeff": float (default 1.0)},
            ...
          ]
        }
    """
    if "slab" not in config or "calculations" not in config:
        raise ValueError("config must have 'slab' and 'calculations' keys")
    slab = Path(config["slab"])
    rows = []
    for entry in config["calculations"]:
        rows.append(compute_one(
            slab=slab,
            slab_ads=Path(entry["slab_ads"]),
            ads_ref=Path(entry["ref"]),
            label=entry["label"],
            coeff=float(entry.get("coeff", 1.0)),
        ))
    return rows


def write_csv(rows: list[dict], path: Path, append: bool = False) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    write_header = not (append and path.exists())
    mode = "a" if append and path.exists() else "w"
    with open(path, mode, newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerows(rows)


def print_table(rows: list[dict]) -> None:
    if not rows:
        return
    header = f"{'label':<20} {'E_ads (eV)':<12} {'coeff':<7} {'slab_ads'}"
    print(header)
    print("-" * len(header))
    for r in rows:
        print(f"{r['label']:<20} {r['E_ads_eV']:<12.4f} {r['coeff']:<7.2f} {r['slab_ads']}")


def main():
    parser = argparse.ArgumentParser(
        description="Compute adsorption energies from relaxed trajectories.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", type=Path,
                        help="JSON config for batch mode (overrides --slab/--slab-ads)")
    parser.add_argument("--slab", type=Path, help="Clean slab trajectory")
    parser.add_argument("--slab-ads", type=Path, help="Slab+adsorbate trajectory")
    parser.add_argument("--ads-ref", type=Path, help="Gas-phase reference trajectory")
    parser.add_argument("--label", default="adsorption", help="Label for the row")
    parser.add_argument("--coeff", type=float, default=1.0,
                        help="Stoichiometric coefficient on the reference (e.g. 0.5 for ½O2)")
    parser.add_argument("--output", type=Path, default=Path("eads.csv"),
                        help="Output CSV path (default: eads.csv)")
    parser.add_argument("--append", action="store_true",
                        help="Append to existing CSV instead of overwriting")
    args = parser.parse_args()

    if args.config:
        config = json.loads(args.config.read_text())
        rows = compute_from_config(config)
    else:
        missing = [n for n in ("slab", "slab_ads", "ads_ref")
                   if getattr(args, n) is None]
        if missing:
            parser.error(f"missing required args: {missing} (or use --config)")
        rows = [compute_one(args.slab, args.slab_ads, args.ads_ref,
                             args.label, args.coeff)]

    write_csv(rows, args.output, append=args.append)
    print_table(rows)
    print(f"\nWrote {len(rows)} row(s) to {args.output}")


if __name__ == "__main__":
    main()
