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

RUN_META_FILENAME = "compcat_run.json"


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


def _ref_label(path: Path) -> str:
    """Pick a readable label for a trajectory: when the filename is generic
    (`opt.traj`, `opt_final.vasp`) use the parent directory's name; otherwise
    use the filename itself."""
    if path.name in ("opt.traj", "opt_final.vasp"):
        return path.parent.name or path.name
    return path.name


def read_run_metadata(traj: Path) -> dict | None:
    """Return the {mlip, uma_task, fmax, optimizer} dict written by
    relaxation-orchestration alongside opt_final.vasp, or None if absent."""
    meta = traj.parent / RUN_META_FILENAME
    if not meta.exists():
        return None
    try:
        return json.loads(meta.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def check_run_consistency(trajs: list[Path], strict: bool = True
                           ) -> tuple[bool, str]:
    """Verify all input trajectories were produced with the same MLIP and
    task head. Returns (ok, message). If any trajectory has no
    compcat_run.json sidecar, treat as 'unknown' and warn (or fail in
    strict mode)."""
    metas = [(t, read_run_metadata(t)) for t in trajs]
    missing = [str(t) for t, m in metas if m is None]
    if missing and strict:
        return (False, "missing compcat_run.json for: " + ", ".join(missing))

    keys = {(m["mlip"], m.get("uma_task", "?"))
            for _, m in metas if m is not None}
    if len(keys) > 1:
        joined = "; ".join(f"{m}|{t}" for m, t in keys)
        return (False, f"inconsistent (mlip, uma_task) across inputs: {joined}")
    if missing:
        return (True, "warning: missing run metadata for "
                       + ", ".join(missing) + " (consistency only partially verified)")
    return (True, "all inputs share (mlip, uma_task)")


def expand_reference(ref_spec, ref_table: dict) -> list[tuple[Path, float]]:
    """Resolve a `ref` field from the JSON config to a list of
    (path, coeff) terms.

    - If `ref_spec` is a list of {"traj": <path>, "coeff": <float>}, use directly.
    - If it's a string and matches a key in ref_table, expand from the table.
    - Else treat the string as a path with coeff=1.0 (back-compat with the
      original single-ref form).
    """
    if isinstance(ref_spec, list):
        return [(Path(t["traj"]), float(t["coeff"])) for t in ref_spec]
    if isinstance(ref_spec, str):
        if ref_spec in ref_table:
            entry = ref_table[ref_spec]
            return [(Path(t["traj"]), float(t["coeff"])) for t in entry]
        return [(Path(ref_spec), 1.0)]
    raise ValueError(f"unsupported reference spec: {ref_spec!r}")


def compute_one(slab: Path, slab_ads: Path, ads_ref: Path,
                label: str, coeff: float = 1.0,
                consistency: bool = True) -> dict:
    """Compute E_ads against a single gas-phase reference trajectory."""
    if consistency:
        ok, msg = check_run_consistency([slab, slab_ads, ads_ref], strict=False)
        if not ok:
            raise ValueError(f"task-head consistency check failed: {msg}")
    e_slab = read_final_energy(slab)
    e_slab_ads = read_final_energy(slab_ads)
    e_ref = read_final_energy(ads_ref)
    e_ads = adsorption_energy(e_slab, e_slab_ads, e_ref, coeff)
    return {
        "label": label,
        "slab": str(slab),
        "slab_ads": str(slab_ads),
        "ads_ref": str(ads_ref),
        "ref_summary": f"{coeff:g} × {_ref_label(ads_ref)}",
        "coeff": coeff,
        "E_slab_eV": round(e_slab, 6),
        "E_slab_ads_eV": round(e_slab_ads, 6),
        "E_ref_eV": round(e_ref, 6),
        "E_ads_eV": round(e_ads, 6),
    }


def compute_one_composite(slab: Path, slab_ads: Path,
                           ref_terms: list[tuple[Path, float]],
                           label: str, consistency: bool = True) -> dict:
    """Compute E_ads against a linear combination of gas-phase references.

    For OER intermediates referenced via the computational hydrogen electrode:
        OH:  [(H2O, 1.0), (H2, -0.5)]
        OOH: [(H2O, 2.0), (H2, -1.5)]
        O:   [(O2,  0.5)]
    """
    all_trajs = [slab, slab_ads, *(t for t, _ in ref_terms)]
    if consistency:
        ok, msg = check_run_consistency(all_trajs, strict=False)
        if not ok:
            raise ValueError(f"task-head consistency check failed: {msg}")

    e_slab = read_final_energy(slab)
    e_slab_ads = read_final_energy(slab_ads)
    e_ref_total = sum(coeff * read_final_energy(t) for t, coeff in ref_terms)
    e_ads = e_slab_ads - e_slab - e_ref_total
    summary = " + ".join(f"{c:g}·{_ref_label(t)}" for t, c in ref_terms)
    return {
        "label": label,
        "slab": str(slab),
        "slab_ads": str(slab_ads),
        "ads_ref": "; ".join(str(t) for t, _ in ref_terms),
        "ref_summary": summary,
        "coeff": "composite",
        "E_slab_eV": round(e_slab, 6),
        "E_slab_ads_eV": round(e_slab_ads, 6),
        "E_ref_eV": round(e_ref_total, 6),
        "E_ads_eV": round(e_ads, 6),
    }


def compute_from_config(config: dict, consistency: bool = True) -> list[dict]:
    """Run all calculations declared in a JSON config.

    Schema:
        {
          "slab": "<path>",
          "references": {                # optional; named composite refs
            "OH":  [{"traj": "H2O.traj", "coeff": 1.0},
                    {"traj": "H2.traj",  "coeff": -0.5}],
            "OOH": [{"traj": "H2O.traj", "coeff": 2.0},
                    {"traj": "H2.traj",  "coeff": -1.5}],
            "O":   [{"traj": "O2.traj",  "coeff": 0.5}]
          },
          "calculations": [
            # `ref` is either:
            #   - a key from the references table (composite),
            #   - a single path with optional coeff (single-ref form).
            {"label": "OH_top_0", "slab_ads": "OH/Pt_top_0/opt.traj",
             "ref": "OH"},
            {"label": "CO_top_0", "slab_ads": "CO/Pt_top_0/opt.traj",
             "ref": "CO_gas/opt.traj", "coeff": 1.0}
          ]
        }
    """
    if "slab" not in config or "calculations" not in config:
        raise ValueError("config must have 'slab' and 'calculations' keys")
    slab = Path(config["slab"])
    ref_table = config.get("references", {})
    rows = []
    for entry in config["calculations"]:
        terms = expand_reference(entry["ref"], ref_table)
        # Single term with coeff=1.0 from a path → use simple form so the
        # output row stays compatible with single-mode CSV
        if len(terms) == 1 and entry.get("coeff", 1.0) != 1.0:
            ref_path, _ = terms[0]
            rows.append(compute_one(
                slab=slab, slab_ads=Path(entry["slab_ads"]),
                ads_ref=ref_path, label=entry["label"],
                coeff=float(entry["coeff"]),
                consistency=consistency,
            ))
        elif len(terms) == 1:
            ref_path, coeff = terms[0]
            rows.append(compute_one(
                slab=slab, slab_ads=Path(entry["slab_ads"]),
                ads_ref=ref_path, label=entry["label"],
                coeff=coeff, consistency=consistency,
            ))
        else:
            rows.append(compute_one_composite(
                slab=slab, slab_ads=Path(entry["slab_ads"]),
                ref_terms=terms, label=entry["label"],
                consistency=consistency,
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
        coeff_str = r['coeff'] if isinstance(r['coeff'], str) else f"{r['coeff']:.2f}"
        print(f"{r['label']:<20} {r['E_ads_eV']:<12.4f} {coeff_str:<7} {r['slab_ads']}")


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
    parser.add_argument("--no-consistency-check", action="store_true",
                        help="Skip the task-head consistency check across inputs "
                             "(use only when the inputs predate the run-metadata convention)")
    args = parser.parse_args()
    consistency = not args.no_consistency_check

    if args.config:
        config = json.loads(args.config.read_text())
        rows = compute_from_config(config, consistency=consistency)
    else:
        missing = [n for n in ("slab", "slab_ads", "ads_ref")
                   if getattr(args, n) is None]
        if missing:
            parser.error(f"missing required args: {missing} (or use --config)")
        rows = [compute_one(args.slab, args.slab_ads, args.ads_ref,
                             args.label, args.coeff, consistency=consistency)]

    write_csv(rows, args.output, append=args.append)
    print_table(rows)
    print(f"\nWrote {len(rows)} row(s) to {args.output}")


if __name__ == "__main__":
    main()
