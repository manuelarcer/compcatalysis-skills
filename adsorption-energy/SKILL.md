---
name: adsorption-energy
description: >
  Compute adsorption energies from relaxed trajectories produced by
  MLIP optimization. Single-shot or batch via JSON config.
  Triggers: "adsorption energy", "compute Eads", "binding energy",
  "rank adsorption sites", "compute E_ads from trajectories",
  "OER intermediate energies".
  Use when: user has relaxed trajectories of a clean slab, slab+adsorbate,
  and a gas-phase reference, and wants the adsorption energy E_ads.
  NOT for: running the relaxations (use relaxation-orchestration), NEB
  barriers (use neb-search), or zero-point / free-energy corrections.
metadata: {"requires": {"packages": ["ase"], "python": "3.11+"}}
---

# Adsorption Energy

Compute adsorption energies from existing relaxed trajectories. This is a
**post-processing skill** — it reads `opt.traj` files and applies the
adsorption-energy formula. It does not run any MLIP.

## Formula

```
E_ads(X) = E(slab + X) - E(slab) - coeff · E(X_ref)
```

with `coeff = 1.0` for whole-molecule references (CO, CO₂, H₂O, NH₃) and
`coeff = 0.5` when referencing half of a diatomic (O ↔ ½O₂, H ↔ ½H₂).

For OER intermediates referenced to H₂O / H₂ via the computational
hydrogen electrode (CHE), use the **named composite-reference** form in
the JSON config (see "Batch mode" below). One-shot CLI mode supports only
single-traj references with a scalar coefficient.

| Adsorbate | Reference (CHE) | Form |
|---|---|---|
| O\* | ½ O₂ | single-ref, `--coeff 0.5 --ads-ref O2_gas/opt.traj` |
| OH\* | H₂O − ½ H₂ | composite, JSON config only |
| OOH\* | 2·H₂O − 1.5·H₂ | composite, JSON config only |
| CO\* | CO | single-ref, `--coeff 1.0 --ads-ref CO_gas/opt.traj` |

Negative `E_ads` = exothermic adsorption.

## Prerequisites

- `ase` installed.
- Trajectories from `relaxation-orchestration` (or any ASE-readable file
  that contains a `SinglePointCalculator` with `potential_energy`):
  - one clean-slab traj
  - one slab+adsorbate traj per configuration
  - one gas-phase reference traj per adsorbate species
- All trajectories must come from runs with the **same MLIP and same
  `--uma-task`**. Mixing task heads gives meaningless `E_ads` values.

## Scripts

### compute_eads.py — Single or batch adsorption-energy calculation

**Location:** `scripts/compute_eads.py`

#### Single mode

```bash
python compute_eads.py \
    --slab clean/opt.traj \
    --slab-ads CO_top_0/opt.traj \
    --ads-ref CO_gas/opt.traj \
    --label CO_top_0 \
    --output eads.csv

# O* with ½O2 reference, appended to the same CSV
python compute_eads.py \
    --slab clean/opt.traj \
    --slab-ads O_fcc_0/opt.traj \
    --ads-ref O2_gas/opt.traj \
    --coeff 0.5 --label O_fcc_0 \
    --output eads.csv --append
```

#### Batch mode (JSON config)

```bash
python compute_eads.py --config eads.json --output eads.csv
```

`eads.json` (with named composite references for OER intermediates):

```json
{
  "slab": "clean/opt.traj",
  "references": {
    "OH":  [{"traj": "gas/H2O/opt.traj", "coeff":  1.0},
            {"traj": "gas/H2/opt.traj",  "coeff": -0.5}],
    "OOH": [{"traj": "gas/H2O/opt.traj", "coeff":  2.0},
            {"traj": "gas/H2/opt.traj",  "coeff": -1.5}],
    "O":   [{"traj": "gas/O2/opt.traj",  "coeff":  0.5}]
  },
  "calculations": [
    {"label": "CO_top_0", "slab_ads": "CO/Pt_top_0/opt.traj",
     "ref": "gas/CO/opt.traj"},
    {"label": "OH_top_0", "slab_ads": "OH/Pt_top_0/opt.traj",
     "ref": "OH"},
    {"label": "OOH_top_0", "slab_ads": "OOH/Pt_top_0/opt.traj",
     "ref": "OOH"},
    {"label": "O_fcc",    "slab_ads": "O/fcc/opt.traj",
     "ref": "O"}
  ]
}
```

The `ref` field is resolved as:
- if it matches a key in `references`, expand to that composite formula;
- otherwise treat as a path to a single trajectory (back-compat with
  pre-composite configs).

#### Arguments

| Argument | Description | Default |
|---|---|---|
| `--config` | JSON config for batch mode | — |
| `--slab` | Clean slab trajectory | — |
| `--slab-ads` | Slab+adsorbate trajectory | — |
| `--ads-ref` | Gas-phase reference trajectory | — |
| `--label` | Row label in output CSV | `adsorption` |
| `--coeff` | Stoichiometric coefficient on the reference | `1.0` |
| `--output` | Output CSV path | `eads.csv` |
| `--append` | Append rather than overwrite | off |
| `--no-consistency-check` | Skip the cross-trajectory MLIP+task-head consistency check | off |

`--config` is mutually exclusive with the per-trajectory flags.

#### Output CSV columns

`label, slab, slab_ads, ads_ref, ref_summary, coeff, E_slab_eV, E_slab_ads_eV, E_ref_eV, E_ads_eV`

`coeff` is the scalar for single-ref calculations and the literal string
`composite` for composite-reference calculations; in the composite case the
human-readable formula is in `ref_summary` (e.g. `1·H2O + -0.5·H2`).

#### Task-head consistency check

The script reads the `compcat_run.json` sidecar produced by
`relaxation-orchestration` next to each input trajectory and refuses if the
MLIP model or UMA task head differs across inputs. This catches the silent
"clean slab relaxed with `omat`, slab+ads relaxed with `oc20`" mistake at
the boundary instead of letting it produce a plausible-but-wrong energy.

If your trajectories predate the sidecar convention (or come from a
non-`mlip_platform` source), the check warns rather than failing. To skip
it entirely, pass `--no-consistency-check`.

## Pitfalls

- **Task-head mismatch invalidates E_ads.** A clean slab relaxed with
  `omat` and a slab+ads relaxed with `oc20` give individually-reasonable
  energies but a wrong adsorption energy. Verify all source trajectories
  came from the same `--uma-task`.
- **Selective dynamics inconsistencies.** If the clean slab was free to
  relax but the slab+ads had bottom layers fixed (or vice versa), the
  bulk reference is not consistent. Use the same constraints in both.
- **Reactive co-adsorbate placements.** If the slab+ads relaxation
  spontaneously dissociated or migrated, `E_ads` reflects the final
  product, not the named adsorbate. Inspect `opt_final.vasp` before
  trusting the energy. (See adsorbate-placement pitfalls for spacing
  guidance to avoid this.)
- **Half-reference for atomic adsorbates only.** Use `--coeff 0.5` only
  when the adsorbate is half of a diatomic in the gas phase (O*/½O₂,
  H*/½H₂). Don't use it for "half of a molecule" of CO or H₂O.

## Agent Workflow Integration

This skill is **step 7** in the catalysis pipeline (post-relaxation):

1. Fetch bulk
2. Optimize bulk
3. Generate surfaces
4. Optimize clean slabs
5. Place adsorbates
6. Optimize slab+adsorbate (via relaxation-orchestration)
7. **Compute adsorption energies** (this skill)

For a typical OER site-screening run, the agent should:

1. Walk the `placements/` tree to find every `opt.traj` produced by
   relaxation-orchestration.
2. Group by adsorbate (subdir name).
3. Build an `eads.json` mapping each `slab_ads` to the matching gas-phase
   reference, with `coeff=0.5` for atomic O.
4. Invoke `compute_eads.py --config eads.json --output eads.csv`.
5. Sort the CSV by `E_ads_eV` to rank sites.
