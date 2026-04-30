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
hydrogen electrode you usually need a **derived reference**:

| Adsorbate | Reference | Builder |
|---|---|---|
| O\* | ½O₂ | `--coeff 0.5 --ads-ref O2_gas/opt.traj` |
| OH\* | H₂O − ½H₂ | compute manually: pass H₂O as ref with coeff 1.0, then subtract 0.5·E(H₂) yourself, OR use a config with two passes |
| OOH\* | 2·H₂O − 1.5·H₂ | derived; can't be expressed with a single `(coeff, ref)` pair — use external math or a custom reference traj |
| CO\* | CO | `--coeff 1.0 --ads-ref CO_gas/opt.traj` |

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

`eads.json`:

```json
{
  "slab": "clean/opt.traj",
  "calculations": [
    {"label": "CO_top_0", "slab_ads": "CO/Pt_top_0/opt.traj",
     "ref": "CO_gas/opt.traj"},
    {"label": "CO_top_1", "slab_ads": "CO/Pt_top_1/opt.traj",
     "ref": "CO_gas/opt.traj"},
    {"label": "O_fcc",    "slab_ads": "O/fcc/opt.traj",
     "ref": "O2_gas/opt.traj", "coeff": 0.5}
  ]
}
```

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

`--config` is mutually exclusive with the per-trajectory flags.

#### Output CSV columns

`label, slab, slab_ads, ads_ref, coeff, E_slab_eV, E_slab_ads_eV, E_ref_eV, E_ads_eV`

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
