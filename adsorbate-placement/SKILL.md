---
name: adsorbate-placement
description: >
  Place adsorbate molecules on surface ontop sites of a slab.
  Triggers: "place adsorbate", "add adsorbate", "put OH on surface",
  "adsorb O on slab", "generate adsorption structures", "adsorbate placement",
  "place adsorbate on surface", "ontop sites", "adsorb molecule on slab".
  Use when: user has an optimized clean slab and wants to place adsorbates
  at surface ontop sites to generate input structures for relaxation.
  Supports built-in OER intermediates (O, OH, OOH) and arbitrary adsorbates from files.
  NOT for: slab generation, relaxation, or adsorption energy calculation.
metadata: {"requires": {"packages": ["ase", "numpy"], "python": "3.11+"}}
---

# Adsorbate Placement

Place adsorbate molecules on surface ontop sites of optimized slabs. Supports both built-in OER intermediates (O, OH, OOH) and arbitrary adsorbates loaded from structure files.

## Prerequisites

- Python packages: `ase`, `numpy`
- Input: an optimized clean slab structure file (VASP POSCAR, CIF, XYZ)
- Slab should have selective dynamics (FixAtoms) set from surface generation

## Scripts

### place_adsorbates.py — Place adsorbates on surface sites

**Location:** `scripts/place_adsorbates.py`

#### Usage

```bash
# Built-in OER adsorbates on CoOOH
python place_adsorbates.py --slab clean/t2/opt_final.vasp --adsorbates O OH OOH --sites Co O

# Arbitrary adsorbate from file — CO binding through C (atom index 0)
python place_adsorbates.py --slab slab.vasp --adsorbate-file CO.xyz --binding-atom 0 --sites Pt

# CCH binding through first C
python place_adsorbates.py --slab slab.vasp --adsorbate-file CCH.xyz --binding-atom 0 --sites Co

# Multiple file-based adsorbates
python place_adsorbates.py --slab slab.vasp --adsorbate-file CO.xyz HCOO.xyz --binding-atom 0 0 --sites Pt

# Mix built-in and file-based
python place_adsorbates.py --slab slab.vasp --adsorbates O OH --adsorbate-file CO.xyz --binding-atom 0 --sites Pt

# Override auto-computed height
python place_adsorbates.py --slab slab.vasp --adsorbates OH --sites Co --height 2.0

# Narrower surface depth for flat metal surfaces
python place_adsorbates.py --slab slab.vasp --adsorbates O --sites Pt --depth 1.0
```

#### Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--slab` | Input clean slab structure file | — (required) |
| `--adsorbates` | Built-in adsorbate(s): `O`, `OH`, `OOH` | — |
| `--adsorbate-file` | Adsorbate structure file(s) (XYZ, VASP, CIF) | — |
| `--binding-atom` | Index of binding atom in each adsorbate file (0-based) | `0` |
| `--sites` | Surface species to adsorb onto (e.g. `Co O Pt`) | — (required) |
| `--output-dir` | Root output directory | `.` |
| `--height` | Height above surface site (Å). Auto if omitted. | auto |
| `--depth` | Depth from top for surface detection (Å) | `2.0` |

At least one of `--adsorbates` or `--adsorbate-file` is required.

#### Adsorbate Input

**Built-in presets** (`--adsorbates`): O, OH, OOH with standard geometries. The binding atom is always O (appropriate for OER intermediates).

**Arbitrary adsorbates** (`--adsorbate-file`): Any molecule in a format ASE can read (XYZ, VASP, CIF, etc.). Use `--binding-atom` to specify which atom (0-based index) bonds to the surface. The molecule is recentered on the binding atom and oriented upward. If no `--binding-atom` is given, defaults to atom 0.

Examples of binding atom specification:
- **CO**: C binds to surface → `--binding-atom 0` (if C is first in XYZ)
- **CCH**: first C binds → `--binding-atom 0`
- **HCOO**: O binds → `--binding-atom 2` (depends on atom order in file)

#### Auto Height from Covalent Radii

When `--height` is not specified, the script computes the initial placement height as:

```
height = r_cov(surface_atom) + r_cov(binding_atom)
```

using ASE's covalent radii table. This gives physically reasonable starting geometries that the MLIP relaxation will refine. Examples:

| Surface | Binding atom | Auto height |
|---------|-------------|-------------|
| Co (1.26 Å) | O (0.66 Å) | 1.92 Å |
| Co (1.26 Å) | C (0.76 Å) | 2.02 Å |
| Pt (1.36 Å) | O (0.66 Å) | 2.02 Å |
| Pt (1.36 Å) | C (0.76 Å) | 2.12 Å |
| O (0.66 Å)  | O (0.66 Å) | 1.32 Å |

Use `--height` to override when you know a better value or if the auto-height produces atom overlaps.

#### Output Structure

```
output-dir/
  O/
    Co_top_0/input.vasp
    Co_top_1/input.vasp
    O_top_0/input.vasp
  OH/
    Co_top_0/input.vasp
    ...
  CO/                          # from --adsorbate-file CO.xyz
    Co_top_0/input.vasp
    ...
```

Each `input.vasp` is a complete slab+adsorbate structure ready for MLIP relaxation.

#### Site Detection

1. Finds all atoms within `--depth` Å of the topmost atom
2. Groups by chemical species, filters by `--sites`
3. Sorts by fractional (x, y) coordinates for reproducible `{species}_top_{i}` labeling
4. Only species explicitly listed in `--sites` are included — H atoms at the surface are skipped unless you pass `--sites H`

#### Built-in Adsorbate Geometries

- **O**: Single atom at height above site
- **OH**: O at height, H at +0.98 Å straight up from O
- **OOH**: O1 at height, O2 at 1.21 Å tilted 60° from vertical, H at 0.98 Å from O2 along O1→O2

#### Selective Dynamics

Clean slab FixAtoms constraints are preserved. Adsorbate atoms are always free (not included in FixAtoms).

#### Design Choices

- **Ontop sites only** (monodentate). Bridge/hollow are future extensions — for ionic/oxide surfaces, bridge sites between atoms of opposite charge don't make physical sense for single-atom bonding.
- **No asetools dependency**: standalone implementation using ASE + numpy. The asetools `SurfaceAnalyzer` was designed for metallic fcc hollow-site chemistry and doesn't apply to general ontop placement.

## Pitfalls

These are MLIP-specific gotchas observed when relaxing the structures produced by this skill. They are not bugs in the placement script — they are reminders to validate the placement *after the downstream relaxation*. Concrete cutoffs (separation distance, height) depend on the surface, the MLIP, and the task head; the only reliable diagnostic is post-relaxation inspection.

- **Reactive co-adsorbates can spontaneously react during relaxation.** When placing two reactive species on the same slab (e.g. CO\* and O\* in a lateral-interaction study), the relaxation may form a new bond, dissociate one species, or migrate one onto the other's site. Diagnostic: open the relaxed `opt_final.vasp` and verify the species count matches what you placed. If a CO\* + O\* setup produced a CO₂\*, the placement was too close. There is no universal "safe" separation — concrete cutoffs depend on the surface, the MLIP backend, and the task head. Iterate by stepping the separation outward (e.g. 1NN → 2NN → 3NN) until the relaxed structure preserves both species. Document the spacing used.
- **Bent / multidentate species are height-sensitive.** Bidentate adsorbates (e.g. bent CO₂ binding through both O atoms, peroxide species, side-bound molecules) are stable only at low placement heights. The auto-height from covalent radii (single-bond geometry) often starts them too high and they desorb to the gas phase during relaxation. Diagnostic: after relaxation, check that the binding atoms have not risen by more than ~1 Å above the surface and that all intended surface bonds are present. If they have desorbed, override `--height` with a smaller value and rerun the relaxation. As before, the right number depends on the system; iterate.
- **Multi-atom adsorbates from ASE programmatically.** If you build adsorbates ad-hoc in Python instead of using a file, use `ase.build.molecule('CO')` to get an `Atoms` object — passing the string `'CO'` to `ase.build.add_adsorbate` raises `KeyError`. The `--adsorbate-file` path in this script avoids the issue by reading a structure file directly.

In short: the placement script puts atoms where you tell it; the MLIP decides whether that placement was sensible. Inspect the relaxed structure before trusting any energy.

## Agent Workflow Integration

This skill is the **fifth step** in a catalysis workflow:

1. **Fetch bulk** → bulk POSCAR
2. **Optimize bulk** → relax with MLIP (`optimize run --structure bulk.vasp --mlip uma-s-1p2 --uma-task omat`)
3. **Generate surfaces** → slab POSCARs with constraints
4. **Optimize clean slabs** → relaxed clean slabs
5. **Place adsorbates** (this skill) → slab+adsorbate input structures
6. **Relax slab+adsorbate** → MLIP optimization with `oc20` task head
7. **Compute adsorption energies** → rank sites

### Agent Decision Logic

- **Height**: auto-computed from covalent radii by default. Override with `--height` if needed (e.g., for large adsorbates, or if atoms overlap).
- **Depth**: 2.0 Å captures the corrugated top layer of most oxide surfaces. Decrease to 1.0 Å for flat metals, increase for heavily corrugated surfaces.
- **Site species**: Decide based on the chemistry. For oxides, typically adsorb on the metal cation (e.g. `--sites Co Ti`). Include anion sites (`--sites Co O`) if studying oxygen-on-oxygen binding, vacancy filling, or if the user requests it. For metals, use `--sites Pt Cu` etc. Skip H unless explicitly needed.
- **Binding atom**: For file-based adsorbates, ask the user which atom binds if not obvious. If not specified, defaults to atom index 0 (first atom in the file).
- **Relaxation**: Use `optimize run --mlip uma-s-1p2 --uma-task oc22` for oxide surfaces (CoOOH, Co3O4, TiO2, etc.) or `--uma-task oc20` for metal surfaces (Pt, Cu, etc.). Use `--uma-task oc25` for solid-liquid interface simulations with explicit solvent. Do NOT use `omat` for adsorbate-on-surface systems.
