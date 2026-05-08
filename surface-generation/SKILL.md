---
name: surface-generation
description: >
  Generate surface slabs from bulk crystal structures using pymatgen SlabGenerator.
  Triggers: "generate surface", "create slab", "surface slab", "Miller index",
  "cut surface", "cleave surface", "make slab from bulk", "enumerate terminations",
  "find all surfaces", "surface facets", "surface energy", "termination audit".
  Use when: user has a bulk structure and wants to generate surface slabs for a specific
  Miller index or enumerate all symmetrically distinct surfaces up to a max index.
  NOT for: adsorbate placement, relaxation, or adsorption energy (those are separate skills).
metadata: {"requires": {"packages": ["pymatgen", "ase"], "python": "3.11+"}}
---

# Surface Slab Generation

Generate surface slabs from optimized bulk crystal structures using pymatgen's SlabGenerator.

## Prerequisites

- Python packages: `pymatgen`, `ase`
- Input: a bulk structure file (VASP POSCAR, CIF, XYZ, or JSON)
- **Use conventional cell** for best results (primitive cells can produce unexpected surface orientations)

## Scripts

### 1. generate_slabs.py — Slab generation

**Location:** `scripts/generate_slabs.py`

#### Usage

```bash
# Generate (111) surface slab from bulk Pt
python generate_slabs.py --input Pt_bulk.vasp --miller 1 1 1 --output Pt_111.vasp

# Generate with specific thickness and vacuum
python generate_slabs.py --input Pt_bulk.vasp --miller 1 1 1 --min-slab-size 12.0 --min-vacuum-size 15.0

# Enumerate all terminations for a given Miller index
python generate_slabs.py --input TiO2_bulk.vasp --miller 1 1 0 --all-terminations --output-dir TiO2_110_slabs/

# Generate all symmetrically distinct surfaces up to max Miller index
python generate_slabs.py --input Pt_bulk.vasp --max-index 2 --output-dir Pt_surfaces/

# List available surfaces without generating files
python generate_slabs.py --input Pt_bulk.vasp --max-index 2 --list-surfaces

# Force symmetric slabs (top = bottom surface)
python generate_slabs.py --input Pt_bulk.vasp --miller 1 1 1 --symmetrize

# Disable selective dynamics (all atoms free)
python generate_slabs.py --input Pt_bulk.vasp --miller 1 1 1 --fix-bottom 0

# Fix only 30% of bottom atoms
python generate_slabs.py --input Pt_bulk.vasp --miller 1 1 1 --fix-bottom 0.3

# Use finer ftol for complex oxides
python generate_slabs.py --input Co3O4_bulk.vasp --miller 1 1 0 --ftol 0.05 --all-terminations
```

#### Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--input` | Input bulk structure file | — (required) |
| `--miller` | Miller indices h k l | — |
| `--max-index` | Max Miller index for enumerating all surfaces | — |
| `--min-slab-size` | Minimum slab thickness (Å) | `10.0` |
| `--min-vacuum-size` | Minimum vacuum thickness (Å) | `15.0` |
| `--all-terminations` | Generate all unique terminations | `False` (most stable only) |
| `--symmetrize` | Pass `symmetrize=True` to pymatgen SlabGenerator (does not always succeed) | `False` |
| `--require-symmetric` | Drop terminations whose slab is not symmetric (top != bottom) | `False` |
| `--require-stoichiometric` | Drop terminations whose composition does not match bulk | `False` |
| `--center-slab` | Center slab in vacuum | `True` |
| `--fix-bottom` | Fraction of atoms to fix at bottom | `0.5` |
| `--output` | Output filename (single slab) | auto-generated |
| `--output-dir` | Output directory (multiple slabs) | `.` |
| `--format` | Output format: `vasp`, `cif`, `xyz`, `json` | `vasp` |
| `--list-surfaces` | List surfaces and exit (no files written) | `False` |
| `--ftol` | Tolerance for termination clustering (Å) | `0.1` |
| `--max-broken-bonds` | Max broken bonds when cleaving | `0` |
| `--bonds` | Bond pairs to preserve, e.g. `Ti-O:2.0` | — |

#### Termination Detection — Important

This script uses a **catalysis-aware** termination filter instead of pymatgen's default symmetry filter. The difference:

- **pymatgen default**: considers slab A equivalent to slab B if one can be transformed into the other (including top↔bottom flip). This is correct for crystallography but **wrong for catalysis** — a slab with O-rich top / Co-rich bottom is NOT equivalent to Co-rich top / O-rich bottom.
- **Our approach**: disables pymatgen's symmetry filter, then deduplicates by **top-surface structural signature** (composition + z-layering). This ensures all catalytically distinct terminations are found.

For complex oxides (spinels, perovskites, etc.), use `--ftol 0.05` to resolve closely-spaced sub-layers. Run the termination audit script first to check.

#### Selective Dynamics

By default, the bottom 50% of atoms (by z-coordinate) are fixed (`--fix-bottom 0.5`). This embeds `FixAtoms` constraints in the POSCAR via selective dynamics, which ASE respects during optimization. Set `--fix-bottom 0` to disable.

### 2. termination_audit.py — Verify termination completeness

**Location:** `scripts/termination_audit.py`

Run this **before** generating slabs for a new material to verify that the chosen `ftol` finds all possible terminations.

```bash
# Audit Co3O4(110) terminations
python termination_audit.py --input bulk.vasp --miller 1 1 0

# Custom ftol sweep range
python termination_audit.py --input bulk.vasp --miller 1 1 0 --ftol-range 0.01 0.5 20
```

Output:
- Oriented unit cell analysis: all atomic planes, gaps, species
- Flags gaps smaller than 0.3 Å that may cause terminations to be lost
- Sweeps ftol and shows how termination count varies
- Recommends optimal ftol if default misses terminations

### 3. surface_energy.py — Compare termination stability

**Location:** `scripts/surface_energy.py`

After optimizing all terminations, compute surface energies to rank them.

```bash
python surface_energy.py --bulk-traj bulk/opt.traj --slab-dirs clean/t0 clean/t1 clean/t2 clean/t3 --output clean/surfenergy.log
```

Formula: `E_surf = (E_slab - (N_slab / N_bulk) * E_bulk) / (2 * A)`

Output: `surfenergy.log` with a comparison table (eV/Å² and J/m²) and identification of the most stable termination.

### Output

- Structure file(s) in requested format, with selective dynamics
- Summary: Miller index, termination shift, surface area, atom count, fixed/free, symmetry
- Filenames: `{formula}_{hkl}_t{termination_index}.vasp`

### Behavior

1. **Single Miller index**: generates the first termination unless `--all-terminations`
2. **All terminations**: generates one file per unique termination (catalysis-aware deduplication)
3. **Max index scan**: generates all symmetrically distinct surfaces up to `--max-index`
4. **Bond preservation**: when `--bonds` is set, avoids cleaving through those bonds
5. **Symmetrization**: when `--symmetrize`, ensures top and bottom surfaces are equivalent

### Error Handling

- Missing input file → clear error with suggestion
- No slabs found (all break bonds) → reports and suggests relaxing constraints
- Invalid Miller index → reports valid format

## Agent Workflow Integration

This skill is typically the **third step** in a catalysis workflow:

1. **Fetch bulk** → produces bulk POSCAR
2. **Optimize bulk** → relax with MLIP (`optimize run --structure bulk.vasp --mlip uma-s-1p2 --uma-task omat`)
3. **Audit terminations** → `termination_audit.py` to check ftol
4. **Generate surfaces** (this skill) → produces slab POSCAR(s) with constraints
5. **Optimize each termination** → one folder per termination under `clean/`
6. **Compare surface energies** → `surface_energy.py` writes `surfenergy.log`
7. **Place adsorbates** → on the most stable (or all relevant) termination(s)
8. **Relax slab + adsorbate** → MLIP optimization
9. **Compute properties** → adsorption energies

### Agent Decision Logic

- **Which cell?** Always use conventional cell input for surface generation
- **Slab thickness**: 7–10 Å minimum for screening, 15+ Å for production
- **Vacuum**: 15 Å standard, 20 Å if computing work functions or dipole corrections
- **ftol**: Use 0.1 for simple metals, 0.05 for complex oxides. Always audit first.
- **Constraints**: Default fix bottom 50%. Adjust if user specifies layer count.
- **Common catalysis surfaces**: Pt(111), Cu(111), TiO2(110), CeO2(111), Co3O4(110)
- **Oxides**: Use `--bonds Ti-O:2.0` (or similar) to avoid cleaving through metal-oxygen bonds
- **Termination selection**: Optimize ALL terminations, compute surface energies, use `surfenergy.log` to decide
- **Symmetric slabs**: Preferred for adsorption energy calculations to avoid dipole artifacts
