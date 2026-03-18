---
name: materials-project
description: >
  Fetch bulk crystal structures from Materials Project for computational catalysis workflows.
  Triggers: "fetch bulk", "get structure from Materials Project", "download crystal structure",
  "MP structure", "bulk POSCAR", "get me [element] bulk", "materials project [formula]",
  "polymorph", "most stable structure of [formula]".
  Use when: user wants a bulk crystal structure by formula (Pt, TiO2) or Materials Project ID (mp-126),
  needs to compare polymorphs, or is starting a catalysis workflow that requires a bulk input.
  NOT for: surface generation, adsorbate placement, or MLIP relaxation (those are separate skills).
metadata: {"requires": {"packages": ["mp-api", "pymatgen"], "env": ["MP_API_KEY"], "python": "3.11+"}}
---

# Materials Project Bulk Structure Retrieval

Fetch and prepare bulk crystal structures from the Materials Project database for use in computational catalysis workflows.

## Prerequisites

- Python package: `mp-api` (`pip install mp-api`)
- Environment variable: `MP_API_KEY` (get from https://next-gen.materialsproject.org/api)

## Script

**Location:** `scripts/fetch_bulk.py`

### Usage

```bash
# Fetch by formula (returns most stable polymorph)
python fetch_bulk.py --formula Pt --output Pt_bulk.vasp

# Fetch by Material ID
python fetch_bulk.py --mp-id mp-126 --output Pt_bulk.vasp

# Fetch conventional cell (default is primitive)
python fetch_bulk.py --formula Pt --conventional --output Pt_bulk_conv.vasp

# Search with filters
python fetch_bulk.py --formula Pt --max-ehull 0.025 --output Pt_bulk.vasp

# Output as CIF
python fetch_bulk.py --formula Cu --output Cu_bulk.cif --format cif

# Fetch multiple elements (e.g., for alloy references)
python fetch_bulk.py --formula Pt Cu Ni --output-dir bulk_structures/

# List available polymorphs without downloading
python fetch_bulk.py --formula TiO2 --list-polymorphs
```

### Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--formula` | Chemical formula(s) to search | — |
| `--mp-id` | Materials Project ID (e.g., `mp-126`) | — |
| `--output` | Output filename | `{formula}_bulk.vasp` |
| `--output-dir` | Output directory for multiple structures | `.` |
| `--format` | Output format: `vasp`, `cif`, `xyz`, `json` | `vasp` |
| `--conventional` | Use conventional standard cell | `False` (primitive) |
| `--max-ehull` | Max energy above hull in eV/atom | `0.025` |
| `--list-polymorphs` | List polymorphs and exit | `False` |
| `--api-key` | MP API key (overrides env var) | `$MP_API_KEY` |

### Output

- Structure file(s) in requested format
- Summary printed to stdout: formula, space group, lattice parameters, number of atoms, energy above hull
- When `--list-polymorphs`: table of all matching entries with material IDs, space groups, and energies

### Behavior

1. **Single formula**: fetches the thermodynamically most stable structure (lowest energy above hull ≤ `--max-ehull`)
2. **Material ID**: fetches that exact entry
3. **Multiple formulas**: fetches the most stable structure for each, saves to `--output-dir`
4. **Conventional cell**: applies `SpacegroupAnalyzer.get_conventional_standard_structure()`
5. **Polymorph listing**: queries all entries matching the formula and prints a summary table sorted by stability

### Error Handling

- Missing `MP_API_KEY` → clear error message with setup instructions
- No results found → reports "no structures found" with suggestions (check formula, try relaxed search)
- Multiple polymorphs within threshold → picks lowest energy, warns about alternatives

## Agent Workflow Integration

This skill is typically the **first step** in a catalysis workflow:

1. **Fetch bulk** (this skill) → produces bulk POSCAR
2. **Optimize bulk** → relax with MLIP (`optimize run --structure bulk.vasp --mlip uma-s-1p1 --uma-task omat`)
3. **Generate surface** → pymatgen SlabGenerator on relaxed bulk
4. **Place adsorbates** → enumerate sites, add species
5. **Relax slab + adsorbate** → MLIP optimization with constrained bottom layers
6. **Compute properties** → adsorption energies, barriers

### Agent Decision Logic

- **Which cell?** Use `--conventional` when you need to generate surfaces (SlabGenerator works better with conventional cells for some structures). Use primitive for bulk property calculations.
- **Which task head for optimization?** After fetching, bulk optimization should use `--uma-task omat` (materials task head).
- **Alloy bulks**: For binary/ternary alloys, search by full formula (e.g., `PtNi`) or build manually with pymatgen.
