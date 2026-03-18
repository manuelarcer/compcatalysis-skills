# compcatalysis-skills

Agent skills for computational heterogeneous catalysis workflows.

These skills provide AI agent orchestration for end-to-end catalysis simulations, from structure retrieval to property calculation. They are designed to work with:

- **[mlip-platform](https://github.com/manuelarcer/mlip-platform)** — MLIP-based optimization, MD, and NEB
- **[asetools](https://github.com/manuelarcer/asetools)** — Structure analysis, adsorbate placement, workflow management
- **[ASE](https://wiki.fysik.dtu.dk/ase/)** — Atomic Simulation Environment
- **[pymatgen](https://pymatgen.org/)** — Structure manipulation, Materials Project API

## Skills

| Skill | Description | Status |
|-------|-------------|--------|
| `materials-project` | Fetch and prepare bulk structures from Materials Project | 🚧 In progress |
| `surface-gen` | Generate surface slabs from bulk structures | Planned |
| `adsorbate-placement` | Enumerate and place adsorbates on surfaces | Planned |
| `relaxation` | Run MLIP relaxations via mlip-platform | Planned |
| `adsorption-energy` | Compute and rank adsorption energies | Planned |
| `neb-workflow` | Set up and run NEB calculations | Planned |

## Setup

### Prerequisites

```bash
pip install mp-api pymatgen ase
```

### Materials Project API Key

Several skills require a Materials Project API key. Get one at: https://next-gen.materialsproject.org/api

Set it as an environment variable:
```bash
export MP_API_KEY="your_key_here"
```

Or add to your shell profile for persistence.

## Usage

Each skill directory contains a `SKILL.md` that an AI agent reads to understand how to use the skill. The scripts in `scripts/` can also be used standalone from the command line.

## License

MIT
