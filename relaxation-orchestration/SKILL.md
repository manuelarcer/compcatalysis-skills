---
name: relaxation-orchestration
description: >
  Run MLIP-backed geometry optimization on bulk, slab, gas-phase, and
  slab+adsorbate structures, with consistent task-head conventions for
  catalysis workflows.
  Triggers: "relax structure", "optimize structure", "geometry optimization",
  "MLIP relaxation", "minimize forces", "relax slab", "relax adsorbate",
  "batch optimize", "optimize all placements".
  Use when: user has one or more structure files and needs to run MLIP
  geometry optimization. Wraps the `mlip_platform` `optimize run` CLI and
  adds batch + resume support across an adsorbate-placement output tree.
  NOT for: NEB / transition states (use neb-search), adsorption-energy
  post-processing (use adsorption-energy), or running DFT directly.
metadata: {"requires": {"packages": ["mlip-platform", "ase"], "python": "3.11+"}}
---

# Relaxation Orchestration

Thin wrapper around `mlip_platform`'s `optimize run` command for catalysis
workflows. Two modes:

1. **Single structure** — direct CLI invocation, equivalent to calling
   `optimize run` yourself.
2. **Batch over a tree** — scans a directory for input structures (e.g. the
   output of the adsorbate-placement skill) and relaxes them all with
   consistent settings, with `--resume` support.

## Prerequisites

- `mlip-platform` installed and a UMA / SevenNet / MACE checkpoint available.
  Check with `optimize run --help`.
- Input structures already prepared (POSCAR with selective dynamics for
  slabs; this skill does NOT add constraints — that's the job of
  surface-generation / adsorbate-placement).

## Scripts

### batch_relax.py — Single or batch geometry optimization

**Location:** `scripts/batch_relax.py`

#### Usage

```bash
# Relax a single structure — same effect as `optimize run --structure ...`
python batch_relax.py --structure slab.vasp --mlip uma-s-1p1 --uma-task oc20 --fmax 0.03

# Relax everything under a tree (e.g. adsorbate-placement output)
python batch_relax.py --tree placements/ --pattern 'input.vasp' \
    --mlip uma-s-1p1 --uma-task oc20 --fmax 0.03

# Resume an interrupted batch — skip already-converged structures
python batch_relax.py --tree placements/ --resume

# Dry run to inspect the work plan
python batch_relax.py --tree placements/ --dry-run
```

#### Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--structure` | Single input structure path | — |
| `--tree` | Directory tree to scan for inputs | — |
| `--pattern` | Filename pattern under `--tree` | `input.vasp` |
| `--mlip` | MLIP model (`uma-s-1p1`, `uma-s-1p2`, `mace`, `7net-mf-ompa`, `auto`) | `auto` |
| `--uma-task` | UMA task head: `omat`, `oc20`, `omol`, `odac`, or `auto` | `auto` |
| `--optimizer` | `bfgs`, `lbfgs`, `fire`, ... | `bfgs` |
| `--fmax` | Force convergence threshold (eV/Å) | `0.03` |
| `--max-steps` | Maximum optimization steps | `300` |
| `--resume` | Skip structures with an existing converged run | off |
| `--dry-run` | List pending structures and exit | off |
| `--summary` | Summary CSV filename | `batch_relax.csv` |

Exactly one of `--structure` / `--tree` is required.

#### Output (per structure, written next to the input)

- `opt.traj` — full ASE trajectory
- `opt_final.vasp` — relaxed structure
- `opt_convergence.csv`, `opt_convergence.png` — per-step energy & fmax
- `opt.log` — optimizer log
- `opt_params.txt` — settings + converged flag (used by `--resume`)

Aggregate: `batch_relax.csv` summarizing converged/duration per structure.

## Conventions

These conventions encode lessons from prior MLIP-based catalysis runs.
Override only when the user has a specific reason to.

### Task-head selection

The UMA task head determines which energy reference the model uses. **All
calculations that contribute to a single energy comparison must use the same
task head**, or the energies are not subtractable.

| System | Task head | Status | Notes |
|---|---|---|---|
| Bulk crystal optimization | `omat` | supported | General materials, no surfaces |
| Metal slab + adsorbate + gas-phase reference (Pt, Cu, Ni, …) | `oc20` | supported | OC20 is the canonical surface-chemistry head |
| Oxide slab + adsorbate + gas-phase reference (CoOOH, Co₃O₄, TiO₂) | `oc22` | requires `mlip_platform` to expose it | use `oc20` until then and verify against literature |
| Solid–liquid interface with explicit solvent | `oc25` | requires `mlip_platform` to expose it | only when explicit water/electrolyte is present |
| Isolated molecules in vacuum | `omol` | supported | Not used when the molecule is a reference for an adsorption energy on a slab — use the surface task head instead |
| OC20 datasets / DAC | `odac` | supported | rare in standard catalysis workflows |

For any one workflow (`clean slab` + `slab+ads` + `gas-phase ads ref`), pick
one task head and use it for all three.

### `--uma-task auto` (default)

The script default is `auto`. It only picks a task automatically when it
can be **certain**:

- All inputs are bulk crystals (no vacuum direction) → `omat`
- All inputs are isolated molecules in a vacuum box (≥2 vacuum directions
  and ≤30 atoms) → `omol`

Anything else (slabs, adsorbate-on-slab, mixed batches) → the script
**refuses with a non-zero exit and prints the decision table**. The agent
must pick explicitly. This prevents the silent-wrong-energy class of bug
where `oc20` defaults silently get used on bulk crystals (or vice versa).

### Run-metadata sidecar

Every successful relaxation writes `compcat_run.json` next to
`opt_final.vasp` capturing `{mlip, uma_task, fmax, optimizer}`. The
adsorption-energy skill reads these to verify all input trajectories share
the same MLIP and task head, and refuses if they don't.

### Convergence thresholds

- **Production**: `--fmax 0.03 eV/Å`. Tighter than the mlip_platform default
  of 0.05; this is what produced quantitative agreement with literature in
  prior work.
- **Endpoint preparation for NEB**: `--fmax 0.01` if you can afford it — NEB
  endpoint quality is the largest single source of error in transition-state
  searches.
- **Difficult / strained structures**: drop to `0.05` only if 0.03 fails to
  converge in `max_steps`. Document the looser threshold.

### Optimizer choice

This skill defaults to `bfgs` to match `mlip_platform`'s default — that's
the platform-validated choice and we don't override it here.

- Switch to `--optimizer fire` when `bfgs` fails to converge or oscillates
  on noisy MLIP forces (more common on large systems and softer surfaces).
- For very small molecules (≤4 atoms) `bfgs` is usually fastest.
- `lbfgs` can be a good middle ground.

If you don't know which to pick, start with the default and only switch
when the run fails to converge in `max_steps`.

## Pitfalls

- **Mixing task heads silently invalidates energies.** A clean slab relaxed
  with `omat` is not subtractable from a slab+ads relaxed with `oc20`. The
  energies look reasonable individually but the adsorption energy is wrong
  by hundreds of meV. Always check `opt_params.txt` records the task you
  expected.
- **Selective dynamics is preserved from the input.** This skill does not
  add or modify constraints. If your slab+ads structure is missing
  `Selective dynamics`, the bottom of the slab will move during relaxation,
  invalidating the surface-energy and adsorption-energy comparisons. Fix
  upstream in surface-generation / adsorbate-placement.
- **Re-running on an already-relaxed `opt_final.vasp` is fine** but
  meaningless — it will converge in 0–2 steps and the new files overwrite
  the old. Use `--resume` to avoid this in batch mode.
- **`omat` for adsorbate-on-surface gives wrong adsorption energies.** This
  was confirmed in the CC-presentation MLIP runs: `oc20` is required for
  surface chemistry on Pt, even though `omat` gives sensible-looking
  individual energies.

## Agent Workflow Integration

This skill is steps **2, 4, and 6** in the catalysis pipeline:

1. Fetch bulk → bulk POSCAR
2. **Optimize bulk** (this skill, `--uma-task omat`) → relaxed bulk
3. Generate surfaces → slab POSCARs
4. **Optimize clean slabs** (this skill, `--uma-task oc20`/`oc22`) → relaxed slabs
5. Place adsorbates → slab+adsorbate inputs
6. **Optimize slab+adsorbate** (this skill, batch over the placement tree)
7. Compute adsorption energies → use the adsorption-energy skill

For step 6, point `--tree` at the adsorbate-placement output directory; the
default `--pattern input.vasp` already matches that skill's output naming.

### Decision logic

- **Picking `--mlip`**: `auto` is the default — it delegates to
  `mlip_platform.detect_mlip()`, so this skill tracks platform improvements
  automatically. Pass an explicit name (`uma-s-1p1`, `uma-s-1p2`, `mace`,
  `7net-mf-ompa`) when you need a specific checkpoint for reproducibility
  or for cross-MLIP comparisons.
- **Picking `--uma-task`**: see the table above. When unsure, ask whether
  the system is a metal, oxide, or solid–liquid interface.
- **Batch vs single**: if the user has more than two structures, prefer
  `--tree` so the summary CSV captures the run.
