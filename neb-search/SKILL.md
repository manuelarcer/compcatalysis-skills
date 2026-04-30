---
name: neb-search
description: >
  Find transition states and reaction barriers via climbing-image NEB,
  with pre- and post-run validators that catch the most common setup
  failures (atom-identity swaps, IS/FS geometric mismatches, excessive
  diffusion, inter-image jumps, non-monotonic bond evolution).
  Triggers: "NEB", "transition state", "find barrier", "reaction barrier",
  "climbing image", "minimum energy path", "MEP", "TS search", "kinetic
  barrier", "compute activation energy".
  Use when: user has optimized initial- and final-state slab structures
  and wants the reaction barrier between them. Wraps `mlip_platform`'s
  CustomNEB and adds defensive validation around it.
  NOT for: dimer-method TS refinement, IRC, vibrational analysis, or
  TS-only single-point calculations.
metadata: {"requires": {"packages": ["mlip-platform", "ase", "scipy", "numpy"], "python": "3.11+"}}
---

# NEB Search (Transition States & Reaction Barriers)

NEB (Nudged Elastic Band) with climbing image, plus the validation utilities
that turn most "barrierless 0.0 eV" and "phantom 4 Å diffusion" surprises
into clear pre-flight errors.

## Prerequisites

- `mlip-platform` installed for the `run` subcommand. The validate
  subcommands have no MLIP dependency.
- `ase`, `numpy`, `scipy` (scipy is used by the >2-atom crossing check).
- **Optimized** IS and FS endpoints. Endpoint quality is the single largest
  source of error in NEB — relax both with `--fmax 0.01` if you can.

## Scripts

### neb_workflow.py — three subcommands

**Location:** `scripts/neb_workflow.py`
**Helper:** `scripts/utils/neb_validation.py` (ports the validators developed
in prior MLIP catalysis runs)

#### 1. `validate-setup` — pre-flight check on IS/FS

```bash
python neb_workflow.py validate-setup \
    --initial IS_opt.vasp --final FS_opt.vasp \
    --adsorbate-indices 36 37
```

Checks:
- atom count + chemical-symbol match across IS/FS;
- per-atom displacement magnitudes (flags > `--max-disp` Å, default 3 Å — i.e.
  excessive diffusion that's not a direct reaction);
- atom-identity / crossing detection (Hungarian assignment for ≥3 atoms);
- displacement-direction analysis: parallel = co-diffusion (likely wrong),
  antiparallel = dissociation (correct);
- displacement asymmetry (suggests IS not centered between FS sites);
- substrate relaxation (large substrate motion = IS and FS surfaces don't match).

Exit code 0 = pass, 2 = issues found.

#### 2. `validate-path` — post-flight check on `A2B.traj`

```bash
python neb_workflow.py validate-path \
    --trajectory neb_v1/A2B.traj \
    --adsorbate-indices 36 37 \
    --monotonic-pairs 36-37
```

Checks:
- per-image inter-image jump (default `--max-jump 1.0` Å);
- adsorbate crossings along the path;
- monotonic evolution of the bond-distance pairs you flag (e.g. for
  dissociation, the breaking-bond distance should be monotonic).

#### 3. `run` — validate → NEB → validate

```bash
python neb_workflow.py run \
    --initial IS_opt.vasp --final FS_opt.vasp \
    --adsorbate-indices 36 37 \
    --num-images 7 \
    --mlip uma-s-1p1 --uma-task oc20 --fmax 0.05 \
    --max-steps 600 \
    --barrier-estimate 0.7 \
    --monotonic-pairs 36-37 \
    --output-dir neb_v1/
```

Runs `validate-setup` first; aborts if it fails (override with `--force`).
Calls `mlip_platform.core.neb.CustomNEB` with IDPP interpolation and climb
on (use `--no-climb` to disable). Runs `validate-path` after.

Outputs (in `--output-dir`): `A2B.traj`, `A2B_full.traj`, `neb.log`,
`neb_convergence.csv`, `neb_convergence.png`, `neb_data.csv`,
`neb_energy.png`, `neb_parameters.txt`, `00/POSCAR`...`N/POSCAR`.

### Key arguments for `run`

| Argument | Description | Default |
|---|---|---|
| `--initial` / `--final` | Optimized IS / FS structure files | required |
| `--adsorbate-indices` | Indices of moving atoms (used for validation) | required |
| `--num-images` | Total images including endpoints | `7` |
| `--mlip` | MLIP model | `uma-s-1p1` |
| `--uma-task` | UMA task head | `oc20` |
| `--fmax` | Force convergence on band | `0.05` |
| `--max-steps` | NEB steps cap | `600` |
| `--no-climb` | Disable climbing image | climb on by default |
| `--barrier-estimate` | Rough eV estimate; warns if `--num-images` is too coarse | — |
| `--monotonic-pairs` | e.g. `36-37` for breaking-bond monotonicity | — |
| `--skip-pre-validate` / `--skip-post-validate` / `--force` | Escape hatches | off |

## Conventions and Heuristics

These encode lessons from prior MLIP NEB runs.

### Image count vs. barrier

Approximately one image per 0.1 eV of barrier; minimum 5. The CO\* + O\* →
CO₂\* run (0.84 eV) used 7 images and converged cleanly. The lower-barrier
O₂ dissociation (0.63 eV) needed 11 images — coarser grids skipped the TS.
Pass `--barrier-estimate 0.7` and the script will flag if `--num-images` is
coarse for that estimate.

### Endpoint preparation, especially for dissociation

For dissociation reactions, the IS must be **geometrically centered between
the FS sites**, not independently optimized. A canonical mistake: pick the
lowest-energy O₂\* IS (bridge site) and the lowest-energy 2O\* FS (two fcc
sites that aren't symmetric about the IS), and the NEB shows 3 Å of
spurious diffusion plus a barrierless 0.0 eV "barrier". Always:

1. Pick the FS sites first (e.g. two adjacent fcc sites).
2. Build the IS at the midpoint of those FS sites (not at a separate
   minimum), oriented along the FS-site axis.
3. Optimize that IS, but only with mild constraints if the structure
   wants to migrate elsewhere.

### Atom-order preservation across IS/FS

`add_adsorbate` and `extend` calls must produce the same atom ordering in
IS and FS. If atom 37 is "CO's O" in IS but "bridge O" in FS, the NEB will
show CO₂ rotating non-physically through intermediate images. The
crossing-detection in `validate-setup` catches this (and reports the
optimal mapping) but it's much cheaper to fix the construction.

### Climb is on by default

The reported barrier from a non-climbing NEB is a lower bound. Always run
with climb unless you're doing a very fast scoping pass.

### Task-head consistency

NEB endpoints, the band itself, and any reference energies must all be
computed with the same `--uma-task`. `oc20` for metal surfaces, `oc22` for
oxides, `oc25` for solid–liquid interfaces. (See relaxation-orchestration's
"Task-head selection" table.)

## Pitfalls

- **"Barrierless 0.0 eV" almost always means bad endpoints.** Run
  `validate-setup` before believing it. Look for displacement angle <60°
  between adsorbate atoms (= co-diffusion, not dissociation).
- **Big inter-image jumps post-run** indicate the band is too coarse or
  the IS/FS pair was geometrically incompatible. Add images, or rebuild
  endpoints.
- **Reactive co-adsorbates can rearrange during endpoint optimization.**
  CO\* + O\* at 1NN spontaneously formed CO₂\* during plain relaxation on
  Pt(111)/oc20, so the "IS" you fed to NEB may not be what you intended.
  Inspect `IS_opt.vasp` carefully or place co-adsorbates ≥2NN apart (see
  the adsorbate-placement skill).
- **scipy is required for ≥3-atom crossing detection.** Without it the
  pairwise check still runs (no scipy dependency for 2-atom case).

## Agent Workflow Integration

This skill is a **branch off step 6** in the catalysis pipeline: instead of
computing equilibrium adsorption energies, you compute reaction barriers.

1. Fetch bulk
2. Optimize bulk
3. Generate surfaces
4. Optimize clean slabs
5. Place adsorbates → IS and FS configurations
6. **Optimize IS and FS endpoints with tight fmax (0.01–0.02)** via
   relaxation-orchestration
7. **Run NEB with climb** (this skill, `run` subcommand) to get the barrier
8. Optionally: sample additional barriers with the same skill,
   then aggregate into a kinetics file consumed by KMC or a microkinetic model.

### Decision logic

- **Endpoint optimization first.** If `IS_opt.vasp` and `FS_opt.vasp` don't
  exist, run relaxation-orchestration on both before invoking this skill.
- **Number of images.** Default 7. Pass `--barrier-estimate` and let the
  warning guide you, or step up to 11+ for sub-eV barriers.
- **Adsorbate indices.** These are the indices of the atoms that *move*
  during the reaction. The substrate is excluded (the validator complains
  if the substrate moves > 0.3 Å).
- **Monotonic pairs.** For dissociation A–B → A + B, pass
  `--monotonic-pairs <idx_A>-<idx_B>`. The breaking-bond distance should
  grow monotonically along the path.
