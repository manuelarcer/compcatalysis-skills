# Skill behavioral evals

Layer 2 of the harness. Each scenario dispatches a fresh `claude -p` subagent
with only the target skill's `SKILL.md` injected via `--append-system-prompt`,
runs it in an isolated temp directory, and grades the produced artifacts
deterministically. This complements the `tests/` directory (Layer 1), which
unit-tests the underlying scripts directly.

## Prereqs

- `claude` CLI on `$PATH` (check: `claude --version`)
- Python env with `pyyaml`, `pymatgen`, `ase` (the repo's `catplat` conda env works)
- For the materials-project scenario: `MP_API_KEY` env var

## Running

```bash
# list discovered scenarios
python evals/runner.py --list

# validate scenarios and stage workdirs without calling the agent
python evals/runner.py --dry-run

# run everything (costs tokens!)
python evals/runner.py

# run one skill
python evals/runner.py --skill surface-generation

# run one scenario
python evals/runner.py --scenario evals/scenarios/adsorbate-placement/OH_on_Pt111.yaml
```

Results land in `evals/results/` as JSON (gitignored). The runner exits non-zero
if any non-skipped scenario fails.

## Scenario schema

```yaml
name: <string>                     # optional; defaults to file stem
skill: <skill-directory-name>      # required — SKILL.md loaded from <skill>/SKILL.md
requires_env: [VAR1, VAR2]         # optional — skip scenario if any are unset
setup:                             # optional — files copied into the workdir before run
  - from: <path-relative-to-repo>
    to: <path-relative-to-workdir>
prompt: |                          # required — natural-language task for the agent
  ...
expect:
  artifacts:
    - path: "file.vasp"            # OR path_glob: "pattern/*.vasp"
      min_count: 1                 # only for path_glob
      checks:
        - exists
        - pymatgen_parseable
        - ase_parseable
        - has_selective_dynamics
        - composition: "Pt"
        - contains_species: "O"
```

New checks live in `evals/grading.py` — register them in the `CHECKS` dict
and they become usable by name in any scenario.

## Verifying the harness discriminates

To prove the harness actually detects regressions, try this: deliberately
edit a `SKILL.md` to omit a critical flag (e.g., remove the mention of
`--fix-bottom` from `surface-generation/SKILL.md`), rerun that scenario, and
confirm the `has_selective_dynamics` check fails. Revert the edit afterward.
