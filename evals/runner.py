#!/usr/bin/env python3
"""Skill behavioral eval runner.

Dispatches a fresh `claude -p` subagent per scenario with the target skill's
SKILL.md injected via --append-system-prompt, runs it in an isolated temp
working directory, and grades produced artifacts deterministically.

Usage:
    python evals/runner.py --list
    python evals/runner.py --skill surface-generation
    python evals/runner.py --scenario evals/scenarios/surface-generation/Pt111_basic.yaml
    python evals/runner.py --dry-run                 # validate scenarios without calling the agent
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from grading import grade_scenario  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
SCENARIOS_DIR = REPO_ROOT / "evals" / "scenarios"
FIXTURES_DIR = REPO_ROOT / "evals" / "fixtures"
RESULTS_DIR = REPO_ROOT / "evals" / "results"


def load_scenario(path: Path) -> dict:
    with open(path) as f:
        sc = yaml.safe_load(f)
    sc["_path"] = str(path)
    sc.setdefault("name", path.stem)
    if "skill" not in sc:
        raise ValueError(f"{path}: 'skill' field is required")
    if "prompt" not in sc:
        raise ValueError(f"{path}: 'prompt' field is required")
    return sc


def discover_scenarios(skill: str | None = None) -> list[Path]:
    paths = sorted(SCENARIOS_DIR.rglob("*.yaml"))
    if skill:
        paths = [p for p in paths if p.parent.name == skill]
    return paths


def prepare_workdir(scenario: dict) -> Path:
    workdir = Path(tempfile.mkdtemp(prefix=f"eval_{scenario['name']}_"))
    for seed in scenario.get("setup", []):
        src = (REPO_ROOT / seed["from"]).resolve()
        dst = workdir / seed["to"]
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    return workdir


def build_system_prompt(scenario: dict, python_bin: str) -> str:
    skill_md = (REPO_ROOT / scenario["skill"] / "SKILL.md").read_text()
    return (
        f"You are being evaluated on the `{scenario['skill']}` skill. "
        f"The full SKILL.md content follows between <skill> tags.\n\n"
        f"<skill name=\"{scenario['skill']}\">\n{skill_md}\n</skill>\n\n"
        f"Execution context:\n"
        f"- Repo root: {REPO_ROOT}\n"
        f"- Python interpreter for running scripts: {python_bin}\n"
        f"- Your working directory is a fresh temp directory. Put all output "
        f"artifacts there (relative paths in your final answer).\n"
        f"- Use Bash to invoke the scripts referenced in the SKILL.md.\n"
    )


def check_requirements(scenario: dict) -> tuple[bool, str]:
    for var in scenario.get("requires_env", []):
        if not os.environ.get(var):
            return False, f"required env var not set: {var}"
    return True, ""


def run_agent(scenario: dict, workdir: Path, python_bin: str,
              model: str, budget_usd: float) -> dict:
    system_prompt = build_system_prompt(scenario, python_bin)
    cmd = [
        "claude", "-p", scenario["prompt"],
        "--append-system-prompt", system_prompt,
        "--add-dir", str(REPO_ROOT),
        "--dangerously-skip-permissions",
        "--output-format", "json",
        "--no-session-persistence",
        "--model", model,
        "--max-budget-usd", str(budget_usd),
        "--disable-slash-commands",
    ]
    env = {**os.environ, "COMPCAT_PYTHON": python_bin}
    t0 = time.time()
    proc = subprocess.run(
        cmd, cwd=str(workdir), env=env,
        capture_output=True, text=True, timeout=600,
    )
    elapsed = time.time() - t0

    info = {"return_code": proc.returncode, "duration_s": elapsed,
            "stderr_tail": proc.stderr[-500:] if proc.stderr else ""}
    try:
        info["agent_result"] = json.loads(proc.stdout)
    except json.JSONDecodeError:
        info["agent_result"] = None
        info["stdout_tail"] = proc.stdout[-500:]
    return info


def run_scenario(scenario: dict, *, python_bin: str, model: str,
                 budget_usd: float, dry_run: bool) -> dict:
    ok, reason = check_requirements(scenario)
    if not ok:
        return {"scenario": scenario["name"], "skill": scenario["skill"],
                "skipped": True, "reason": reason}

    workdir = prepare_workdir(scenario)
    result = {"scenario": scenario["name"], "skill": scenario["skill"],
              "workdir": str(workdir), "dry_run": dry_run}

    if dry_run:
        result["system_prompt_len"] = len(build_system_prompt(scenario, python_bin))
        result["setup_files"] = [str(p.relative_to(workdir))
                                  for p in workdir.rglob("*") if p.is_file()]
        result["grading_keys"] = [a.get("path") or a.get("path_glob")
                                   for a in scenario.get("expect", {}).get("artifacts", [])]
        return result

    result["agent"] = run_agent(scenario, workdir, python_bin, model, budget_usd)
    agent_result = result["agent"].get("agent_result") or {}
    agent_text = agent_result.get("result") or ""
    result["grading"] = grade_scenario(workdir, scenario.get("expect", {}), agent_text)
    result["pass"] = result["grading"]["pass"] and result["agent"]["return_code"] == 0
    return result


def write_result(result: dict) -> Path:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    out = RESULTS_DIR / f"{ts}_{result['scenario']}.json"
    with open(out, "w") as f:
        json.dump(result, f, indent=2, default=str)
    return out


def print_summary(results: list[dict]) -> None:
    print(f"\n{'scenario':<32} {'skill':<22} {'status':<12} {'notes'}")
    print("-" * 90)
    for r in results:
        if r.get("skipped"):
            status = "SKIP"
            notes = r.get("reason", "")
        elif r.get("dry_run"):
            status = "DRY-OK"
            notes = f"setup={r.get('setup_files', [])}"
        else:
            status = "PASS" if r["pass"] else "FAIL"
            notes = f"rc={r['agent']['return_code']} t={r['agent']['duration_s']:.1f}s"
        print(f"{r['scenario']:<32} {r['skill']:<22} {status:<12} {notes}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--skill", help="Only run scenarios for this skill")
    p.add_argument("--scenario", help="Run a specific scenario YAML path")
    p.add_argument("--list", action="store_true", help="List discovered scenarios and exit")
    p.add_argument("--dry-run", action="store_true",
                   help="Validate scenarios and stage workdirs without calling the agent")
    p.add_argument("--model", default="sonnet", help="Model (default: sonnet)")
    p.add_argument("--budget-usd", type=float, default=2.0, help="Per-scenario USD cap")
    p.add_argument("--python-bin", default=sys.executable,
                   help="Python interpreter for the agent to run scripts with")
    args = p.parse_args()

    if args.scenario:
        paths = [Path(args.scenario).resolve()]
    else:
        paths = discover_scenarios(args.skill)

    if args.list:
        for path in paths:
            print(path.relative_to(REPO_ROOT))
        return

    if not paths:
        print("No scenarios matched.", file=sys.stderr)
        sys.exit(1)

    results = []
    for path in paths:
        print(f"\n>>> {path.relative_to(REPO_ROOT)}")
        scenario = load_scenario(path)
        result = run_scenario(scenario, python_bin=args.python_bin,
                              model=args.model, budget_usd=args.budget_usd,
                              dry_run=args.dry_run)
        results.append(result)
        if not args.dry_run:
            out = write_result(result)
            print(f"  → {out.relative_to(REPO_ROOT)}")

    print_summary(results)
    if any((not r.get("skipped") and not r.get("dry_run") and not r.get("pass"))
           for r in results):
        sys.exit(1)


if __name__ == "__main__":
    main()
