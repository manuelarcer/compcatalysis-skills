"""Unit tests for relaxation-orchestration/scripts/batch_relax.py.

The MLIP itself is not exercised — `relax_one` is monkeypatched. The tests
cover the orchestration logic: structure discovery, resume detection,
summary CSV writing, and CLI argument plumbing.
"""

import csv
import sys
import types
from pathlib import Path

import pytest

import batch_relax


# ---- find_structures --------------------------------------------------------

def test_find_structures_walks_tree(tmp_path):
    (tmp_path / "OH" / "Pt_top_0").mkdir(parents=True)
    (tmp_path / "OH" / "Pt_top_0" / "input.vasp").write_text("x")
    (tmp_path / "OH" / "Pt_top_1").mkdir(parents=True)
    (tmp_path / "OH" / "Pt_top_1" / "input.vasp").write_text("x")
    (tmp_path / "noise.txt").write_text("y")

    found = batch_relax.find_structures(tmp_path, "input.vasp")
    assert len(found) == 2
    assert all(p.name == "input.vasp" for p in found)


def test_find_structures_missing_tree_exits(tmp_path, capsys):
    with pytest.raises(SystemExit):
        batch_relax.find_structures(tmp_path / "nope", "input.vasp")
    assert "tree not found" in capsys.readouterr().err


# ---- is_already_converged ---------------------------------------------------

def _seed_converged_dir(d: Path, converged: bool = True):
    d.mkdir(parents=True, exist_ok=True)
    (d / "input.vasp").write_text("x")
    (d / "opt_final.vasp").write_text("x")
    (d / "opt_convergence.csv").write_text("step,energy\n1,0\n")
    (d / "opt_params.txt").write_text(
        f"Geometry Optimization Parameters\nConverged:         {converged}\n"
    )
    return d / "input.vasp"


def test_is_already_converged_true(tmp_path):
    s = _seed_converged_dir(tmp_path / "ok")
    assert batch_relax.is_already_converged(s) is True


def test_is_already_converged_false_when_not_converged(tmp_path):
    s = _seed_converged_dir(tmp_path / "bad", converged=False)
    assert batch_relax.is_already_converged(s) is False


def test_is_already_converged_false_when_files_missing(tmp_path):
    d = tmp_path / "partial"
    d.mkdir()
    (d / "input.vasp").write_text("x")
    (d / "opt_final.vasp").write_text("x")  # missing csv + params
    assert batch_relax.is_already_converged(d / "input.vasp") is False


# ---- write_summary ----------------------------------------------------------

def test_write_summary_roundtrip(tmp_path):
    rows = [
        {"structure": "a.vasp", "converged": True, "duration_s": 1.0},
        {"structure": "b.vasp", "converged": False, "duration_s": 2.0},
    ]
    out = tmp_path / "sum.csv"
    batch_relax.write_summary(rows, out)
    with open(out) as f:
        loaded = list(csv.DictReader(f))
    assert len(loaded) == 2
    assert loaded[0]["structure"] == "a.vasp"
    assert loaded[1]["converged"] == "False"


def test_write_summary_noop_on_empty_rows(tmp_path):
    out = tmp_path / "x.csv"
    batch_relax.write_summary([], out)
    assert not out.exists()


# ---- main entry: dry-run, resume, batch -------------------------------------

def _build_fake_mlip_modules(monkeypatch):
    """Inject lightweight fake mlip_platform + ase modules so `relax_one` runs
    without the real MLIP machinery. Each fake call writes the success markers
    that `is_already_converged` looks for, so tests can verify side effects."""
    fake_atoms = types.SimpleNamespace()

    fake_ase_io = types.ModuleType("ase.io")
    fake_ase_io.read = lambda path: fake_atoms
    fake_ase = types.ModuleType("ase")
    fake_ase.io = fake_ase_io

    fake_utils = types.ModuleType("mlip_platform.cli.utils")
    fake_utils.setup_calculator = lambda atoms, mlip, uma_task: atoms

    def fake_run_optimization(atoms, optimizer, fmax, max_steps, output_dir,
                               model_name, verbose, trajectory="opt.traj",
                               logfile="opt.log"):
        out = Path(output_dir)
        (out / "opt.traj").write_text("x")
        (out / "opt_final.vasp").write_text("x")
        (out / "opt_convergence.csv").write_text("step\n1\n")
        (out / "opt_params.txt").write_text("Converged:         True\n")
        return True

    fake_optimize = types.ModuleType("mlip_platform.core.optimize")
    fake_optimize.run_optimization = fake_run_optimization

    fake_core = types.ModuleType("mlip_platform.core")
    fake_core.optimize = fake_optimize

    fake_cli = types.ModuleType("mlip_platform.cli")
    fake_cli.utils = fake_utils

    fake_root = types.ModuleType("mlip_platform")
    fake_root.cli = fake_cli
    fake_root.core = fake_core

    monkeypatch.setitem(sys.modules, "ase", fake_ase)
    monkeypatch.setitem(sys.modules, "ase.io", fake_ase_io)
    monkeypatch.setitem(sys.modules, "mlip_platform", fake_root)
    monkeypatch.setitem(sys.modules, "mlip_platform.cli", fake_cli)
    monkeypatch.setitem(sys.modules, "mlip_platform.cli.utils", fake_utils)
    monkeypatch.setitem(sys.modules, "mlip_platform.core", fake_core)
    monkeypatch.setitem(sys.modules, "mlip_platform.core.optimize", fake_optimize)


def test_main_dry_run_lists_pending(tmp_path, monkeypatch, capsys):
    (tmp_path / "OH" / "Pt_top_0").mkdir(parents=True)
    (tmp_path / "OH" / "Pt_top_0" / "input.vasp").write_text("x")

    monkeypatch.setattr(sys, "argv",
                        ["batch_relax.py", "--tree", str(tmp_path), "--dry-run"])
    batch_relax.main()
    out = capsys.readouterr().out
    assert "Structures total: 1" in out
    assert "input.vasp" in out


def test_main_resume_skips_converged(tmp_path, monkeypatch, capsys):
    pending = tmp_path / "OH" / "Pt_top_0"
    pending.mkdir(parents=True)
    (pending / "input.vasp").write_text("x")
    _seed_converged_dir(tmp_path / "OH" / "Pt_top_1", converged=True)

    _build_fake_mlip_modules(monkeypatch)
    monkeypatch.setattr(sys, "argv",
                        ["batch_relax.py", "--tree", str(tmp_path),
                         "--resume", "--summary", "out.csv"])
    batch_relax.main()
    out = capsys.readouterr().out
    assert "skipped (resume): 1" in out
    summary = tmp_path / "out.csv"
    assert summary.exists()
    rows = list(csv.DictReader(open(summary)))
    assert len(rows) == 1
    assert rows[0]["converged"] == "True"
