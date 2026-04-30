"""Unit tests for relaxation-orchestration/scripts/batch_relax.py.

The MLIP itself is not exercised — `relax_one` is monkeypatched. The tests
cover the orchestration logic: structure discovery, resume detection,
summary CSV writing, and CLI argument plumbing.
"""

import csv
import json
import sys
import types
from pathlib import Path

import pytest
from ase import Atoms
from ase.io import write

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

def _seed_converged_dir(d: Path, *, last_fmax: float = 0.01, fmax_thr: float = 0.03):
    d.mkdir(parents=True, exist_ok=True)
    (d / "input.vasp").write_text("x")
    (d / "opt_final.vasp").write_text("x")
    (d / "opt_convergence.csv").write_text(
        f"step,energy(eV),fmax(eV/Å)\n1,0.0,0.5\n10,-1.0,{last_fmax}\n"
    )
    (d / batch_relax.RUN_META_FILENAME).write_text(json.dumps({
        "mlip": "uma-s-1p1", "uma_task": "oc20",
        "fmax": fmax_thr, "optimizer": "fire",
    }))
    return d / "input.vasp"


def test_is_already_converged_true(tmp_path):
    s = _seed_converged_dir(tmp_path / "ok", last_fmax=0.01, fmax_thr=0.03)
    assert batch_relax.is_already_converged(s, fmax=0.03) is True


def test_is_already_converged_false_when_above_threshold(tmp_path):
    s = _seed_converged_dir(tmp_path / "bad", last_fmax=0.10, fmax_thr=0.03)
    assert batch_relax.is_already_converged(s, fmax=0.03) is False


def test_is_already_converged_false_when_files_missing(tmp_path):
    d = tmp_path / "partial"
    d.mkdir()
    (d / "input.vasp").write_text("x")
    (d / "opt_final.vasp").write_text("x")  # missing csv
    assert batch_relax.is_already_converged(d / "input.vasp", fmax=0.03) is False


def test_is_already_converged_uses_recorded_fmax_when_none_passed(tmp_path):
    s = _seed_converged_dir(tmp_path / "auto", last_fmax=0.005, fmax_thr=0.03)
    # fmax=None → fall back to the recorded threshold (0.03 ≥ 0.005)
    assert batch_relax.is_already_converged(s) is True


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
    """Inject a lightweight fake mlip_platform so `relax_one` runs without the
    real MLIP. ASE is left as the real package so structure classification
    (vacuum detection) works on real .vasp inputs."""
    fake_utils = types.ModuleType("mlip_platform.cli.utils")
    fake_utils.setup_calculator = lambda atoms, mlip, uma_task: atoms

    def fake_run_optimization(atoms, optimizer, fmax, max_steps, output_dir,
                               model_name, verbose, trajectory="opt.traj",
                               logfile="opt.log"):
        out = Path(output_dir)
        (out / "opt.traj").write_text("x")
        (out / "opt_final.vasp").write_text("x")
        (out / "opt_convergence.csv").write_text(
            "step,energy(eV),fmax(eV/Å)\n1,0.0,0.001\n"
        )
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

    monkeypatch.setitem(sys.modules, "mlip_platform", fake_root)
    monkeypatch.setitem(sys.modules, "mlip_platform.cli", fake_cli)
    monkeypatch.setitem(sys.modules, "mlip_platform.cli.utils", fake_utils)
    monkeypatch.setitem(sys.modules, "mlip_platform.core", fake_core)
    monkeypatch.setitem(sys.modules, "mlip_platform.core.optimize", fake_optimize)


def test_main_dry_run_lists_pending(tmp_path, monkeypatch, capsys):
    d = tmp_path / "OH" / "Pt_top_0"
    d.mkdir(parents=True)
    write(str(d / "input.vasp"), _bulk_atoms(), format="vasp")

    monkeypatch.setattr(sys, "argv",
                        ["batch_relax.py", "--tree", str(tmp_path),
                         "--uma-task", "oc20", "--dry-run"])
    batch_relax.main()
    out = capsys.readouterr().out
    assert "Structures total: 1" in out
    assert "input.vasp" in out


def test_main_resume_skips_converged(tmp_path, monkeypatch, capsys):
    pending = tmp_path / "OH" / "Pt_top_0"
    pending.mkdir(parents=True)
    write(str(pending / "input.vasp"), _bulk_atoms(), format="vasp")
    converged_dir = tmp_path / "OH" / "Pt_top_1"
    converged_dir.mkdir(parents=True)
    write(str(converged_dir / "input.vasp"), _bulk_atoms(), format="vasp")
    _seed_converged_dir(converged_dir, last_fmax=0.01, fmax_thr=0.03)

    _build_fake_mlip_modules(monkeypatch)
    monkeypatch.setattr(sys, "argv",
                        ["batch_relax.py", "--tree", str(tmp_path),
                         "--uma-task", "oc20",
                         "--resume", "--summary", "out.csv"])
    batch_relax.main()
    out = capsys.readouterr().out
    assert "skipped (resume): 1" in out
    summary = tmp_path / "out.csv"
    assert summary.exists()
    rows = list(csv.DictReader(open(summary)))
    assert len(rows) == 1
    assert rows[0]["converged"] == "True"


# ---- structure classification + auto task inference ------------------------

def _bulk_atoms():
    return Atoms("Pt4", positions=[[0,0,0],[2,2,0],[2,0,2],[0,2,2]],
                  cell=[4, 4, 4], pbc=True)


def _slab_atoms():
    """Bottom-anchored 4-atom Pt slab in a 4×4×30 box (vacuum on c)."""
    return Atoms("Pt4",
                  positions=[[0,0,0],[2,2,0],[2,0,2.3],[0,2,2.3]],
                  cell=[4, 4, 30], pbc=True)


def _molecule_atoms():
    return Atoms("CO", positions=[[10, 10, 10], [10, 10, 11.13]],
                  cell=[20, 20, 20], pbc=True)


def test_classify_bulk_slab_molecule():
    assert batch_relax.classify_structure(_bulk_atoms()) == "bulk"
    assert batch_relax.classify_structure(_slab_atoms()) == "slab"
    assert batch_relax.classify_structure(_molecule_atoms()) == "molecule"


def test_infer_uma_task_bulks(tmp_path):
    p = tmp_path / "bulk.vasp"
    write(str(p), _bulk_atoms(), format="vasp")
    assert batch_relax.infer_uma_task([p]) == "omat"


def test_infer_uma_task_slabs_returns_none(tmp_path):
    p = tmp_path / "slab.vasp"
    write(str(p), _slab_atoms(), format="vasp")
    assert batch_relax.infer_uma_task([p]) is None


def test_infer_uma_task_molecules(tmp_path):
    p = tmp_path / "mol.vasp"
    write(str(p), _molecule_atoms(), format="vasp")
    assert batch_relax.infer_uma_task([p]) == "omol"


def test_infer_uma_task_mixed_returns_none(tmp_path):
    pb = tmp_path / "bulk.vasp"
    ps = tmp_path / "slab.vasp"
    write(str(pb), _bulk_atoms(), format="vasp")
    write(str(ps), _slab_atoms(), format="vasp")
    assert batch_relax.infer_uma_task([pb, ps]) is None


def test_main_auto_task_refuses_for_slab(tmp_path, monkeypatch, capsys):
    p = tmp_path / "input.vasp"
    write(str(p), _slab_atoms(), format="vasp")
    _build_fake_mlip_modules(monkeypatch)

    monkeypatch.setattr(sys, "argv",
                        ["batch_relax.py", "--structure", str(p),
                         "--dry-run"])  # default --uma-task=auto
    with pytest.raises(SystemExit) as excinfo:
        batch_relax.main()
    assert excinfo.value.code == 2
    err = capsys.readouterr().err
    assert "could not safely infer --uma-task" in err
    assert "oc20" in err


def test_main_auto_task_picks_omat_for_bulk(tmp_path, monkeypatch, capsys):
    p = tmp_path / "input.vasp"
    write(str(p), _bulk_atoms(), format="vasp")
    _build_fake_mlip_modules(monkeypatch)

    monkeypatch.setattr(sys, "argv",
                        ["batch_relax.py", "--structure", str(p), "--dry-run"])
    batch_relax.main()
    out = capsys.readouterr().out
    assert "Auto-inferred --uma-task = omat" in out


# ---- run metadata sidecar ---------------------------------------------------

def test_write_run_metadata(tmp_path):
    out_dir = tmp_path / "out"
    batch_relax.write_run_metadata(
        out_dir, mlip="uma-s-1p1", uma_task="oc20",
        fmax=0.03, optimizer="fire",
    )
    meta_path = out_dir / batch_relax.RUN_META_FILENAME
    assert meta_path.exists()
    data = json.loads(meta_path.read_text())
    assert data == {"mlip": "uma-s-1p1", "uma_task": "oc20",
                    "fmax": 0.03, "optimizer": "fire"}


def test_relax_one_writes_metadata(tmp_path, monkeypatch):
    p = tmp_path / "input.vasp"
    write(str(p), _bulk_atoms(), format="vasp")
    _build_fake_mlip_modules(monkeypatch)

    batch_relax.relax_one(p, mlip="uma-s-1p1", uma_task="oc20",
                          optimizer="fire", fmax=0.03,
                          max_steps=10, verbose=False)
    meta = json.loads((p.parent / batch_relax.RUN_META_FILENAME).read_text())
    assert meta["uma_task"] == "oc20"
    assert meta["mlip"] == "uma-s-1p1"
