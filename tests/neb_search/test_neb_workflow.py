"""Unit tests for neb-search/scripts/neb_workflow.py.

The MLIP path is exercised via a fake CustomNEB so the test stays hermetic.
The validate-setup and validate-path subcommands run against real ASE
structures.
"""

import sys
import types
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import Trajectory, write

import neb_workflow


def _build_pair(tmp_path):
    """Write IS/FS POSCARs that should pass validate_neb_setup."""
    base = [[0, 0, 0], [4, 0, 0], [2, 3.46, 0]]
    is_atoms = Atoms(
        symbols=["Pt", "Pt", "Pt", "O", "O"],
        positions=base + [[6, 6, 6], [7.4, 6, 6]],
        cell=[20, 20, 20], pbc=True,
    )
    fs_atoms = Atoms(
        symbols=["Pt", "Pt", "Pt", "O", "O"],
        positions=base + [[5, 6, 6], [8.4, 6, 6]],
        cell=[20, 20, 20], pbc=True,
    )
    is_path = tmp_path / "IS.vasp"
    fs_path = tmp_path / "FS.vasp"
    write(str(is_path), is_atoms, format="vasp")
    write(str(fs_path), fs_atoms, format="vasp")
    return is_path, fs_path


def test_parse_pair():
    assert neb_workflow.parse_pair("36-37") == (36, 37)


def test_validate_setup_cli_pass(tmp_path, monkeypatch, capsys):
    is_p, fs_p = _build_pair(tmp_path)
    monkeypatch.setattr(sys, "argv", [
        "neb_workflow.py", "validate-setup",
        "--initial", str(is_p), "--final", str(fs_p),
        "--adsorbate-indices", "3", "4",
        "--quiet",
    ])
    rc = neb_workflow.main()
    assert rc == 0


def test_validate_setup_cli_fails_on_swap(tmp_path, monkeypatch):
    base = [[0, 0, 0], [4, 0, 0]]
    is_atoms = Atoms(
        symbols=["Pt", "Pt", "O", "O"],
        positions=base + [[6, 6, 6], [9, 6, 6]],
        cell=[20, 20, 20], pbc=True,
    )
    fs_atoms = Atoms(
        symbols=["Pt", "Pt", "O", "O"],
        positions=base + [[9, 6, 6], [6, 6, 6]],
        cell=[20, 20, 20], pbc=True,
    )
    is_p = tmp_path / "IS.vasp"
    fs_p = tmp_path / "FS.vasp"
    write(str(is_p), is_atoms, format="vasp")
    write(str(fs_p), fs_atoms, format="vasp")

    monkeypatch.setattr(sys, "argv", [
        "neb_workflow.py", "validate-setup",
        "--initial", str(is_p), "--final", str(fs_p),
        "--adsorbate-indices", "2", "3",
        "--quiet",
    ])
    rc = neb_workflow.main()
    assert rc == 2  # validation fails


def test_validate_path_cli_pass(tmp_path, monkeypatch):
    base = [[0, 0, 0], [4, 0, 0]]
    images = []
    for d_offset in [0.0, 0.5, 1.0, 1.5, 2.0]:
        atoms = Atoms(
            symbols=["Pt", "Pt", "O", "O"],
            positions=base + [[6 - 0.3 * d_offset, 6, 6],
                              [7.4 + 0.5 * d_offset, 6, 6]],
            cell=[20, 20, 20], pbc=True,
        )
        atoms.calc = SinglePointCalculator(atoms, energy=-10 + 0.5 * d_offset)
        images.append(atoms)
    traj_path = tmp_path / "A2B.traj"
    with Trajectory(str(traj_path), "w") as t:
        for img in images:
            t.write(img)

    monkeypatch.setattr(sys, "argv", [
        "neb_workflow.py", "validate-path",
        "--trajectory", str(traj_path),
        "--adsorbate-indices", "2", "3",
        "--monotonic-pairs", "2-3",
        "--quiet",
    ])
    rc = neb_workflow.main()
    assert rc == 0


def _inject_fake_mlip(monkeypatch):
    """Inject a fake mlip_platform so `run` exercises the orchestration path
    without launching a real MLIP."""
    fake_neb_module = types.ModuleType("mlip_platform.core.neb")

    class FakeCustomNEB:
        def __init__(self, initial, final, num_images, mlip, uma_task,
                     fmax, output_dir):
            self.output_dir = Path(output_dir)
            self.initial = initial
            self.final = final
            self.num_images = num_images

        def interpolate_idpp(self):
            pass

        def run_neb(self, climb=True, max_steps=600):
            # Write A2B.traj with the same atom layout as the IS so
            # adsorbate-indices stay valid.
            base = [[0, 0, 0], [4, 0, 0], [2, 3.46, 0]]
            self.output_dir.mkdir(parents=True, exist_ok=True)
            traj_path = self.output_dir / "A2B.traj"
            with Trajectory(str(traj_path), "w") as t:
                for d in [0.0, 0.5, 1.0]:
                    a = Atoms(
                        symbols=["Pt", "Pt", "Pt", "O", "O"],
                        positions=base + [[6 - 0.2 * d, 6, 6],
                                          [7.4 + 0.4 * d, 6, 6]],
                        cell=[20, 20, 20], pbc=True,
                    )
                    a.calc = SinglePointCalculator(a, energy=-10.0 + 0.3 * d)
                    t.write(a)

    fake_neb_module.CustomNEB = FakeCustomNEB
    fake_core = types.ModuleType("mlip_platform.core")
    fake_core.neb = fake_neb_module
    fake_root = types.ModuleType("mlip_platform")
    fake_root.core = fake_core
    monkeypatch.setitem(sys.modules, "mlip_platform", fake_root)
    monkeypatch.setitem(sys.modules, "mlip_platform.core", fake_core)
    monkeypatch.setitem(sys.modules, "mlip_platform.core.neb", fake_neb_module)


def test_run_subcommand_end_to_end(tmp_path, monkeypatch, capsys):
    is_p, fs_p = _build_pair(tmp_path)
    out_dir = tmp_path / "neb_out"
    _inject_fake_mlip(monkeypatch)

    monkeypatch.setattr(sys, "argv", [
        "neb_workflow.py", "run",
        "--initial", str(is_p), "--final", str(fs_p),
        "--adsorbate-indices", "3", "4",
        "--num-images", "5",
        "--output-dir", str(out_dir),
    ])
    rc = neb_workflow.main()
    assert rc == 0
    assert (out_dir / "A2B.traj").exists()


def test_run_aborts_on_bad_endpoints(tmp_path, monkeypatch, capsys):
    """Pre-flight failure should block the NEB unless --force is passed."""
    base = [[0, 0, 0], [4, 0, 0]]
    is_atoms = Atoms(
        symbols=["Pt", "Pt", "O"],
        positions=base + [[6, 6, 6]],
        cell=[20, 20, 20], pbc=True,
    )
    fs_atoms = Atoms(
        symbols=["Pt", "Pt", "O"],
        positions=base + [[12, 6, 6]],  # 6 Å excessive displacement
        cell=[20, 20, 20], pbc=True,
    )
    is_p = tmp_path / "IS.vasp"
    fs_p = tmp_path / "FS.vasp"
    write(str(is_p), is_atoms, format="vasp")
    write(str(fs_p), fs_atoms, format="vasp")

    _inject_fake_mlip(monkeypatch)
    monkeypatch.setattr(sys, "argv", [
        "neb_workflow.py", "run",
        "--initial", str(is_p), "--final", str(fs_p),
        "--adsorbate-indices", "2",
        "--num-images", "5",
        "--output-dir", str(tmp_path / "out"),
        "--max-disp", "3.0",
    ])
    rc = neb_workflow.main()
    assert rc == 2
    # NEB should not have been invoked
    assert not (tmp_path / "out" / "A2B.traj").exists()


# ---- new generality fixes (#7, #8) -----------------------------------------

def test_derive_num_images_heuristic():
    assert neb_workflow.derive_num_images(0.4) == 5  # floor
    assert neb_workflow.derive_num_images(0.7) == 7
    assert neb_workflow.derive_num_images(1.5) == 15


def test_run_refuses_without_num_images_or_barrier(tmp_path, monkeypatch, capsys):
    is_p, fs_p = _build_pair(tmp_path)
    _inject_fake_mlip(monkeypatch)

    monkeypatch.setattr(sys, "argv", [
        "neb_workflow.py", "run",
        "--initial", str(is_p), "--final", str(fs_p),
        "--adsorbate-indices", "3", "4",
        "--output-dir", str(tmp_path / "out"),
    ])
    rc = neb_workflow.main()
    assert rc == 2
    err = capsys.readouterr().err
    assert "no silent default" in err
    assert "--num-images" in err and "--barrier-estimate" in err


def test_run_derives_num_images_from_barrier_estimate(tmp_path, monkeypatch, capsys):
    is_p, fs_p = _build_pair(tmp_path)
    _inject_fake_mlip(monkeypatch)

    monkeypatch.setattr(sys, "argv", [
        "neb_workflow.py", "run",
        "--initial", str(is_p), "--final", str(fs_p),
        "--adsorbate-indices", "3", "4",
        "--barrier-estimate", "0.7",
        "--output-dir", str(tmp_path / "out"),
    ])
    rc = neb_workflow.main()
    assert rc == 0
    assert "Derived --num-images = 7" in capsys.readouterr().out


def test_run_warns_when_explicit_images_too_coarse(tmp_path, monkeypatch, capsys):
    is_p, fs_p = _build_pair(tmp_path)
    _inject_fake_mlip(monkeypatch)

    monkeypatch.setattr(sys, "argv", [
        "neb_workflow.py", "run",
        "--initial", str(is_p), "--final", str(fs_p),
        "--adsorbate-indices", "3", "4",
        "--num-images", "5",
        "--barrier-estimate", "1.5",   # heuristic suggests 15
        "--output-dir", str(tmp_path / "out"),
    ])
    rc = neb_workflow.main()
    assert rc == 0
    assert "may be coarse" in capsys.readouterr().out


def test_run_diffusion_reaction_type_relaxes_max_disp(tmp_path, monkeypatch, capsys):
    """A 5 Å direct displacement should pass when --reaction-type=diffusion
    (max-disp 8) but fail under direct (max-disp 3)."""
    base = [[0, 0, 0], [4, 0, 0]]
    is_atoms = Atoms(symbols=["Pt", "Pt", "O"], positions=base + [[6, 6, 6]],
                      cell=[20, 20, 20], pbc=True)
    fs_atoms = Atoms(symbols=["Pt", "Pt", "O"], positions=base + [[11, 6, 6]],
                      cell=[20, 20, 20], pbc=True)  # 5 Å motion
    is_p = tmp_path / "IS.vasp"
    fs_p = tmp_path / "FS.vasp"
    write(str(is_p), is_atoms, format="vasp")
    write(str(fs_p), fs_atoms, format="vasp")

    _inject_fake_mlip(monkeypatch)

    # direct → fail
    monkeypatch.setattr(sys, "argv", [
        "neb_workflow.py", "run",
        "--initial", str(is_p), "--final", str(fs_p),
        "--adsorbate-indices", "2",
        "--num-images", "5",
        "--reaction-type", "direct",
        "--output-dir", str(tmp_path / "out_direct"),
    ])
    assert neb_workflow.main() == 2

    # diffusion → pass
    monkeypatch.setattr(sys, "argv", [
        "neb_workflow.py", "run",
        "--initial", str(is_p), "--final", str(fs_p),
        "--adsorbate-indices", "2",
        "--num-images", "5",
        "--reaction-type", "diffusion",
        "--output-dir", str(tmp_path / "out_diff"),
    ])
    assert neb_workflow.main() == 0


def test_import_mlip_platform_friendly_error(monkeypatch, capsys):
    import builtins
    real_import = builtins.__import__

    def blocked_import(name, *args, **kwargs):
        if name.startswith("mlip_platform"):
            raise ImportError("blocked for test")
        return real_import(name, *args, **kwargs)

    for k in list(sys.modules):
        if k.startswith("mlip_platform"):
            monkeypatch.delitem(sys.modules, k, raising=False)
    monkeypatch.setattr(builtins, "__import__", blocked_import)

    with pytest.raises(SystemExit) as excinfo:
        neb_workflow._import_mlip_platform()
    assert excinfo.value.code == 2
    err = capsys.readouterr().err
    assert "mlip-platform is not available" in err
