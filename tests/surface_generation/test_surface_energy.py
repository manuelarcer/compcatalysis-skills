"""Unit tests for surface-generation/scripts/surface_energy.py."""

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write

import surface_energy as se


def _write_traj(path, atoms, energy):
    """Attach a single-point energy to atoms and write as ASE trajectory."""
    atoms.calc = SinglePointCalculator(atoms, energy=energy)
    write(str(path), atoms)


@pytest.fixture
def bulk_traj(tmp_path):
    atoms = Atoms("Pt4",
                  positions=[[0, 0, 0], [2, 2, 0], [2, 0, 2], [0, 2, 2]],
                  cell=[4, 4, 4], pbc=True)
    path = tmp_path / "bulk_opt.traj"
    _write_traj(path, atoms, energy=-24.0)
    return path


@pytest.fixture
def slab_traj(tmp_path):
    atoms = Atoms("Pt8",
                  positions=[[i % 2 * 2, (i // 2) % 2 * 2, i * 2.3] for i in range(8)],
                  cell=[4, 4, 30], pbc=True)
    d = tmp_path / "t0"
    d.mkdir()
    path = d / "opt.traj"
    _write_traj(path, atoms, energy=-47.0)
    return path


def test_compute_surface_energy_math():
    # E_surf = (E_slab - (N_slab/N_bulk)*E_bulk) / (2*A)
    e_ev_a2, e_j_m2 = se.compute_surface_energy(
        e_slab=-47.0, n_slab=8, e_bulk=-24.0, n_bulk=4, area=16.0,
    )
    expected_ev_a2 = (-47.0 - (8 / 4) * -24.0) / (2 * 16.0)
    assert e_ev_a2 == pytest.approx(expected_ev_a2)
    assert e_j_m2 == pytest.approx(expected_ev_a2 * 16.0217663)


def test_read_surface_area_orthorhombic():
    atoms = Atoms("Pt", positions=[[0, 0, 0]], cell=[3.0, 4.0, 20.0], pbc=True)
    assert se.read_surface_area(atoms) == pytest.approx(12.0)


def test_read_energy_from_traj_missing_exits(tmp_path, capsys):
    with pytest.raises(SystemExit):
        se.read_energy_from_traj(tmp_path / "nope.traj")
    assert "Trajectory not found" in capsys.readouterr().err


def test_read_energy_from_traj_roundtrip(bulk_traj):
    energy, atoms = se.read_energy_from_traj(bulk_traj)
    assert energy == pytest.approx(-24.0)
    assert len(atoms) == 4


def test_find_traj_in_dir_explicit_name(slab_traj):
    found = se.find_traj_in_dir(slab_traj.parent, traj_name="opt.traj")
    assert found == slab_traj


def test_find_traj_in_dir_glob_fallback(tmp_path):
    d = tmp_path / "t1"
    d.mkdir()
    p = d / "custom.traj"
    p.write_bytes(b"x")
    assert se.find_traj_in_dir(d, traj_name="nonexistent.traj") == p


def test_find_traj_in_dir_none_when_empty(tmp_path):
    d = tmp_path / "empty"
    d.mkdir()
    assert se.find_traj_in_dir(d) is None


def test_compute_slab_thickness_simple_layers():
    # Three layers at z = 0, 2, 4 (gap 2 Å > tol 0.5)
    atoms = Atoms("Pt6",
                  positions=[[0, 0, 0], [1, 0, 0],
                             [0, 0, 2], [1, 0, 2],
                             [0, 0, 4], [1, 0, 4]],
                  cell=[2, 2, 30], pbc=True)
    assert se.compute_slab_thickness(atoms) == pytest.approx(4.0)


def test_compute_slab_thickness_robust_to_dangling_atom():
    # Two clean layers at z=0,2; one protruding H at z=5
    # Layer-mean thickness should be 5 - 0 = 5 (top "layer" is single atom).
    # The point: a single dangling atom doesn't appear in the bottom layer
    # mean, so the bottom mean stays at 0 instead of being averaged up.
    atoms = Atoms("Pt4H",
                  positions=[[0, 0, 0], [1, 0, 0],
                             [0, 0, 2], [1, 0, 2],
                             [0, 0, 5]],
                  cell=[2, 2, 30], pbc=True)
    # bottom layer = {0, 0} -> mean 0; top layer = {5} -> mean 5
    assert se.compute_slab_thickness(atoms) == pytest.approx(5.0)


def test_compute_slab_thickness_buckled_layer_groups_together():
    # A buckled layer (z spread 0.3 Å) should group as one layer
    atoms = Atoms("Pt4",
                  positions=[[0, 0, 0.0], [1, 0, 0.3],
                             [0, 0, 2.0], [1, 0, 2.2]],
                  cell=[2, 2, 30], pbc=True)
    # bottom mean = 0.15, top mean = 2.1 -> 1.95
    assert se.compute_slab_thickness(atoms) == pytest.approx(1.95)


def test_compute_slab_thickness_single_layer_falls_back_to_span():
    atoms = Atoms("Pt2",
                  positions=[[0, 0, 0.0], [1, 0, 0.4]],
                  cell=[2, 2, 30], pbc=True)
    assert se.compute_slab_thickness(atoms) == pytest.approx(0.4)


def test_write_surfenergy_log_includes_thickness_column(tmp_path):
    bulk_info = {"path": "bulk.traj", "energy": -24.0, "n_atoms": 4}
    slab_results = [
        {"name": "t0", "e_slab": -47.0, "n_atoms": 8, "area": 16.0,
         "thickness": 6.5,
         "e_surf_ev_a2": 0.03, "e_surf_j_m2": 0.48, "path": "x", "traj": "y"},
    ]
    out = tmp_path / "surfenergy.log"
    se.write_surfenergy_log(out, bulk_info, slab_results)
    text = out.read_text()
    assert "Thickness (Å)" in text
    assert "6.5000" in text
    assert "z_top_layer" in text


def test_write_surfenergy_log_sorts_and_marks_stable(tmp_path, capsys):
    bulk_info = {"path": "bulk.traj", "energy": -24.0, "n_atoms": 4}
    slab_results = [
        {"name": "t1", "e_slab": -40.0, "n_atoms": 8, "area": 16.0,
         "e_surf_ev_a2": 0.5, "e_surf_j_m2": 8.01, "path": "x", "traj": "y"},
        {"name": "t0", "e_slab": -47.0, "n_atoms": 8, "area": 16.0,
         "e_surf_ev_a2": 0.03, "e_surf_j_m2": 0.48, "path": "x", "traj": "y"},
    ]
    out = tmp_path / "surfenergy.log"
    best = se.write_surfenergy_log(out, bulk_info, slab_results)

    text = out.read_text()
    # Lower-energy termination comes first, and is marked "<<"
    t0_line = next(ln for ln in text.splitlines() if ln.startswith("t0"))
    assert "<<" in t0_line
    assert best["name"] == "t0"
    assert "Most stable termination: t0" in text
    # Content also echoed to stdout
    assert "SURFACE ENERGY COMPARISON" in capsys.readouterr().out
