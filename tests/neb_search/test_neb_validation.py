"""Unit tests for neb-search/scripts/utils/neb_validation.py."""

import sys
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import Trajectory, write

# Add the skill's scripts/ path so `from utils.neb_validation import ...` works
SCRIPTS = Path(__file__).resolve().parents[2] / "neb-search" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from utils.neb_validation import validate_neb_setup, validate_neb_path  # noqa: E402


def _slab_with_adsorbates(positions, cell=20.0, symbols=None):
    """Build a 6-atom Pt slab plus arbitrary adsorbate atoms."""
    base_pt = [
        [0, 0, 0], [4, 0, 0], [2, 3.46, 0],
        [0, 0, 4], [4, 0, 4], [2, 3.46, 4],
    ]
    syms = ["Pt"] * len(base_pt) + (symbols or [])
    pos = base_pt + list(positions)
    return Atoms(symbols=syms, positions=pos,
                  cell=[cell, cell, cell], pbc=True)


# ---- validate_neb_setup -----------------------------------------------------

def test_setup_valid_dissociation():
    """Two O atoms moving apart from a midpoint = good dissociation IS/FS."""
    initial = _slab_with_adsorbates(
        [[6, 6, 6], [7.4, 6, 6]],  # O-O at 1.4 Å, midpoint (6.7, 6, 6)
        symbols=["O", "O"],
    )
    final = _slab_with_adsorbates(
        [[5, 6, 6], [8.4, 6, 6]],  # O atoms pulled apart, same midpoint
        symbols=["O", "O"],
    )
    result = validate_neb_setup(initial, final, adsorbate_indices=[6, 7],
                                  verbose=False)
    assert result["valid"]
    assert result["crossing"] is False


def test_setup_detects_atom_count_mismatch():
    initial = _slab_with_adsorbates([[6, 6, 6]], symbols=["O"])
    final = _slab_with_adsorbates([[6, 6, 6], [7, 6, 6]], symbols=["O", "O"])
    result = validate_neb_setup(initial, final, adsorbate_indices=[6],
                                  verbose=False)
    assert not result["valid"]
    assert any("count mismatch" in i for i in result["issues"])


def test_setup_detects_chemical_symbols_mismatch():
    initial = _slab_with_adsorbates([[6, 6, 6]], symbols=["O"])
    final = _slab_with_adsorbates([[6, 6, 6]], symbols=["N"])
    result = validate_neb_setup(initial, final, adsorbate_indices=[6],
                                  verbose=False)
    assert not result["valid"]
    assert any("symbols mismatch" in i for i in result["issues"])


def test_setup_detects_excessive_displacement():
    initial = _slab_with_adsorbates([[6, 6, 6]], symbols=["O"])
    final = _slab_with_adsorbates([[10, 6, 6]], symbols=["O"])  # 4 Å
    result = validate_neb_setup(initial, final, adsorbate_indices=[6],
                                  max_disp=3.0, verbose=False)
    assert not result["valid"]
    assert any("excessive diffusion" in i for i in result["issues"])


def test_setup_detects_co_diffusion():
    """Both atoms moving in the same direction = diffusion, not reaction."""
    initial = _slab_with_adsorbates(
        [[6, 6, 6], [7.4, 6, 6]], symbols=["O", "O"],
    )
    final = _slab_with_adsorbates(
        [[8, 6, 6], [9.4, 6, 6]],  # both shifted +2 Å in x
        symbols=["O", "O"],
    )
    result = validate_neb_setup(initial, final, adsorbate_indices=[6, 7],
                                  verbose=False)
    assert not result["valid"]
    assert any("similar directions" in i or "diffusion" in i.lower()
               for i in result["issues"])


def test_setup_detects_atom_swap():
    """Two atoms swap positions → must be flagged."""
    initial = _slab_with_adsorbates(
        [[6, 6, 6], [9, 6, 6]], symbols=["O", "O"],
    )
    final = _slab_with_adsorbates(
        [[9, 6, 6], [6, 6, 6]],  # swapped
        symbols=["O", "O"],
    )
    result = validate_neb_setup(initial, final, adsorbate_indices=[6, 7],
                                  verbose=False)
    assert result["crossing"] is True
    assert any("SWAP" in i for i in result["issues"])


# ---- validate_neb_path ------------------------------------------------------

def _path_traj(tmp_path, positions_per_image, energies, symbols=None):
    """Write a series of images to a trajectory file with single-point energies."""
    traj_path = tmp_path / "A2B.traj"
    with Trajectory(str(traj_path), "w") as t:
        for pos, e in zip(positions_per_image, energies):
            atoms = _slab_with_adsorbates(pos, symbols=symbols)
            atoms.calc = SinglePointCalculator(atoms, energy=e)
            t.write(atoms)
    return traj_path


def test_path_smooth_dissociation_valid(tmp_path):
    # Five images: O-O at 1.4 Å gradually growing to 4 Å, energies have a barrier
    images_pos = [
        [[6, 6, 6], [7.4, 6, 6]],
        [[5.7, 6, 6], [8.0, 6, 6]],
        [[5.4, 6, 6], [8.6, 6, 6]],
        [[5.1, 6, 6], [9.1, 6, 6]],
        [[4.9, 6, 6], [9.6, 6, 6]],
    ]
    energies = [-10.0, -9.5, -9.3, -9.5, -10.2]
    traj = _path_traj(tmp_path, images_pos, energies, symbols=["O", "O"])

    result = validate_neb_path(str(traj), adsorbate_indices=[6, 7],
                                expected_monotonic_pairs=[(6, 7)],
                                max_jump=1.0, verbose=False)
    assert result["valid"]
    assert result["max_jump"] < 1.0


def test_path_detects_inter_image_jump(tmp_path):
    images_pos = [
        [[6, 6, 6]],
        [[6.5, 6, 6]],
        [[10, 6, 6]],  # 3.5 Å jump from previous
        [[10.1, 6, 6]],
    ]
    traj = _path_traj(tmp_path, images_pos, [-10, -9, -8, -7], symbols=["O"])
    result = validate_neb_path(str(traj), adsorbate_indices=[6],
                                max_jump=1.0, verbose=False)
    assert not result["valid"]
    assert any("jumps" in i for i in result["issues"])


def test_path_detects_non_monotonic_distance(tmp_path):
    # d(6,7) supposed to grow monotonically, but image 2 dips back
    images_pos = [
        [[6, 6, 6], [7.4, 6, 6]],   # d = 1.4
        [[5.7, 6, 6], [8.4, 6, 6]],  # d = 2.7
        [[5.8, 6, 6], [7.7, 6, 6]],  # d = 1.9 ← drops
        [[5.0, 6, 6], [9.4, 6, 6]],  # d = 4.4
    ]
    energies = [-10, -9.5, -9.5, -10]
    traj = _path_traj(tmp_path, images_pos, energies, symbols=["O", "O"])
    result = validate_neb_path(str(traj), adsorbate_indices=[6, 7],
                                expected_monotonic_pairs=[(6, 7)],
                                max_jump=2.0, verbose=False)
    assert not result["valid"]
    assert any("monotonically" in i for i in result["issues"])
