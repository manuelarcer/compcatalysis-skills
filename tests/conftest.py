"""Shared fixtures for compcatalysis-skills tests.

Structures are built programmatically (pymatgen/ase) rather than committed as
binary fixtures — keeps the repo small and the fixtures self-documenting.
"""

import pytest


@pytest.fixture
def pt_bulk():
    """Minimal fcc Pt conventional bulk (4 atoms)."""
    from pymatgen.core import Lattice, Structure

    a = 3.924
    lattice = Lattice.cubic(a)
    species = ["Pt"] * 4
    coords = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ]
    return Structure(lattice, species, coords)


@pytest.fixture
def pt_bulk_file(pt_bulk, tmp_path):
    """fcc Pt bulk written as POSCAR; returns the path."""
    path = tmp_path / "Pt_bulk.vasp"
    pt_bulk.to(filename=str(path), fmt="poscar")
    return path


@pytest.fixture
def pt111_slab(pt_bulk):
    """A single Pt(111) slab with selective dynamics applied on the bottom half."""
    from pymatgen.core.surface import SlabGenerator

    gen = SlabGenerator(
        pt_bulk,
        miller_index=(1, 1, 1),
        min_slab_size=8.0,
        min_vacuum_size=12.0,
        center_slab=True,
        lll_reduce=True,
        primitive=True,
    )
    slabs = gen.get_slabs(filter_out_sym_slabs=False)
    return slabs[0]


@pytest.fixture
def pt111_slab_file(pt111_slab, tmp_path):
    """Pt(111) slab written as POSCAR; returns the path."""
    path = tmp_path / "Pt111_slab.vasp"
    pt111_slab.to(filename=str(path), fmt="poscar")
    return path
