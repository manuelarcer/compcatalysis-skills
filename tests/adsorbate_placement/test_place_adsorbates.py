"""Unit tests for adsorbate-placement/scripts/place_adsorbates.py.

Uses real ASE on a tiny Pt(111) fixture slab.
"""

import numpy as np
import pytest
from ase import Atom, Atoms
from ase.constraints import FixAtoms
from ase.io import read, write

import place_adsorbates as pa


def test_estimate_height_sums_covalent_radii():
    from ase.data import atomic_numbers, covalent_radii

    expected = covalent_radii[atomic_numbers["Pt"]] + covalent_radii[atomic_numbers["O"]]
    assert pa.estimate_height("Pt", "O") == pytest.approx(expected)


def test_build_O_returns_single_atom():
    atoms, binding = pa.build_O()
    assert binding == "O"
    assert len(atoms) == 1
    assert atoms[0][0] == "O"
    np.testing.assert_allclose(atoms[0][1], [0, 0, 0])


def test_build_OH_geometry():
    atoms, binding = pa.build_OH()
    assert binding == "O"
    assert [a[0] for a in atoms] == ["O", "H"]
    np.testing.assert_allclose(atoms[0][1], [0, 0, 0])
    # H is directly above O at ~0.98 Å
    assert atoms[1][1][2] == pytest.approx(0.98)


def test_build_OOH_geometry():
    atoms, binding = pa.build_OOH()
    assert binding == "O"
    assert [a[0] for a in atoms] == ["O", "O", "H"]
    np.testing.assert_allclose(atoms[0][1], [0, 0, 0])
    # second O is tilted (has non-zero x component)
    assert atoms[1][1][0] > 0


def _slab_to_ase(pt111_slab, tmp_path):
    path = tmp_path / "slab.vasp"
    pt111_slab.to(filename=str(path), fmt="poscar")
    return read(str(path))


def test_find_surface_sites_returns_pt_sites(pt111_slab, tmp_path):
    slab = _slab_to_ase(pt111_slab, tmp_path)
    sites = pa.find_surface_sites(slab, depth=2.0, site_species=["Pt"])
    assert "Pt" in sites and len(sites["Pt"]) >= 1
    first = sites["Pt"][0]
    assert first["label"].startswith("Pt_top_")
    assert first["symbol"] == "Pt"


def test_find_surface_sites_filters_by_species(pt111_slab, tmp_path):
    slab = _slab_to_ase(pt111_slab, tmp_path)
    sites = pa.find_surface_sites(slab, depth=2.0, site_species=["Au"])
    assert sites == {}  # no Au on a Pt slab


def test_place_adsorbate_appends_atoms_and_preserves_constraints(pt111_slab, tmp_path):
    slab = _slab_to_ase(pt111_slab, tmp_path)
    fixed_indices = [0, 1]
    slab.set_constraint(FixAtoms(indices=fixed_indices))

    ads_atoms, _ = pa.build_OH()
    placed = pa.place_adsorbate(slab, site_index=len(slab) - 1, adsorbate_atoms=ads_atoms, height=2.0)

    assert len(placed) == len(slab) + 2  # O + H appended
    # Constraint preserved, only original indices fixed
    surviving = set()
    for c in placed.constraints:
        if isinstance(c, FixAtoms):
            surviving.update(c.index)
    assert surviving == set(fixed_indices)

    # New atoms sit above the site
    site_z = slab.positions[-1, 2]
    assert placed.positions[-2, 2] == pytest.approx(site_z + 2.0)  # O at base
    assert placed.positions[-1, 2] > placed.positions[-2, 2]        # H above O


def test_load_adsorbate_from_file_recenters(tmp_path):
    mol = Atoms("OH", positions=[[1.0, 2.0, 3.0], [1.0, 2.0, 3.98]])
    path = tmp_path / "OH.xyz"
    write(str(path), mol)

    atoms, binding = pa.load_adsorbate_from_file(path, binding_atom_index=0)

    assert binding == "O"
    # Binding atom at origin
    np.testing.assert_allclose(atoms[0][1], [0, 0, 0], atol=1e-6)
    # Second atom retains relative offset
    np.testing.assert_allclose(atoms[1][1], [0, 0, 0.98], atol=1e-6)


def test_load_adsorbate_flips_if_pointing_down(tmp_path):
    # Binding atom at top; other atom below it → should be flipped upward
    mol = Atoms("OH", positions=[[0.0, 0.0, 1.0], [0.0, 0.0, 0.0]])
    path = tmp_path / "OH_down.xyz"
    write(str(path), mol)

    atoms, _ = pa.load_adsorbate_from_file(path, binding_atom_index=0)
    # After flip, non-binding atom sits at positive z
    assert atoms[1][1][2] > 0


def test_load_adsorbate_out_of_range_exits(tmp_path, capsys):
    mol = Atoms("O", positions=[[0.0, 0.0, 0.0]])
    path = tmp_path / "O.xyz"
    write(str(path), mol)
    with pytest.raises(SystemExit):
        pa.load_adsorbate_from_file(path, binding_atom_index=5)
    assert "out of range" in capsys.readouterr().err


def test_generate_all_placements_writes_tree(pt111_slab, tmp_path):
    slab = _slab_to_ase(pt111_slab, tmp_path)
    ads_atoms, binding = pa.build_OH()
    specs = [{"name": "OH", "atoms": ads_atoms, "binding_symbol": binding}]

    out_dir = tmp_path / "out"
    results = pa.generate_all_placements(
        slab, adsorbate_specs=specs, site_species=["Pt"],
        output_dir=out_dir, depth=2.0,
    )

    assert len(results) >= 1
    for r in results:
        assert r["adsorbate"] == "OH"
        assert r["site_species"] == "Pt"
        assert r["site_label"].startswith("Pt_top_")
        p = tmp_path / "out" / "OH" / r["site_label"] / "input.vasp"
        assert p.exists() and p.stat().st_size > 0
