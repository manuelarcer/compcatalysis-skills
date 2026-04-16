"""Unit tests for surface-generation/scripts/generate_slabs.py.

Uses real pymatgen calls on a tiny fcc Pt fixture — fast and deterministic.
"""

import pytest

import generate_slabs


def test_parse_bonds_none_returns_none():
    assert generate_slabs.parse_bonds(None) is None
    assert generate_slabs.parse_bonds([]) is None


def test_parse_bonds_parses_pair():
    from pymatgen.core import Element

    bonds = generate_slabs.parse_bonds(["Ti-O:2.0"])
    assert bonds == {(Element("Ti"), Element("O")): 2.0}


def test_parse_bonds_invalid_format_exits(capsys):
    with pytest.raises(SystemExit):
        generate_slabs.parse_bonds(["malformed"])
    assert "Invalid bond format" in capsys.readouterr().err


def test_load_structure_missing_file_exits(tmp_path, capsys):
    with pytest.raises(SystemExit):
        generate_slabs.load_structure(tmp_path / "nope.vasp")
    assert "Input file not found" in capsys.readouterr().err


def test_load_structure_roundtrip(pt_bulk_file):
    structure = generate_slabs.load_structure(pt_bulk_file)
    assert structure.num_sites == 4
    assert str(structure[0].specie) == "Pt"


def test_generate_slabs_for_miller_pt111(pt_bulk):
    slabs = generate_slabs.generate_slabs_for_miller(
        pt_bulk, miller_index=(1, 1, 1),
        min_slab_size=8.0, min_vacuum_size=12.0,
    )
    assert len(slabs) >= 1
    for slab in slabs:
        assert all(str(s.specie) == "Pt" for s in slab.sites)


def test_top_surface_signature_distinguishes_species(pt111_slab):
    sig = generate_slabs.get_top_surface_signature(pt111_slab)
    assert "Pt" in sig


def test_apply_selective_dynamics_fixes_bottom(pt111_slab):
    slab, n_fixed = generate_slabs.apply_selective_dynamics(pt111_slab, fix_bottom=0.5)
    sd = slab.site_properties["selective_dynamics"]
    assert n_fixed > 0
    assert n_fixed < len(sd)  # not everything fixed
    # fixed atoms → [False, False, False]; free → [True, True, True]
    assert any(flags == [False, False, False] for flags in sd)
    assert any(flags == [True, True, True] for flags in sd)


def test_apply_selective_dynamics_zero_skips(pt111_slab):
    result = generate_slabs.apply_selective_dynamics(pt111_slab, fix_bottom=0)
    # fix_bottom=0 returns the slab unchanged (not a tuple)
    assert result is pt111_slab


def test_slab_filename_format():
    assert generate_slabs.slab_filename("Pt", (1, 1, 1), 0) == "Pt_111_t0.vasp"
    assert generate_slabs.slab_filename("TiO2", (1, 1, 0), 2, fmt="cif") == "TiO2_110_t2.cif"


def test_miller_display():
    assert generate_slabs.miller_display((1, 1, 1)) == "(111)"
    assert generate_slabs.miller_str((1, 0, 0)) == "100"
