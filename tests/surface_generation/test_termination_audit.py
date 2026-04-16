"""Unit tests for surface-generation/scripts/termination_audit.py."""

import pytest

import termination_audit as ta


def test_load_structure_missing_file_exits(tmp_path, capsys):
    with pytest.raises(SystemExit):
        ta.load_structure(tmp_path / "nope.vasp")
    assert "File not found" in capsys.readouterr().err


def test_get_ouc_and_planes_pt111(pt_bulk):
    ouc, sorted_z, sorted_species, planes, plane_z, plane_species, gaps = \
        ta.get_ouc_and_planes(pt_bulk, (1, 1, 1))

    assert ouc.num_sites == len(sorted_z) == len(sorted_species)
    assert len(plane_z) == len(planes) == len(plane_species)
    assert len(planes) >= 1
    # Gaps list length matches plane count (includes periodic wrap)
    assert len(gaps) == len(planes)
    # Gaps sum to c-axis length (modulo the periodic wrap construction)
    import numpy as np
    assert np.isclose(sum(gaps), ouc.lattice.c, atol=1e-6)
    # Plane z values are monotonically increasing
    assert all(plane_z[i] < plane_z[i + 1] for i in range(len(plane_z) - 1))
    # All species should be Pt
    assert all("Pt" in sp for sp in plane_species)


def test_sweep_ftol_pt111(pt_bulk):
    results = ta.sweep_ftol(pt_bulk, (1, 1, 1), ftol_min=0.05, ftol_max=0.2, n_steps=3)
    assert len(results) == 3
    for r in results:
        assert set(r) == {"ftol", "n_terminations", "shifts"}
        assert r["n_terminations"] >= 1
        assert len(r["shifts"]) == r["n_terminations"]


def test_audit_returns_summary(pt_bulk, capsys):
    summary = ta.audit(
        pt_bulk, (1, 1, 1),
        ftol_min=0.05, ftol_max=0.2, n_ftol_steps=3,
    )
    assert summary["miller_index"] == (1, 1, 1)
    assert summary["n_planes"] >= 1
    assert summary["max_terminations"] >= 1
    assert len(summary["sweep"]) == 3
    # Function prints to stdout — sanity-check it ran
    out = capsys.readouterr().out
    assert "TERMINATION AUDIT" in out and "Summary" in out
