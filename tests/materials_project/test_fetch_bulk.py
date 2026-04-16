"""Unit tests for materials-project/scripts/fetch_bulk.py.

The mp-api client is mocked — these tests never hit the network. A single
live smoke test is marked @pytest.mark.network for on-demand runs.
"""

import json
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

import fetch_bulk


def test_get_api_key_from_env(monkeypatch):
    monkeypatch.setenv("MP_API_KEY", "sentinel-key")
    assert fetch_bulk.get_api_key() == "sentinel-key"


def test_get_api_key_cli_overrides_env(monkeypatch):
    monkeypatch.setenv("MP_API_KEY", "env-key")
    assert fetch_bulk.get_api_key("cli-key") == "cli-key"


def test_get_api_key_missing_exits(monkeypatch, capsys):
    monkeypatch.delenv("MP_API_KEY", raising=False)
    with pytest.raises(SystemExit):
        fetch_bulk.get_api_key()
    assert "No Materials Project API key" in capsys.readouterr().err


def test_write_structure_vasp_and_json(pt_bulk, tmp_path):
    vasp_path = tmp_path / "Pt.vasp"
    fetch_bulk.write_structure(pt_bulk, vasp_path, fmt="vasp")
    assert vasp_path.exists() and vasp_path.stat().st_size > 0
    assert "Pt" in vasp_path.read_text()

    json_path = tmp_path / "Pt.json"
    fetch_bulk.write_structure(pt_bulk, json_path, fmt="json")
    data = json.loads(json_path.read_text())
    assert "lattice" in data and "sites" in data


def test_write_structure_rejects_unknown_format(pt_bulk, tmp_path):
    with pytest.raises(ValueError, match="Unknown format"):
        fetch_bulk.write_structure(pt_bulk, tmp_path / "x.bin", fmt="bin")


def test_to_conventional_returns_structure(pt_bulk):
    conv = fetch_bulk.to_conventional(pt_bulk)
    # fcc Pt conventional cell has 4 atoms
    assert conv.num_sites == 4
    assert all(str(s.specie) == "Pt" for s in conv.sites)


def _fake_doc(pt_structure, mp_id="mp-126", ehull=0.0):
    return SimpleNamespace(
        material_id=mp_id,
        formula_pretty="Pt",
        structure=pt_structure,
        symmetry=SimpleNamespace(symbol="Fm-3m", number=225),
        energy_above_hull=ehull,
        formation_energy_per_atom=0.0,
        band_gap=0.0,
        nsites=pt_structure.num_sites,
    )


def test_fetch_by_formula_sorts_by_ehull(pt_bulk):
    unstable = _fake_doc(pt_bulk, mp_id="mp-unstable", ehull=0.02)
    stable = _fake_doc(pt_bulk, mp_id="mp-stable", ehull=0.0)

    fake_mpr = MagicMock()
    fake_mpr.__enter__.return_value = fake_mpr
    fake_mpr.__exit__.return_value = False
    fake_mpr.materials.summary.search.return_value = [unstable, stable]

    with patch("mp_api.client.MPRester", return_value=fake_mpr):
        docs = fetch_bulk.fetch_by_formula("fake-key", "Pt")

    assert [d.material_id for d in docs] == ["mp-stable", "mp-unstable"]


def test_fetch_by_mpid_exits_on_empty(pt_bulk, capsys):
    fake_mpr = MagicMock()
    fake_mpr.__enter__.return_value = fake_mpr
    fake_mpr.__exit__.return_value = False
    fake_mpr.materials.summary.search.return_value = []

    with patch("mp_api.client.MPRester", return_value=fake_mpr):
        with pytest.raises(SystemExit):
            fetch_bulk.fetch_by_mpid("fake-key", "mp-does-not-exist")
    assert "No structure found" in capsys.readouterr().err


@pytest.mark.network
def test_fetch_pt_live():
    """Smoke test against the real MP API. Skipped unless MP_API_KEY is set."""
    import os

    if not os.environ.get("MP_API_KEY"):
        pytest.skip("MP_API_KEY not set")
    docs = fetch_bulk.fetch_by_formula(os.environ["MP_API_KEY"], "Pt")
    assert docs and docs[0].formula_pretty == "Pt"
