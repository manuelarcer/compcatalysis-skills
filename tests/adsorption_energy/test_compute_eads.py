"""Unit tests for adsorption-energy/scripts/compute_eads.py.

Trajectories are synthesized with ASE's SinglePointCalculator so the tests
are hermetic — no MLIP, no real DFT.
"""

import csv
import json
import sys

import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write

import compute_eads


def _traj(path, energy):
    a = Atoms("Pt", positions=[[0, 0, 0]], cell=[5, 5, 5], pbc=True)
    a.calc = SinglePointCalculator(a, energy=energy)
    write(str(path), a)
    return path


def test_adsorption_energy_math():
    # E_ads = E(slab+ads) - E(slab) - coeff*E(ref)
    assert compute_eads.adsorption_energy(-100.0, -103.0, -2.0, 1.0) == pytest.approx(-1.0)
    assert compute_eads.adsorption_energy(-100.0, -101.0, -2.0, 0.5) == pytest.approx(0.0)


def test_read_final_energy(tmp_path):
    p = _traj(tmp_path / "x.traj", energy=-7.5)
    assert compute_eads.read_final_energy(p) == pytest.approx(-7.5)


def test_read_final_energy_missing(tmp_path):
    with pytest.raises(FileNotFoundError):
        compute_eads.read_final_energy(tmp_path / "nope.traj")


def test_compute_one(tmp_path):
    slab = _traj(tmp_path / "slab.traj", -100.0)
    slab_ads = _traj(tmp_path / "ads.traj", -103.5)
    ref = _traj(tmp_path / "ref.traj", -2.0)
    row = compute_eads.compute_one(slab, slab_ads, ref, label="CO_top_0", coeff=1.0)
    assert row["label"] == "CO_top_0"
    assert row["E_ads_eV"] == pytest.approx(-1.5)


def test_compute_one_half_o2(tmp_path):
    slab = _traj(tmp_path / "slab.traj", -100.0)
    slab_ads = _traj(tmp_path / "Oads.traj", -101.5)
    o2 = _traj(tmp_path / "O2.traj", -2.0)
    row = compute_eads.compute_one(slab, slab_ads, o2, label="O_fcc", coeff=0.5)
    # E_ads = -101.5 - (-100.0) - 0.5*(-2.0) = -101.5 + 100.0 + 1.0 = -0.5
    assert row["E_ads_eV"] == pytest.approx(-0.5)
    assert row["coeff"] == 0.5


def test_compute_from_config(tmp_path):
    slab = _traj(tmp_path / "slab.traj", -100.0)
    co_ads = _traj(tmp_path / "co.traj", -103.0)
    o_ads = _traj(tmp_path / "o.traj", -101.5)
    co_ref = _traj(tmp_path / "co_gas.traj", -2.0)
    o2_ref = _traj(tmp_path / "o2.traj", -2.0)

    cfg = {
        "slab": str(slab),
        "calculations": [
            {"label": "CO", "slab_ads": str(co_ads), "ref": str(co_ref)},
            {"label": "O", "slab_ads": str(o_ads), "ref": str(o2_ref), "coeff": 0.5},
        ],
    }
    rows = compute_eads.compute_from_config(cfg)
    assert [r["label"] for r in rows] == ["CO", "O"]
    assert rows[0]["E_ads_eV"] == pytest.approx(-1.0)
    assert rows[1]["E_ads_eV"] == pytest.approx(-0.5)


def test_compute_from_config_missing_keys():
    with pytest.raises(ValueError, match="must have"):
        compute_eads.compute_from_config({"slab": "x"})


def test_write_csv_overwrite_then_append(tmp_path):
    rows1 = [{"label": "a", "E_ads_eV": -1.0, "coeff": 1.0}]
    rows2 = [{"label": "b", "E_ads_eV": -0.5, "coeff": 0.5}]
    out = tmp_path / "eads.csv"
    compute_eads.write_csv(rows1, out)
    compute_eads.write_csv(rows2, out, append=True)
    with open(out) as f:
        loaded = list(csv.DictReader(f))
    assert [r["label"] for r in loaded] == ["a", "b"]


def test_main_single(tmp_path, monkeypatch, capsys):
    slab = _traj(tmp_path / "slab.traj", -100.0)
    ads = _traj(tmp_path / "ads.traj", -103.0)
    ref = _traj(tmp_path / "ref.traj", -2.0)
    out = tmp_path / "eads.csv"

    monkeypatch.setattr(sys, "argv", [
        "compute_eads.py",
        "--slab", str(slab),
        "--slab-ads", str(ads),
        "--ads-ref", str(ref),
        "--label", "CO_top_0",
        "--output", str(out),
    ])
    compute_eads.main()
    rows = list(csv.DictReader(open(out)))
    assert len(rows) == 1
    assert rows[0]["label"] == "CO_top_0"
    assert float(rows[0]["E_ads_eV"]) == pytest.approx(-1.0)


def test_main_config(tmp_path, monkeypatch):
    slab = _traj(tmp_path / "slab.traj", -100.0)
    ads = _traj(tmp_path / "co.traj", -103.0)
    ref = _traj(tmp_path / "co_gas.traj", -2.0)
    cfg_path = tmp_path / "cfg.json"
    cfg_path.write_text(json.dumps({
        "slab": str(slab),
        "calculations": [
            {"label": "CO", "slab_ads": str(ads), "ref": str(ref)},
        ],
    }))
    out = tmp_path / "eads.csv"
    monkeypatch.setattr(sys, "argv", [
        "compute_eads.py", "--config", str(cfg_path), "--output", str(out),
    ])
    compute_eads.main()
    rows = list(csv.DictReader(open(out)))
    assert rows[0]["label"] == "CO"
