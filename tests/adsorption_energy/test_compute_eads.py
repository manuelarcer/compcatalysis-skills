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


def _traj(path, energy, *, mlip="uma-s-1p1", uma_task="oc20", with_meta=True):
    """Write a single-frame ASE traj with `energy` and (by default) a sibling
    compcat_run.json so the consistency check passes."""
    a = Atoms("Pt", positions=[[0, 0, 0]], cell=[5, 5, 5], pbc=True)
    a.calc = SinglePointCalculator(a, energy=energy)
    write(str(path), a)
    if with_meta:
        meta_path = path.parent / compute_eads.RUN_META_FILENAME
        meta_path.write_text(json.dumps({
            "mlip": mlip, "uma_task": uma_task,
            "fmax": 0.03, "optimizer": "fire",
        }))
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


def _mk_dir(tmp_path, name):
    d = tmp_path / name
    d.mkdir(parents=True, exist_ok=True)
    return d


def test_compute_one(tmp_path):
    slab = _traj(_mk_dir(tmp_path, "clean") / "opt.traj", -100.0)
    slab_ads = _traj(_mk_dir(tmp_path, "CO") / "opt.traj", -103.5)
    ref = _traj(_mk_dir(tmp_path, "CO_gas") / "opt.traj", -2.0)
    row = compute_eads.compute_one(slab, slab_ads, ref, label="CO_top_0", coeff=1.0)
    assert row["label"] == "CO_top_0"
    assert row["E_ads_eV"] == pytest.approx(-1.5)
    assert "CO_gas" in row["ref_summary"]


def test_compute_one_half_o2(tmp_path):
    slab = _traj(_mk_dir(tmp_path, "clean") / "opt.traj", -100.0)
    slab_ads = _traj(_mk_dir(tmp_path, "O") / "opt.traj", -101.5)
    o2 = _traj(_mk_dir(tmp_path, "O2_gas") / "opt.traj", -2.0)
    row = compute_eads.compute_one(slab, slab_ads, o2, label="O_fcc", coeff=0.5)
    # E_ads = -101.5 - (-100.0) - 0.5*(-2.0) = -101.5 + 100.0 + 1.0 = -0.5
    assert row["E_ads_eV"] == pytest.approx(-0.5)
    assert row["coeff"] == 0.5


# ---- task-head consistency check -------------------------------------------

def test_consistency_check_pass(tmp_path):
    slab = _traj(_mk_dir(tmp_path, "a") / "opt.traj", -100.0)
    ads = _traj(_mk_dir(tmp_path, "b") / "opt.traj", -101.0)
    ref = _traj(_mk_dir(tmp_path, "c") / "opt.traj", -2.0)
    ok, msg = compute_eads.check_run_consistency([slab, ads, ref])
    assert ok and "share" in msg


def test_consistency_check_fails_on_task_mismatch(tmp_path):
    slab = _traj(_mk_dir(tmp_path, "a") / "opt.traj", -100.0, uma_task="oc20")
    ads = _traj(_mk_dir(tmp_path, "b") / "opt.traj", -101.0, uma_task="omat")
    ref = _traj(_mk_dir(tmp_path, "c") / "opt.traj", -2.0, uma_task="oc20")
    ok, msg = compute_eads.check_run_consistency([slab, ads, ref])
    assert not ok and "inconsistent" in msg


def test_consistency_check_fails_on_mlip_mismatch(tmp_path):
    slab = _traj(_mk_dir(tmp_path, "a") / "opt.traj", -100.0, mlip="uma-s-1p1")
    ads = _traj(_mk_dir(tmp_path, "b") / "opt.traj", -101.0, mlip="mace")
    ref = _traj(_mk_dir(tmp_path, "c") / "opt.traj", -2.0, mlip="uma-s-1p1")
    ok, _ = compute_eads.check_run_consistency([slab, ads, ref])
    assert not ok


def test_consistency_check_warns_on_missing_metadata(tmp_path):
    slab = _traj(_mk_dir(tmp_path, "a") / "opt.traj", -100.0)
    ads = _traj(_mk_dir(tmp_path, "b") / "opt.traj", -101.0, with_meta=False)
    ref = _traj(_mk_dir(tmp_path, "c") / "opt.traj", -2.0)
    ok, msg = compute_eads.check_run_consistency([slab, ads, ref], strict=False)
    assert ok and "warning" in msg


def test_compute_one_raises_on_task_mismatch(tmp_path):
    slab = _traj(_mk_dir(tmp_path, "a") / "opt.traj", -100.0, uma_task="oc20")
    ads = _traj(_mk_dir(tmp_path, "b") / "opt.traj", -101.0, uma_task="omat")
    ref = _traj(_mk_dir(tmp_path, "c") / "opt.traj", -2.0, uma_task="oc20")
    with pytest.raises(ValueError, match="consistency"):
        compute_eads.compute_one(slab, ads, ref, label="x")


def test_compute_one_skips_check_when_disabled(tmp_path):
    slab = _traj(_mk_dir(tmp_path, "a") / "opt.traj", -100.0, uma_task="oc20")
    ads = _traj(_mk_dir(tmp_path, "b") / "opt.traj", -101.0, uma_task="omat")
    ref = _traj(_mk_dir(tmp_path, "c") / "opt.traj", -2.0, uma_task="oc20")
    row = compute_eads.compute_one(slab, ads, ref, label="x", consistency=False)
    assert "E_ads_eV" in row


# ---- composite OER references ----------------------------------------------

def test_expand_reference_named(tmp_path):
    table = {
        "OH": [{"traj": "H2O.traj", "coeff": 1.0},
               {"traj": "H2.traj", "coeff": -0.5}],
    }
    terms = compute_eads.expand_reference("OH", table)
    assert [(p.name, c) for p, c in terms] == [("H2O.traj", 1.0), ("H2.traj", -0.5)]


def test_expand_reference_path_fallback(tmp_path):
    terms = compute_eads.expand_reference("CO_gas/opt.traj", {})
    assert len(terms) == 1
    p, c = terms[0]
    assert str(p) == "CO_gas/opt.traj"
    assert c == 1.0


def test_compute_one_composite_OH(tmp_path):
    # Mock OH adsorption: E_slab=-100, E(slab+OH)=-105, E(H2O)=-14, E(H2)=-7
    # E_ads(OH) = -105 - (-100) - (1.0*-14 + -0.5*-7) = -5 - (-14 + 3.5) = -5 - (-10.5) = 5.5
    slab = _traj(_mk_dir(tmp_path, "a") / "opt.traj", -100.0)
    ads = _traj(_mk_dir(tmp_path, "b") / "opt.traj", -105.0)
    h2o = _traj(_mk_dir(tmp_path, "H2O") / "opt.traj", -14.0)
    h2 = _traj(_mk_dir(tmp_path, "H2") / "opt.traj", -7.0)
    row = compute_eads.compute_one_composite(
        slab, ads, [(h2o, 1.0), (h2, -0.5)], label="OH",
    )
    assert row["E_ads_eV"] == pytest.approx(5.5)
    assert row["coeff"] == "composite"
    assert "H2O" in row["ref_summary"] and "H2" in row["ref_summary"]


def test_compute_from_config_with_oer_references(tmp_path):
    slab = _traj(_mk_dir(tmp_path, "clean") / "opt.traj", -100.0)
    oh_ads = _traj(_mk_dir(tmp_path, "OH") / "opt.traj", -105.0)
    o_ads = _traj(_mk_dir(tmp_path, "O") / "opt.traj", -101.5)
    h2o = _traj(_mk_dir(tmp_path, "H2O") / "opt.traj", -14.0)
    h2 = _traj(_mk_dir(tmp_path, "H2") / "opt.traj", -7.0)
    o2 = _traj(_mk_dir(tmp_path, "O2") / "opt.traj", -2.0)

    cfg = {
        "slab": str(slab),
        "references": {
            "OH": [{"traj": str(h2o), "coeff": 1.0},
                   {"traj": str(h2),  "coeff": -0.5}],
            "O":  [{"traj": str(o2),  "coeff": 0.5}],
        },
        "calculations": [
            {"label": "OH_top_0", "slab_ads": str(oh_ads), "ref": "OH"},
            {"label": "O_fcc",    "slab_ads": str(o_ads),  "ref": "O"},
        ],
    }
    rows = compute_eads.compute_from_config(cfg)
    assert {r["label"] for r in rows} == {"OH_top_0", "O_fcc"}
    by_label = {r["label"]: r for r in rows}
    assert by_label["OH_top_0"]["E_ads_eV"] == pytest.approx(5.5)
    assert by_label["O_fcc"]["E_ads_eV"] == pytest.approx(-0.5)


def test_compute_from_config(tmp_path):
    slab = _traj(_mk_dir(tmp_path, "clean") / "opt.traj", -100.0)
    co_ads = _traj(_mk_dir(tmp_path, "CO") / "opt.traj", -103.0)
    o_ads = _traj(_mk_dir(tmp_path, "O") / "opt.traj", -101.5)
    co_ref = _traj(_mk_dir(tmp_path, "CO_gas") / "opt.traj", -2.0)
    o2_ref = _traj(_mk_dir(tmp_path, "O2_gas") / "opt.traj", -2.0)

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
    slab = _traj(_mk_dir(tmp_path, "clean") / "opt.traj", -100.0)
    ads = _traj(_mk_dir(tmp_path, "CO") / "opt.traj", -103.0)
    ref = _traj(_mk_dir(tmp_path, "CO_gas") / "opt.traj", -2.0)
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


def test_print_table_handles_composite_coeff(capsys):
    """Composite-reference rows have coeff='composite' (a string), not a float.
    print_table must not try to format it as a float."""
    rows = [
        {"label": "OH", "E_ads_eV": -4.0, "coeff": "composite", "slab_ads": "OH/opt.traj"},
        {"label": "CO", "E_ads_eV": -1.0, "coeff": 1.0, "slab_ads": "CO/opt.traj"},
    ]
    compute_eads.print_table(rows)
    out = capsys.readouterr().out
    assert "composite" in out
    assert "1.00" in out  # the float coeff still formats correctly


def test_main_config(tmp_path, monkeypatch):
    slab = _traj(_mk_dir(tmp_path, "clean") / "opt.traj", -100.0)
    ads = _traj(_mk_dir(tmp_path, "CO") / "opt.traj", -103.0)
    ref = _traj(_mk_dir(tmp_path, "CO_gas") / "opt.traj", -2.0)
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
