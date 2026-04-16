"""Deterministic artifact checks for skill behavioral evals.

Each check receives (artifact_path: Path, arg) and returns (passed, message).
Register new checks by adding to CHECKS. Scenarios reference them by name.
"""

from pathlib import Path


def check_exists(path: Path, _arg=None):
    if path.exists() and path.stat().st_size > 0:
        return True, f"exists ({path.stat().st_size} B)"
    return False, f"missing or empty: {path}"


def check_pymatgen_parseable(path: Path, _arg=None):
    from pymatgen.core import Structure
    try:
        s = Structure.from_file(str(path))
        return True, f"pymatgen parsed {s.num_sites} sites"
    except Exception as e:
        return False, f"pymatgen failed: {e}"


def check_ase_parseable(path: Path, _arg=None):
    from ase.io import read
    try:
        atoms = read(str(path))
        return True, f"ase parsed {len(atoms)} atoms"
    except Exception as e:
        return False, f"ase failed: {e}"


def check_composition(path: Path, expected: str):
    from pymatgen.core import Structure
    try:
        s = Structure.from_file(str(path))
        actual = s.composition.reduced_formula
        if actual == expected:
            return True, f"composition {actual}"
        return False, f"expected {expected}, got {actual}"
    except Exception as e:
        return False, f"composition check failed: {e}"


def check_contains_species(path: Path, species: str):
    from ase.io import read
    try:
        atoms = read(str(path))
        if species in atoms.get_chemical_symbols():
            return True, f"contains {species}"
        return False, f"no {species} in {sorted(set(atoms.get_chemical_symbols()))}"
    except Exception as e:
        return False, f"species check failed: {e}"


def check_has_selective_dynamics(path: Path, _arg=None):
    """Verify a VASP POSCAR has selective dynamics flags set."""
    try:
        text = Path(path).read_text()
        if "Selective" in text or "selective" in text:
            return True, "selective dynamics present"
        return False, "no selective dynamics block"
    except Exception as e:
        return False, f"read failed: {e}"


CHECKS = {
    "exists": check_exists,
    "pymatgen_parseable": check_pymatgen_parseable,
    "ase_parseable": check_ase_parseable,
    "composition": check_composition,
    "contains_species": check_contains_species,
    "has_selective_dynamics": check_has_selective_dynamics,
}


def grade_artifact(workdir: Path, spec: dict):
    """Grade a single expected artifact against its checks."""
    path_glob = spec.get("path") or spec.get("path_glob")
    matches = sorted(workdir.glob(path_glob)) if path_glob else []

    min_count = spec.get("min_count", 1)
    result = {
        "path": path_glob,
        "matches": [str(p.relative_to(workdir)) for p in matches],
        "min_count": min_count,
        "count_ok": len(matches) >= min_count,
        "checks": [],
    }

    target = matches[0] if matches else (workdir / path_glob)
    for check in spec.get("checks", []):
        if isinstance(check, str):
            name, arg = check, None
        elif isinstance(check, dict) and len(check) == 1:
            name, arg = next(iter(check.items()))
        else:
            result["checks"].append({"name": str(check), "pass": False,
                                      "message": "malformed check spec"})
            continue

        fn = CHECKS.get(name)
        if fn is None:
            result["checks"].append({"name": name, "pass": False,
                                      "message": f"unknown check '{name}'"})
            continue

        if not matches:
            result["checks"].append({"name": name, "pass": False,
                                      "message": "no matching artifact"})
            continue

        try:
            passed, msg = fn(target, arg)
        except Exception as e:
            passed, msg = False, f"check raised: {e}"
        result["checks"].append({"name": name, "pass": passed, "message": msg})

    result["pass"] = result["count_ok"] and all(c["pass"] for c in result["checks"])
    return result


def grade_scenario(workdir: Path, expect: dict):
    """Grade every expected artifact; return aggregate pass + details."""
    artifacts = [grade_artifact(workdir, spec) for spec in expect.get("artifacts", [])]
    return {
        "pass": all(a["pass"] for a in artifacts) if artifacts else False,
        "artifacts": artifacts,
    }
