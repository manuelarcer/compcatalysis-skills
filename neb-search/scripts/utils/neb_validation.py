"""
NEB validation utilities for surface reactions.

Two main functions:
  validate_neb_setup()  - Check IS/FS pair before running NEB
  validate_neb_path()   - Check NEB trajectory after running NEB
"""

import numpy as np
from ase.geometry import find_mic


def validate_neb_setup(initial, final, adsorbate_indices, max_disp=3.0,
                       check_crossing=True, verbose=True):
    """
    Validate IS/FS pair for NEB: check atom mapping, displacements, and crossings.

    Parameters
    ----------
    initial : ase.Atoms
        Initial state (optimized)
    final : ase.Atoms
        Final state (optimized)
    adsorbate_indices : list of int
        Indices of adsorbate atoms (the atoms that move during the reaction)
    max_disp : float
        Maximum allowed displacement per adsorbate atom (A). Default 3.0
    check_crossing : bool
        If True, check whether adsorbate atoms cross paths (swap)
    verbose : bool
        Print detailed report

    Returns
    -------
    dict with keys:
        'valid' : bool - overall pass/fail
        'issues' : list of str - description of each issue found
        'displacements' : dict - per-atom displacement info
        'crossing' : bool or None - True if atoms cross
    """
    issues = []
    result = {'valid': True, 'issues': issues, 'displacements': {}, 'crossing': None}

    # --- Basic checks ---
    if len(initial) != len(final):
        issues.append(f"Atom count mismatch: IS={len(initial)}, FS={len(final)}")
        result['valid'] = False
        return result

    sym_i = initial.get_chemical_symbols()
    sym_f = final.get_chemical_symbols()
    if sym_i != sym_f:
        issues.append(f"Chemical symbols mismatch at indices: "
                      f"{[i for i in range(len(sym_i)) if sym_i[i] != sym_f[i]]}")
        result['valid'] = False
        return result

    # --- Per-atom displacements (MIC) ---
    all_disp = final.get_positions() - initial.get_positions()
    mic_disp, mic_dist = find_mic(all_disp, initial.cell, pbc=True)

    if verbose:
        print("=" * 65)
        print("NEB Setup Validation")
        print("=" * 65)
        print(f"\nAdsorbate atoms: {adsorbate_indices}")
        print(f"{'Atom':>6} {'Elem':>4} {'IS_x':>7} {'IS_y':>7} {'IS_z':>7} "
              f"{'FS_x':>7} {'FS_y':>7} {'FS_z':>7} {'|d|':>6}")
        print("-" * 65)

    for idx in adsorbate_indices:
        pi = initial.positions[idx]
        pf = final.positions[idx]
        d = mic_dist[idx]
        elem = sym_i[idx]

        result['displacements'][idx] = {
            'element': elem,
            'is_pos': pi.copy(),
            'fs_pos': pf.copy(),
            'mic_disp': mic_disp[idx].copy(),
            'distance': float(d),
        }

        if verbose:
            print(f"{idx:>6} {elem:>4} {pi[0]:>7.3f} {pi[1]:>7.3f} {pi[2]:>7.3f} "
                  f"{pf[0]:>7.3f} {pf[1]:>7.3f} {pf[2]:>7.3f} {d:>6.3f}")

        if d > max_disp:
            issues.append(f"Atom {idx} ({elem}) displacement {d:.3f} A > max {max_disp} A "
                          f"-- excessive diffusion, not a direct reaction path")

    # --- Check for adsorbate crossings (atom swap) ---
    if check_crossing and len(adsorbate_indices) >= 2:
        ads_idx = list(adsorbate_indices)
        n = len(ads_idx)

        # Build cost matrix: distance from IS position of atom i to FS position of atom j
        cost = np.zeros((n, n))
        for i_row, idx_i in enumerate(ads_idx):
            for j_col, idx_j in enumerate(ads_idx):
                disp = final.positions[idx_j] - initial.positions[idx_i]
                _, d_mic = find_mic(disp.reshape(1, 3), initial.cell, pbc=True)
                cost[i_row, j_col] = d_mic[0]

        # Check if the identity mapping (i->i) is optimal
        identity_cost = sum(cost[i, i] for i in range(n))

        # For 2 atoms, just check the swap
        if n == 2:
            swap_cost = cost[0, 1] + cost[1, 0]
            if swap_cost < identity_cost:
                result['crossing'] = True
                issues.append(
                    f"ATOM SWAP DETECTED: atoms {ads_idx[0]} and {ads_idx[1]} "
                    f"cross paths. Identity cost={identity_cost:.2f} A, "
                    f"swap cost={swap_cost:.2f} A. "
                    f"Swap the FS atom positions to fix.")
            else:
                result['crossing'] = False
        else:
            # For >2 atoms, use Hungarian algorithm
            try:
                from scipy.optimize import linear_sum_assignment
                row_ind, col_ind = linear_sum_assignment(cost)
                optimal_cost = cost[row_ind, col_ind].sum()
                if not np.array_equal(col_ind, np.arange(n)):
                    result['crossing'] = True
                    mapping = {ads_idx[r]: ads_idx[c] for r, c in zip(row_ind, col_ind)}
                    issues.append(
                        f"ATOM REMAPPING NEEDED: optimal mapping {mapping}, "
                        f"identity cost={identity_cost:.2f}, optimal cost={optimal_cost:.2f}")
                else:
                    result['crossing'] = False
            except ImportError:
                result['crossing'] = None
                issues.append("scipy not available for optimal assignment check")

        if verbose:
            print(f"\n--- Crossing check ---")
            print(f"Cost matrix (IS_row -> FS_col distances):")
            header = "      " + "".join(f"  FS_{ads_idx[j]:>3}" for j in range(n))
            print(header)
            for i_row in range(n):
                row_str = f"IS_{ads_idx[i_row]:>3}"
                for j_col in range(n):
                    marker = " *" if i_row == j_col else "  "
                    row_str += f" {cost[i_row, j_col]:>5.2f}{marker}"
                print(row_str)
            if result['crossing']:
                print(">>> CROSSING DETECTED - atoms swap positions!")
            else:
                print(">>> No crossing - identity mapping is optimal")

    # --- Check displacement vectors for dissociation vs diffusion ---
    if len(adsorbate_indices) >= 2:
        disps = [mic_disp[i] for i in adsorbate_indices]
        for i in range(len(disps)):
            for j in range(i + 1, len(disps)):
                xy_i = disps[i][:2]
                xy_j = disps[j][:2]
                norm_i = np.linalg.norm(xy_i)
                norm_j = np.linalg.norm(xy_j)
                idx_i = adsorbate_indices[i]
                idx_j = adsorbate_indices[j]
                if norm_i > 0.1 and norm_j > 0.1:
                    cos_angle = np.dot(xy_i, xy_j) / (norm_i * norm_j)
                    angle = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))
                    if verbose:
                        print(f"\nAngle between displacements of atoms {idx_i} and {idx_j}: "
                              f"{angle:.1f} deg")
                        if angle < 30:
                            print("  (parallel = both moving same direction = co-diffusion)")
                        elif angle > 120:
                            print("  (antiparallel = moving apart = dissociation, good)")
                        else:
                            print("  (oblique = mixed character)")
                    if angle < 60:
                        issues.append(
                            f"Atoms {idx_i} and {idx_j} move in similar directions "
                            f"(angle={angle:.0f} deg) -- looks like diffusion, "
                            f"not a direct reaction. Consider choosing IS/FS as a "
                            f"matched pair where the IS sits between the FS sites.")

                # Check displacement asymmetry
                ratio = max(norm_i, norm_j) / max(min(norm_i, norm_j), 0.01)
                if ratio > 3.0 and max(norm_i, norm_j) > 1.5:
                    issues.append(
                        f"Asymmetric displacements: atom {idx_i} moves {norm_i:.2f} A "
                        f"vs atom {idx_j} moves {norm_j:.2f} A (ratio {ratio:.1f}x). "
                        f"The IS may not be centered between the FS sites.")

    # --- Surface atom check ---
    substrate_disps = [mic_dist[i] for i in range(len(initial))
                       if i not in adsorbate_indices]
    max_substrate = max(substrate_disps) if substrate_disps else 0
    if max_substrate > 0.3:
        issues.append(
            f"Large substrate relaxation: max displacement {max_substrate:.3f} A. "
            f"This suggests the IS and FS have very different surface geometries.")
    if verbose:
        print(f"\nMax substrate atom displacement: {max_substrate:.4f} A")

    result['valid'] = len(issues) == 0
    if verbose:
        print(f"\n{'PASS' if result['valid'] else 'FAIL'}: "
              f"{len(issues)} issue(s) found")
        for issue in issues:
            print(f"  - {issue}")

    return result


def validate_neb_path(images, adsorbate_indices, expected_monotonic_pairs=None,
                      max_jump=1.0, verbose=True):
    """
    Validate a NEB path: check smoothness, crossings, and monotonicity.

    Parameters
    ----------
    images : list of ase.Atoms or str
        List of NEB images, or path to A2B.traj file
    adsorbate_indices : list of int
        Indices of adsorbate atoms
    expected_monotonic_pairs : list of tuple, optional
        Pairs of atom indices (i, j) where d(i,j) should change monotonically.
        E.g., [(36, 37)] for O-O dissociation.
    max_jump : float
        Maximum allowed displacement of any adsorbate atom between adjacent images (A)
    verbose : bool
        Print detailed table

    Returns
    -------
    dict with keys:
        'valid' : bool
        'issues' : list of str
        'table' : list of dict (per-image data)
        'max_jump' : float (largest inter-image displacement found)
    """
    from ase.io.trajectory import Trajectory

    if isinstance(images, str):
        images = Trajectory(images)

    issues = []
    n_images = len(images)
    table = []

    if verbose:
        print("=" * 75)
        print("NEB Path Validation")
        print("=" * 75)

    # --- Build per-image data ---
    for i, img in enumerate(images):
        row = {'image': i}

        # Energies
        try:
            e = img.get_potential_energy()
            row['energy'] = e
        except Exception:
            row['energy'] = None

        # Adsorbate positions
        for idx in adsorbate_indices:
            row[f'pos_{idx}'] = img.positions[idx].copy()

        # Inter-adsorbate distances
        if expected_monotonic_pairs:
            for (a, b) in expected_monotonic_pairs:
                d = img.get_distance(a, b, mic=True)
                row[f'd_{a}_{b}'] = d

        table.append(row)

    # --- Relative energies ---
    if table[0].get('energy') is not None:
        e0 = table[0]['energy']
        for row in table:
            row['rel_e'] = row['energy'] - e0

    # --- Check inter-image jumps ---
    max_jump_found = 0
    for i in range(1, n_images):
        for idx in adsorbate_indices:
            disp = images[i].positions[idx] - images[i - 1].positions[idx]
            mic_d, mic_dist = find_mic(disp.reshape(1, 3), images[0].cell, pbc=True)
            jump = mic_dist[0]
            if jump > max_jump_found:
                max_jump_found = jump
            if jump > max_jump:
                issues.append(
                    f"Image {i}: atom {idx} jumps {jump:.3f} A from image {i-1} "
                    f"(max allowed: {max_jump} A)")

    # --- Check for adsorbate crossings along path ---
    if len(adsorbate_indices) >= 2:
        for i in range(1, n_images):
            for a_idx, idx_a in enumerate(adsorbate_indices):
                for b_idx, idx_b in enumerate(adsorbate_indices):
                    if a_idx >= b_idx:
                        continue
                    # Check if atom a in image i is closer to atom b's position
                    # in image i-1 than to its own previous position
                    pos_a_now = images[i].positions[idx_a]
                    pos_a_prev = images[i - 1].positions[idx_a]
                    pos_b_prev = images[i - 1].positions[idx_b]

                    d_self = np.linalg.norm(pos_a_now[:2] - pos_a_prev[:2])
                    d_cross = np.linalg.norm(pos_a_now[:2] - pos_b_prev[:2])

                    if d_cross < d_self and d_self > 0.5:
                        issues.append(
                            f"Image {i}: atom {idx_a} appears to swap toward "
                            f"atom {idx_b}'s previous position "
                            f"(d_self={d_self:.2f}, d_cross={d_cross:.2f})")

    # --- Check monotonicity of specified distances ---
    if expected_monotonic_pairs:
        for (a, b) in expected_monotonic_pairs:
            key = f'd_{a}_{b}'
            dists = [row[key] for row in table]

            # Determine direction (increasing or decreasing)
            if dists[-1] > dists[0]:
                direction = 'increasing'
                violations = [(i, dists[i], dists[i - 1])
                              for i in range(1, len(dists))
                              if dists[i] < dists[i - 1] - 0.02]
            else:
                direction = 'decreasing'
                violations = [(i, dists[i], dists[i - 1])
                              for i in range(1, len(dists))
                              if dists[i] > dists[i - 1] + 0.02]

            if violations:
                issues.append(
                    f"d({a},{b}) not monotonically {direction}: "
                    f"violations at images {[v[0] for v in violations]}")

    # --- Print table ---
    if verbose:
        # Header
        cols = [f"{'Img':>3}"]
        for idx in adsorbate_indices:
            elem = images[0].get_chemical_symbols()[idx]
            cols.extend([f"{elem}({idx})_x", f"{elem}({idx})_y", f"{elem}({idx})_z"])
        if expected_monotonic_pairs:
            for (a, b) in expected_monotonic_pairs:
                cols.append(f"d({a},{b})")
        if table[0].get('rel_e') is not None:
            cols.append("rel_E")
        cols.append("Note")

        header = f"{'Img':>3}"
        for idx in adsorbate_indices:
            elem = images[0].get_chemical_symbols()[idx]
            header += f" {elem}({idx})x" + f" {elem}({idx})y" + f" {elem}({idx})z"
        if expected_monotonic_pairs:
            for (a, b) in expected_monotonic_pairs:
                header += f" d({a},{b})"
        if table[0].get('rel_e') is not None:
            header += "   rel_E"
        print(f"\n{header}")
        print("-" * len(header))

        # Find TS
        if table[0].get('rel_e') is not None:
            ts_idx = max(range(n_images), key=lambda i: table[i]['rel_e'])
        else:
            ts_idx = -1

        for i, row in enumerate(table):
            line = f"{i:>3}"
            for idx in adsorbate_indices:
                p = row[f'pos_{idx}']
                line += f" {p[0]:>7.3f} {p[1]:>7.3f} {p[2]:>7.3f}"
            if expected_monotonic_pairs:
                for (a, b) in expected_monotonic_pairs:
                    line += f" {row[f'd_{a}_{b}']:>7.3f}"
            if row.get('rel_e') is not None:
                line += f" {row['rel_e']:>7.4f}"
            note = ""
            if i == 0:
                note = "IS"
            elif i == n_images - 1:
                note = "FS"
            elif i == ts_idx:
                note = "<-TS"
            line += f" {note}"
            print(line)

        print(f"\nMax inter-image jump: {max_jump_found:.3f} A "
              f"(threshold: {max_jump} A)")

    result = {
        'valid': len(issues) == 0,
        'issues': issues,
        'table': table,
        'max_jump': max_jump_found,
    }

    if verbose:
        print(f"\n{'PASS' if result['valid'] else 'FAIL'}: "
              f"{len(issues)} issue(s) found")
        for issue in issues:
            print(f"  - {issue}")

    return result
