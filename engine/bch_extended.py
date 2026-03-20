"""
Extended BCH code search over GF(4) using larger extension fields.

Lengths:
  n=15: GF(16)  -- 15 = 4^2-1
  n=21: GF(64)  -- 21 | (4^3-1)=63
  n=63: GF(64)  -- 63 = 4^3-1
  n=31: GF(1024) -- 31 | (4^5-1)=1023

Enumerate 4-cyclotomic cosets, find HSO defining sets, build GF(4) codes,
puncture/shorten to target lengths, compare against KNOWN_BOUNDS.

Usage:
    cd ~/workspace/qec-search && venv/bin/python3 engine/bch_extended.py
"""

import numpy as np
import galois
import signal
import time
import os
import sys
import functools
from itertools import combinations

print = functools.partial(print, flush=True)
sys.path.insert(0, os.path.dirname(__file__))

from known_codes import KNOWN_BOUNDS
from gf4_codes import (gf4_to_symplectic, GF4, is_hermitian_self_orthogonal,
                        hermitian_inner_product)
from symplectic import StabilizerCode, compute_distance_qldpc

os.nice(15)

GF2 = galois.GF(2)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def cyclotomic_cosets_q(n, q=4):
    visited = set()
    cosets = []
    for i in range(n):
        if i in visited:
            continue
        coset = set()
        j = i
        while j not in coset:
            coset.add(j)
            j = (q * j) % n
        cosets.append(sorted(coset))
        visited.update(coset)
    return cosets


def check_hso(Z, n):
    comp = set(range(n)) - Z
    return all((-2 * j) % n in Z for j in comp)


def find_extension_degree(n, q=4):
    m, power = 1, q
    while m <= 20:
        if (power - 1) % n == 0:
            return m
        m += 1
        power *= q
    return None


def build_gf4_subfield_map(GFext):
    subfield = [int(x) for x in GFext.elements if x**4 == x]
    assert len(subfield) == 4
    gfext_one = int(GFext(1))
    others = [x for x in subfield if x != 0 and x != gfext_one]
    assert len(others) == 2
    a, b = others
    assert int(GFext(a) + GFext(b)) == gfext_one
    if int(GFext(a)**2) == b:
        ext_to_gf4 = {0: 0, gfext_one: 1, a: 2, b: 3}
    else:
        ext_to_gf4 = {0: 0, gfext_one: 1, a: 3, b: 2}
    gf4_to_ext = {v: k for k, v in ext_to_gf4.items()}
    return ext_to_gf4, gf4_to_ext


def build_decomposition_table(GFext, ext_to_gf4, ext_deg):
    beta = GFext.primitive_element
    basis = [beta**i for i in range(ext_deg)]
    gf4_ext_elts = sorted(ext_to_gf4.keys())
    decomp = {}

    if ext_deg == 2:
        for a in gf4_ext_elts:
            for b in gf4_ext_elts:
                val = int(GFext(a) * basis[0] + GFext(b) * basis[1])
                decomp[val] = (ext_to_gf4[a], ext_to_gf4[b])
    elif ext_deg == 3:
        for a in gf4_ext_elts:
            va = GFext(a) * basis[0]
            for b in gf4_ext_elts:
                vab = va + GFext(b) * basis[1]
                for c in gf4_ext_elts:
                    val = int(vab + GFext(c) * basis[2])
                    decomp[val] = (ext_to_gf4[a], ext_to_gf4[b], ext_to_gf4[c])
    elif ext_deg == 5:
        for a in gf4_ext_elts:
            va = GFext(a) * basis[0]
            for b in gf4_ext_elts:
                vab = va + GFext(b) * basis[1]
                for c in gf4_ext_elts:
                    vabc = vab + GFext(c) * basis[2]
                    for d in gf4_ext_elts:
                        vabcd = vabc + GFext(d) * basis[3]
                        for e in gf4_ext_elts:
                            val = int(vabcd + GFext(e) * basis[4])
                            decomp[val] = (ext_to_gf4[a], ext_to_gf4[b], ext_to_gf4[c],
                                           ext_to_gf4[d], ext_to_gf4[e])
    else:
        raise ValueError(f"Degree {ext_deg} not implemented")

    assert len(decomp) == GFext.order
    return decomp


def build_bch_code(Z_list, n_cyclic, GFext, decomp, ext_deg):
    """Build GF(4) code from BCH defining zeros via extension field."""
    group_order = GFext.order - 1
    alpha = GFext.primitive_element ** (group_order // n_cyclic)
    nZ = len(Z_list)

    total_rows = nZ * ext_deg
    H_gf4_check = np.zeros((total_rows, n_cyclic), dtype=int)

    for zi_idx, z in enumerate(Z_list):
        for j in range(n_cyclic):
            val = int(alpha ** (z * j))
            coeffs = decomp[val]
            for m in range(ext_deg):
                H_gf4_check[zi_idx * ext_deg + m, j] = coeffs[m]

    # GF(2) multiplication constraints
    n_bin_rows = total_rows
    C_mat = np.zeros((2 * n_bin_rows, 2 * n_cyclic), dtype=int)
    for i in range(n_bin_rows):
        for j in range(n_cyclic):
            val = H_gf4_check[i, j]
            h0 = val & 1
            h1 = (val >> 1) & 1
            C_mat[2*i, 2*j] = h0
            C_mat[2*i, 2*j+1] = h1
            C_mat[2*i+1, 2*j] = h1
            C_mat[2*i+1, 2*j+1] = (h0 ^ h1)

    try:
        null = GF2(C_mat).null_space()
    except Exception:
        return None
    if null.shape[0] < 1:
        return None

    gf4_rows = []
    seen = set()
    for ni in range(null.shape[0]):
        v = null[ni]
        row = tuple(int(v[2*j]) + 2*int(v[2*j+1]) for j in range(n_cyclic))
        if row not in seen and any(x != 0 for x in row):
            seen.add(row)
            gf4_rows.append(list(row))
    if not gf4_rows:
        return None
    return GF4(np.array(gf4_rows))


# ---------------------------------------------------------------------------
# Evaluate
# ---------------------------------------------------------------------------

def evaluate_code(H_gf4, n, timeout=20):
    if H_gf4 is None or H_gf4.shape[0] < 1:
        return None
    if not is_hermitian_self_orthogonal(H_gf4):
        return None
    sym = gf4_to_symplectic(H_gf4)
    code = StabilizerCode(H=sym, n=n)
    valid, _ = code.is_valid()
    if not valid or code.k < 1:
        return None
    # Scale timeout with k (more logical qubits = harder distance)
    t = min(timeout, max(5, timeout - code.k))
    d = compute_distance_qldpc(code, timeout_sec=t)
    if d is None or d < 2:
        return None
    return (n, code.k, d)


def status_str(n, k, d):
    key = (n, k)
    if key in KNOWN_BOUNDS:
        dl, du = KNOWN_BOUNDS[key]
        if d > dl:
            return f" *** IMPROVEMENT (known d={dl}) ***"
        elif d == dl:
            return f" MATCHES(d={dl})"
        else:
            return f" (vs d={dl})"
    return ""


all_results = {}


def record(n, k, d, H_gf4, label=""):
    key = (n, k)
    if key not in all_results or d > all_results[key][0]:
        all_results[key] = (d, H_gf4)
        s = status_str(n, k, d)
        print(f"  {label}[[{n},{k},{d}]]{s}")
        return True
    return False


# ---------------------------------------------------------------------------
# HSO defining sets
# ---------------------------------------------------------------------------

def get_hso_defining_sets(n, max_enum=18):
    cosets = cyclotomic_cosets_q(n, q=4)
    print(f"  {len(cosets)} cosets mod {n}: sizes {[len(c) for c in cosets]}")

    nc = min(len(cosets), max_enum)
    valid = []
    for r in range(1, nc + 1):
        for combo in combinations(range(nc), r):
            Z = set()
            for ci in combo:
                Z.update(cosets[ci])
            if check_hso(Z, n):
                dim = n - len(Z)
                if 1 <= dim <= n - 1:
                    valid.append(sorted(Z))
    return valid


# ---------------------------------------------------------------------------
# BCH phase
# ---------------------------------------------------------------------------

def phase_bch(n_cyclic, label, dist_timeout=30, skip_dist_above_k=20):
    print(f"\n{'='*60}")
    print(f"BCH length {n_cyclic} ({label})")
    print(f"{'='*60}")

    ext_deg = find_extension_degree(n_cyclic, q=4)
    if ext_deg is None:
        print(f"  No extension field!")
        return {}

    q_ext = 4 ** ext_deg
    print(f"  GF(4^{ext_deg})=GF({q_ext})")

    GFext = galois.GF(q_ext)
    ext_to_gf4, gf4_to_ext = build_gf4_subfield_map(GFext)

    t0 = time.time()
    decomp = build_decomposition_table(GFext, ext_to_gf4, ext_deg)
    print(f"  Decomp table: {time.time()-t0:.1f}s")

    defs = get_hso_defining_sets(n_cyclic)
    print(f"  {len(defs)} HSO defining sets")

    base_codes = {}
    for Z in defs:
        H = build_bch_code(Z, n_cyclic, GFext, decomp, ext_deg)
        if H is None:
            continue
        if not is_hermitian_self_orthogonal(H):
            continue

        base_codes[tuple(Z)] = H

        sym = gf4_to_symplectic(H)
        code = StabilizerCode(H=sym, n=n_cyclic)
        valid, _ = code.is_valid()
        if not valid or code.k < 1:
            continue

        if code.k > skip_dist_above_k:
            print(f"  BCH-{n_cyclic} |Z|={len(Z)} rows={H.shape[0]} k={code.k} (skip dist, k>{skip_dist_above_k})")
            continue

        d = compute_distance_qldpc(code, timeout_sec=dist_timeout)
        if d is not None and d >= 2:
            record(n_cyclic, code.k, d, H, f"BCH-{n_cyclic} |Z|={len(Z)}: ")
        else:
            print(f"  BCH-{n_cyclic} |Z|={len(Z)} k={code.k} d={'timeout' if d is None else d}")

    return base_codes


# ---------------------------------------------------------------------------
# Puncture
# ---------------------------------------------------------------------------

def phase_puncture(base_codes, n_orig, target_range, n_trials=200, time_limit=120):
    targets = [t for t in target_range if t < n_orig]
    if not targets or not base_codes:
        return
    print(f"\n{'='*60}")
    print(f"PUNCTURE {n_orig} -> [{min(targets)}..{max(targets)}]")
    print(f"{'='*60}")

    rng = np.random.default_rng(42 + n_orig)
    t0 = time.time()

    for Z_key, H_full in base_codes.items():
        if H_full is None:
            continue
        for n_target in targets:
            if time.time() - t0 > time_limit:
                print(f"  Time limit ({time_limit}s) reached")
                return
            for _ in range(n_trials):
                cols = sorted(rng.choice(n_orig, n_target, replace=False))
                H_p = H_full[:, cols]
                mask = np.array([np.any(H_p[i] != GF4(0)) for i in range(H_p.shape[0])])
                if not np.any(mask):
                    continue
                H_p = H_p[mask]
                result = evaluate_code(H_p, n_target, timeout=15)
                if result:
                    _, k, d = result
                    record(n_target, k, d, H_p, f"P{n_orig}->{n_target}: ")


# ---------------------------------------------------------------------------
# Shorten
# ---------------------------------------------------------------------------

def phase_shorten(base_codes, n_orig, target_range, n_trials=150, time_limit=120):
    targets = [t for t in target_range if t < n_orig]
    if not targets or not base_codes:
        return
    print(f"\n{'='*60}")
    print(f"SHORTEN {n_orig} -> [{min(targets)}..{max(targets)}]")
    print(f"{'='*60}")

    rng = np.random.default_rng(123 + n_orig)
    t0 = time.time()

    for Z_key, H_full in base_codes.items():
        if H_full is None or H_full.shape[0] < 2:
            continue
        for n_target in targets:
            if time.time() - t0 > time_limit:
                print(f"  Time limit ({time_limit}s) reached")
                return
            n_fix = n_orig - n_target
            for _ in range(n_trials):
                fix_cols = sorted(rng.choice(n_orig, n_fix, replace=False))
                keep_cols = [j for j in range(n_orig) if j not in fix_cols]
                keep_rows = [i for i in range(H_full.shape[0])
                             if all(H_full[i, c] == GF4(0) for c in fix_cols)]
                if len(keep_rows) < 1:
                    continue
                H_s = H_full[np.array(keep_rows)][:, np.array(keep_cols)]
                result = evaluate_code(H_s, n_target, timeout=15)
                if result:
                    _, k, d = result
                    record(n_target, k, d, H_s, f"S{n_orig}->{n_target}: ")


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def print_summary():
    print(f"\n{'='*60}")
    print("EXTENDED BCH SEARCH RESULTS")
    print(f"{'='*60}")

    matches, improvements, other = [], [], []

    for (n, k), (d, _) in sorted(all_results.items()):
        if (n, k) in KNOWN_BOUNDS:
            dl, du = KNOWN_BOUNDS[(n, k)]
            if d > dl:
                improvements.append((n, k, d, dl, du))
            elif d == dl:
                matches.append((n, k, d, dl, du))
            else:
                other.append((n, k, d, dl, du))
        else:
            other.append((n, k, d, None, None))

    if improvements:
        print("\n*** IMPROVEMENTS OVER KNOWN BOUNDS ***")
        for n, k, d, dl, du in improvements:
            print(f"  [[{n},{k},{d}]] -- was d={dl}, upper={du}")

    if matches:
        print(f"\nMATCHES ({len(matches)}):")
        for n, k, d, dl, du in matches:
            opt = " (optimal)" if dl == du else ""
            print(f"  [[{n},{k},{d}]]{opt}")

    if other:
        print(f"\nOther ({len(other)}):")
        for n, k, d, dl, du in other:
            if dl is not None:
                print(f"  [[{n},{k},{d}]] (known: d={dl})")
            else:
                print(f"  [[{n},{k},{d}]]")

    print(f"\nTotal: {len(all_results)} unique (n,k)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    GLOBAL_TIMEOUT = 1500  # 25 minutes

    def alarm_handler(signum, frame):
        print("\n*** GLOBAL TIMEOUT ***")
        print_summary()
        sys.exit(0)

    signal.signal(signal.SIGALRM, alarm_handler)
    signal.alarm(GLOBAL_TIMEOUT)
    t_start = time.time()

    # ===== n=15: GF(16) =====
    codes_15 = phase_bch(15, "GF(16)", dist_timeout=20)
    phase_puncture(codes_15, 15, range(5, 15), n_trials=150, time_limit=60)
    phase_shorten(codes_15, 15, range(5, 15), n_trials=100, time_limit=60)
    print(f"\n  [{time.time()-t_start:.0f}s elapsed]")

    # ===== n=21: GF(64) =====
    codes_21 = phase_bch(21, "GF(64)", dist_timeout=30, skip_dist_above_k=12)
    phase_puncture(codes_21, 21, range(17, 21), n_trials=300, time_limit=90)
    phase_shorten(codes_21, 21, range(17, 21), n_trials=200, time_limit=90)
    phase_puncture(codes_21, 21, range(5, 17), n_trials=80, time_limit=60)
    phase_shorten(codes_21, 21, range(5, 17), n_trials=60, time_limit=60)
    # Extend 21 -> 22..25: just add parity-like columns
    print(f"\n  [{time.time()-t_start:.0f}s elapsed]")

    # ===== n=63: GF(64) -- big codes, focus on puncture/shorten =====
    codes_63 = phase_bch(63, "GF(64), 63=4^3-1", dist_timeout=20, skip_dist_above_k=10)
    # Puncture/shorten 63 -> 17..31 (focus on interesting range)
    phase_puncture(codes_63, 63, range(17, 32), n_trials=250, time_limit=180)
    phase_shorten(codes_63, 63, range(17, 32), n_trials=200, time_limit=180)
    print(f"\n  [{time.time()-t_start:.0f}s elapsed]")

    # ===== n=31: GF(1024) =====
    print(f"\nAttempting n=31 via GF(1024)...")
    try:
        codes_31 = phase_bch(31, "GF(1024)", dist_timeout=30, skip_dist_above_k=10)
        phase_puncture(codes_31, 31, range(17, 31), n_trials=200, time_limit=120)
        phase_shorten(codes_31, 31, range(17, 31), n_trials=150, time_limit=120)
    except Exception as e:
        print(f"  n=31 failed: {e}")

    signal.alarm(0)
    elapsed = time.time() - t_start
    print(f"\nTotal time: {elapsed:.1f}s")
    print_summary()
