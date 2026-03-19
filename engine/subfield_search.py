"""
Subfield subcode search for quantum error-correcting codes.

Strategy: Build cyclic/BCH codes over GF(4) using roots of unity in GF(16).
The GF(16) field structure constrains the code's defining set (via cyclotomic cosets),
enabling construction of Hermitian self-orthogonal codes with good distance.

Constructions:
1. Length-15 BCH codes via GF(16) roots (4-cyclotomic cosets mod 15)
2. Punctured versions for lengths 5..14
3. Shortened versions for lengths 5..14
4. Extended codes for lengths 16..20 via self-orthogonal column extension
5. Random augmentation of best codes

Usage:
    cd ~/workspace/qec-search && venv/bin/python3 engine/subfield_search.py
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
                        hermitian_inner_product, hermitian_conjugate_vec)
from symplectic import StabilizerCode, compute_distance_qldpc

# Set nice priority
os.nice(15)

# ---------------------------------------------------------------------------
# Field setup
# ---------------------------------------------------------------------------

GF16 = galois.GF(16)
GF2 = galois.GF(2)

# GF(4) subfield of GF(16): {x : x^4 = x}
GF4_SUBFIELD = frozenset(int(x) for x in GF16.elements if x**4 == x)
GF16_TO_GF4 = {0: 0, 1: 1, 6: 2, 7: 3}
GF4_TO_GF16 = {v: k for k, v in GF16_TO_GF4.items()}

# Decomposition GF(16) = GF(4) + GF(4)*beta
BETA = GF16.primitive_element
DECOMP_TABLE = {}
for a in [0, 1, 6, 7]:
    for b in [0, 1, 6, 7]:
        x = GF16(a) + GF16(b) * BETA
        DECOMP_TABLE[int(x)] = (GF16_TO_GF4[a], GF16_TO_GF4[b])


def gf4_to_gf2_pair(val):
    """GF(4) element -> (x, y) in GF(2)^2 where val = x + y*omega."""
    v = int(val)
    return ((v & 1), ((v >> 1) & 1))  # 0->(0,0), 1->(1,0), 2->(0,1), 3->(1,1)


def compute_Mh(h_int):
    """For h in GF(16): 4x2 GF(2) matrix mapping GF(4) c to decomposition of h*embed(c)."""
    h = GF16(h_int)
    M = np.zeros((4, 2), dtype=int)
    # c = 1 -> embed = GF16(1)
    prod = h * GF16(1)
    a4, b4 = DECOMP_TABLE[int(prod)]
    M[:, 0] = [*gf4_to_gf2_pair(a4), *gf4_to_gf2_pair(b4)]
    # c = omega -> embed = GF16(6)
    prod = h * GF16(6)
    a4, b4 = DECOMP_TABLE[int(prod)]
    M[:, 1] = [*gf4_to_gf2_pair(a4), *gf4_to_gf2_pair(b4)]
    return M


# ---------------------------------------------------------------------------
# Cyclotomic cosets and HSO check
# ---------------------------------------------------------------------------

def cyclotomic_cosets(n, q=4):
    """Compute q-cyclotomic cosets mod n."""
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


def check_hdc(Z, n):
    """Check if defining set Z gives Hermitian dual-containing cyclic code.
    For C^perp_H <= C: if j not in Z, then -2j mod n must be in Z."""
    comp = set(range(n)) - Z
    return all((-2 * j) % n in Z for j in comp)


# ---------------------------------------------------------------------------
# BCH code builder
# ---------------------------------------------------------------------------

def build_bch_code_gf16(Z_list, n_cyclic=15):
    """Build GF(4) code from BCH defining zeros via GF(16) roots."""
    alpha = GF16.primitive_element
    nZ = len(Z_list)

    # Parity check over GF(16): H[i,j] = alpha^(z_i * j)
    H16 = GF16(np.array([[int(alpha ** (z * j)) for j in range(n_cyclic)] for z in Z_list]))

    # Build GF(2) constraint matrix
    C_mat = np.zeros((4 * nZ, 2 * n_cyclic), dtype=int)
    for zi in range(nZ):
        for j in range(n_cyclic):
            Mh = compute_Mh(int(H16[zi, j]))
            C_mat[4 * zi:4 * zi + 4, 2 * j:2 * j + 2] = Mh

    null = GF2(C_mat).null_space()
    if null.shape[0] < 1:
        return None

    # Convert nullspace to GF(4)
    gf4_rows = []
    for ni in range(null.shape[0]):
        v = null[ni]
        gf4_rows.append([int(v[2*j]) + 2*int(v[2*j+1]) for j in range(n_cyclic)])

    return GF4(np.array(gf4_rows))


# ---------------------------------------------------------------------------
# Evaluate quantum code
# ---------------------------------------------------------------------------

def evaluate_code(H_gf4, n, timeout=15):
    """Return (n, k, d) or None."""
    if H_gf4 is None or H_gf4.shape[0] < 1:
        return None
    if not is_hermitian_self_orthogonal(H_gf4):
        return None
    sym = gf4_to_symplectic(H_gf4)
    code = StabilizerCode(H=sym, n=n)
    valid, _ = code.is_valid()
    if not valid or code.k < 1:
        return None
    d = compute_distance_qldpc(code, timeout_sec=timeout)
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


def record(all_results, n, k, d, H_gf4, label=""):
    """Record result if it improves on what we have."""
    key = (n, k)
    if key not in all_results or d > all_results[key][0]:
        all_results[key] = (d, H_gf4)
        s = status_str(n, k, d)
        print(f"  {label}[[{n},{k},{d}]]{s}")
        return True
    return False


# ---------------------------------------------------------------------------
# Get all HSO defining sets for length 15
# ---------------------------------------------------------------------------

def get_hso_defining_sets():
    """Return all valid HSO defining sets for length-15 BCH codes."""
    cosets = cyclotomic_cosets(15, q=4)
    valid = []
    for r in range(1, len(cosets) + 1):
        for combo in combinations(range(len(cosets)), r):
            Z = set()
            for ci in combo:
                Z.update(cosets[ci])
            if check_hdc(Z, 15):
                dim = 15 - len(Z)
                if 1 <= dim <= 14:
                    valid.append(sorted(Z))
    return valid


# Cache the base codes
BASE_CODES = {}  # Z_tuple -> H_gf4


def get_base_code(Z_list):
    key = tuple(Z_list)
    if key not in BASE_CODES:
        BASE_CODES[key] = build_bch_code_gf16(Z_list, n_cyclic=15)
    return BASE_CODES[key]


# ---------------------------------------------------------------------------
# Phase 1: Length 15 BCH
# ---------------------------------------------------------------------------

def phase1_bch15(all_results):
    print("\n" + "=" * 60)
    print("PHASE 1: BCH CODES FROM GF(16), length 15")
    print("=" * 60)

    defs = get_hso_defining_sets()
    print(f"  {len(defs)} HSO defining sets")

    for Z in defs:
        H = get_base_code(Z)
        result = evaluate_code(H, 15)
        if result:
            n, k, d = result
            record(all_results, n, k, d, H, f"BCH Z={Z}: ")


# ---------------------------------------------------------------------------
# Phase 2: Puncturing (lengths 5..14)
# ---------------------------------------------------------------------------

def phase2_puncture(all_results):
    print("\n" + "=" * 60)
    print("PHASE 2: PUNCTURED BCH (lengths 5..14)")
    print("=" * 60)

    rng = np.random.default_rng(42)
    defs = get_hso_defining_sets()

    for Z in defs:
        H_full = get_base_code(Z)
        if H_full is None:
            continue

        for n_target in range(5, 15):
            n_trials = min(200, max(50, int(3 * 15 * 14 / max(1, n_target))))
            for _ in range(n_trials):
                cols = sorted(rng.choice(15, n_target, replace=False))
                H_p = H_full[:, cols]
                # Remove zero rows
                mask = np.array([np.any(H_p[i] != GF4(0)) for i in range(H_p.shape[0])])
                if not np.any(mask):
                    continue
                H_p = H_p[mask]
                result = evaluate_code(H_p, n_target)
                if result:
                    _, k, d = result
                    record(all_results, n_target, k, d, H_p, f"Punct n={n_target}: ")


# ---------------------------------------------------------------------------
# Phase 3: Shortening (lengths 5..14)
# ---------------------------------------------------------------------------

def phase3_shorten(all_results):
    print("\n" + "=" * 60)
    print("PHASE 3: SHORTENED BCH (lengths 5..14)")
    print("=" * 60)

    rng = np.random.default_rng(123)
    defs = get_hso_defining_sets()

    for Z in defs:
        H_full = get_base_code(Z)
        if H_full is None or H_full.shape[0] < 2:
            continue

        for n_target in range(5, 15):
            n_fix = 15 - n_target
            for _ in range(80):
                fix_cols = sorted(rng.choice(15, n_fix, replace=False))
                keep_cols = [j for j in range(15) if j not in fix_cols]

                keep_rows = [i for i in range(H_full.shape[0])
                             if all(H_full[i, c] == GF4(0) for c in fix_cols)]
                if len(keep_rows) < 1:
                    continue

                H_s = H_full[np.array(keep_rows)][:, np.array(keep_cols)]
                result = evaluate_code(H_s, n_target)
                if result:
                    _, k, d = result
                    record(all_results, n_target, k, d, H_s, f"Short n={n_target}: ")


# ---------------------------------------------------------------------------
# Phase 4: Extension (lengths 16..20)
# ---------------------------------------------------------------------------

def phase4_extend(all_results):
    print("\n" + "=" * 60)
    print("PHASE 4: EXTENDED CODES (lengths 16..20)")
    print("=" * 60)

    rng = np.random.default_rng(77)
    defs = get_hso_defining_sets()

    for Z in defs:
        H_full = get_base_code(Z)
        if H_full is None or H_full.shape[0] < 1:
            continue

        rows, n_orig = H_full.shape

        for n_ext in range(16, 21):
            n_add = n_ext - n_orig

            # Strategy 1: Random columns with HSO preservation
            for _ in range(100):
                extra = GF4(rng.integers(0, 4, size=(rows, n_add)))
                H_ext = GF4(np.hstack([np.array(H_full, dtype=int),
                                        np.array(extra, dtype=int)]))
                result = evaluate_code(H_ext, n_ext)
                if result:
                    _, k, d = result
                    record(all_results, n_ext, k, d, H_ext, f"Ext n={n_ext}: ")

            # Strategy 2: Overall parity column(s)
            if n_add >= 1:
                # Add an overall parity column: p_j = sum of other entries
                parity_col = GF4(np.zeros((rows, 1), dtype=int))
                for i in range(rows):
                    s = GF4(0)
                    for j in range(n_orig):
                        s = s + H_full[i, j]
                    parity_col[i, 0] = s
                if n_add == 1:
                    H_ext2 = GF4(np.hstack([np.array(H_full, dtype=int),
                                             np.array(parity_col, dtype=int)]))
                else:
                    extra2 = GF4(rng.integers(0, 4, size=(rows, n_add - 1)))
                    H_ext2 = GF4(np.hstack([np.array(H_full, dtype=int),
                                             np.array(parity_col, dtype=int),
                                             np.array(extra2, dtype=int)]))
                result = evaluate_code(H_ext2, n_ext)
                if result:
                    _, k, d = result
                    record(all_results, n_ext, k, d, H_ext2, f"ExtParity n={n_ext}: ")

            # Strategy 3: Conjugate parity column
            if n_add >= 1:
                conj_col = GF4(np.zeros((rows, 1), dtype=int))
                for i in range(rows):
                    s = GF4(0)
                    for j in range(n_orig):
                        s = s + H_full[i, j] ** 2  # Hermitian conjugate
                    conj_col[i, 0] = s
                if n_add == 1:
                    H_ext3 = GF4(np.hstack([np.array(H_full, dtype=int),
                                             np.array(conj_col, dtype=int)]))
                else:
                    extra3 = GF4(rng.integers(0, 4, size=(rows, n_add - 1)))
                    H_ext3 = GF4(np.hstack([np.array(H_full, dtype=int),
                                             np.array(conj_col, dtype=int),
                                             np.array(extra3, dtype=int)]))
                result = evaluate_code(H_ext3, n_ext)
                if result:
                    _, k, d = result
                    record(all_results, n_ext, k, d, H_ext3, f"ExtConj n={n_ext}: ")


# ---------------------------------------------------------------------------
# Phase 5: Augment best codes with random rows
# ---------------------------------------------------------------------------

def phase5_augment(all_results):
    print("\n" + "=" * 60)
    print("PHASE 5: AUGMENTATION")
    print("=" * 60)

    rng = np.random.default_rng(999)

    for (n, k_orig), (d_orig, H_orig) in list(all_results.items()):
        if n < 5 or n > 20:
            continue

        rows = list(H_orig)
        for _ in range(500):
            new_row = GF4(rng.integers(0, 4, size=n))
            if np.all(new_row == 0):
                continue
            if hermitian_inner_product(new_row, new_row) != GF4(0):
                continue
            ok = True
            for ex in rows:
                if hermitian_inner_product(new_row, ex) != GF4(0):
                    ok = False
                    break
            if ok:
                rows.append(new_row)

        if len(rows) <= H_orig.shape[0]:
            continue

        H_aug = GF4(np.array([list(r) for r in rows]))
        result = evaluate_code(H_aug, n)
        if result:
            _, k, d = result
            record(all_results, n, k, d, H_aug, f"Aug [[{n},{k_orig}]]->{k}: ")

        # Also try various subsets of rows (drop some to get higher k)
        for n_drop in range(1, min(4, len(rows))):
            for _ in range(20):
                drop_idx = rng.choice(len(rows), n_drop, replace=False)
                keep = [rows[i] for i in range(len(rows)) if i not in drop_idx]
                if len(keep) < 1:
                    continue
                H_sub = GF4(np.array([list(r) for r in keep]))
                result = evaluate_code(H_sub, n)
                if result:
                    _, k, d = result
                    record(all_results, n, k, d, H_sub, f"AugSub n={n}: ")


# ---------------------------------------------------------------------------
# Phase 6: Cross-combine: use extension + puncture combinations
# ---------------------------------------------------------------------------

def phase6_cross(all_results):
    print("\n" + "=" * 60)
    print("PHASE 6: CROSS-COMBINATIONS")
    print("=" * 60)

    rng = np.random.default_rng(555)

    # Take best codes for each (n,k), try modifying them
    for (n, k), (d, H) in list(all_results.items()):
        if n < 5 or n > 20:
            continue

        # Column permutation search
        for _ in range(30):
            perm = rng.permutation(n)
            H_perm = H[:, perm]
            # Apply column conjugation (omega <-> omega^2) to some columns
            conj_cols = rng.choice(n, size=rng.integers(0, n+1), replace=False)
            H_mod = H_perm.copy()
            for c in conj_cols:
                H_mod[:, c] = H_mod[:, c] ** 2
            result = evaluate_code(H_mod, n)
            if result:
                _, k2, d2 = result
                record(all_results, n, k2, d2, H_mod, f"ColConj n={n}: ")

        # Row swap + combine from different base codes
        # Try adding a random self-orth row
        for _ in range(100):
            new_row = GF4(rng.integers(0, 4, size=n))
            if np.all(new_row == 0):
                continue
            if hermitian_inner_product(new_row, new_row) != GF4(0):
                continue
            ok = True
            for i in range(H.shape[0]):
                if hermitian_inner_product(new_row, H[i]) != GF4(0):
                    ok = False
                    break
            if not ok:
                continue

            H_new = GF4(np.vstack([np.array(H, dtype=int), [list(new_row)]]))
            result = evaluate_code(H_new, n)
            if result:
                _, k2, d2 = result
                record(all_results, n, k2, d2, H_new, f"AddRow n={n}: ")
                break


# ---------------------------------------------------------------------------
# Phase 7: Random GF(4) self-orthogonal codes for gaps
# ---------------------------------------------------------------------------

def phase7_random_gf4(all_results):
    """Random search for GF(4) Hermitian self-orthogonal codes.

    Targets specific (n,k) pairs from KNOWN_BOUNDS where we haven't matched yet.
    """
    print("\n" + "=" * 60)
    print("PHASE 7: RANDOM GF(4) SEARCH (gap filling)")
    print("=" * 60)

    from gf4_codes import random_gf4_code

    rng = np.random.default_rng(2026)

    # Identify gaps: known bounds we haven't matched
    gaps = []
    for (n, k), (dl, du) in sorted(KNOWN_BOUNDS.items()):
        if 5 <= n <= 20:
            current = all_results.get((n, k), (0, None))[0]
            if current < dl:
                gaps.append((n, k, dl, du))

    if not gaps:
        print("  No gaps to fill!")
        return

    print(f"  {len(gaps)} gaps to fill")

    for n, k_target, dl, du in gaps:
        t0 = time.time()
        # How many stabilizers for this k? n-k rows.
        n_stab = n - k_target

        best_d = all_results.get((n, k_target), (0, None))[0]
        found = False

        # Try many random codes with target number of rows
        n_trials = 300 if n <= 12 else 150 if n <= 16 else 80
        for trial in range(n_trials):
            if time.time() - t0 > 30:  # 30s per gap
                break

            # Vary rows around target
            target_rows = n_stab + rng.integers(-1, 2)
            target_rows = max(1, min(target_rows, n - 1))

            H = random_gf4_code(n, target_rows=target_rows, rng=rng)
            if H is None:
                continue

            sym = gf4_to_symplectic(H)
            code = StabilizerCode(H=sym, n=n)
            valid, _ = code.is_valid()
            if not valid or code.k < 1:
                continue

            d = compute_distance_qldpc(code, timeout_sec=10)
            if d is None or d < 2:
                continue

            key = (n, code.k)
            cur = all_results.get(key, (0, None))[0]
            if d > cur:
                all_results[key] = (d, H)
                s = status_str(n, code.k, d)
                if 'MATCHES' in s or 'IMPROVEMENT' in s:
                    print(f"  Random [[{n},{code.k},{d}]]{s}")
                    found = True

        if not found and time.time() - t0 < 5:
            # Quick evolutionary refinement
            from gf4_codes import mutate_gf4
            # Seed from best code at this n
            seeds = [(d_val, H_val) for (nn, kk), (d_val, H_val) in all_results.items()
                     if nn == n and d_val >= 2]
            if seeds:
                seeds.sort(key=lambda x: -x[0])
                _, seed_H = seeds[0]
                for _ in range(50):
                    if time.time() - t0 > 20:
                        break
                    child = seed_H.copy()
                    for _ in range(rng.integers(1, 5)):
                        child = mutate_gf4(child, n, rng)
                    if child.shape[0] < 1:
                        continue
                    if not is_hermitian_self_orthogonal(child):
                        continue
                    sym = gf4_to_symplectic(child)
                    code = StabilizerCode(H=sym, n=n)
                    valid, _ = code.is_valid()
                    if not valid or code.k < 1:
                        continue
                    d = compute_distance_qldpc(code, timeout_sec=10)
                    if d is None or d < 2:
                        continue
                    key = (n, code.k)
                    cur = all_results.get(key, (0, None))[0]
                    if d > cur:
                        all_results[key] = (d, child)
                        s = status_str(n, code.k, d)
                        if 'MATCHES' in s or 'IMPROVEMENT' in s:
                            print(f"  Evolve [[{n},{code.k},{d}]]{s}")


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def print_summary(all_results):
    print("\n" + "=" * 60)
    print("SUBFIELD SUBCODE SEARCH RESULTS")
    print("=" * 60)

    matches = []
    improvements = []
    other = []

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
        print(f"\nMATCHES ({len(matches)} codes match known best):")
        for n, k, d, dl, du in matches:
            opt = " (optimal)" if dl == du else ""
            print(f"  [[{n},{k},{d}]]{opt}")

    if other:
        print(f"\nOther codes found ({len(other)}):")
        for n, k, d, dl, du in other:
            if dl is not None:
                print(f"  [[{n},{k},{d}]] (known: d={dl})")
            else:
                print(f"  [[{n},{k},{d}]] (no known bound)")

    total = len(all_results)
    print(f"\nTotal unique (n,k) pairs: {total}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    GLOBAL_TIMEOUT = 1200  # 20 minutes

    all_results = {}

    def alarm_handler(signum, frame):
        print("\n*** GLOBAL TIMEOUT ***")
        print_summary(all_results)
        sys.exit(0)

    signal.signal(signal.SIGALRM, alarm_handler)
    signal.alarm(GLOBAL_TIMEOUT)
    t_start = time.time()

    phase1_bch15(all_results)
    phase2_puncture(all_results)
    phase3_shorten(all_results)
    phase4_extend(all_results)
    phase5_augment(all_results)
    phase6_cross(all_results)
    phase7_random_gf4(all_results)

    signal.alarm(0)
    elapsed = time.time() - t_start
    print(f"\nTotal time: {elapsed:.1f}s")
    print_summary(all_results)
