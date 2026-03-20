"""
Rich factorization constacyclic search.

Key insight: x^n - omega factors into MANY small irreducible polynomials over
GF(4) for certain n values. Rich factorizations = more generator combinations
= better codes. This script:

1. Surveys all odd n from 5 to 50, factoring x^n - omega and x^n - omega^2
2. Scores by "richness" = num_factors / max_factor_degree
3. Runs the constacyclic code search on the top-15 most promising n values

Usage:
    cd ~/workspace/qec-search && venv/bin/python3 engine/rich_factor_search.py
"""

import numpy as np
import galois
import signal
import time
import os
import sys
import functools
from itertools import combinations
from math import gcd

# Nice priority
os.nice(15)

# Flush prints
print = functools.partial(print, flush=True)

sys.path.insert(0, os.path.dirname(__file__))
from symplectic import StabilizerCode, compute_distance_qldpc
from gf4_codes import (gf4_to_symplectic, is_hermitian_self_orthogonal,
                        hermitian_inner_product)
from known_codes import KNOWN_BOUNDS
from codes import gf2_rank

# ---------------------------------------------------------------------------
# GF(4) field setup
# ---------------------------------------------------------------------------
GF4 = galois.GF(4)
GF2 = galois.GF(2)
OMEGA = GF4.primitive_element       # omega (element 2)
OMEGA2 = OMEGA ** 2                 # omega^2 (element 3)

# Timeouts
DIST_TIMEOUT = 30       # per code distance computation
PER_N_TIMEOUT = 300     # per block length (generous for promising n)
FACTOR_TIMEOUT = 60     # factorization timeout


# ---------------------------------------------------------------------------
# Timeout helper
# ---------------------------------------------------------------------------
class AlarmTimeout(Exception):
    pass

def _alarm_handler(signum, frame):
    raise AlarmTimeout()


# ---------------------------------------------------------------------------
# Constacyclic shift operator
# ---------------------------------------------------------------------------

def constacyclic_shift(v, alpha):
    """Apply constacyclic shift: (c_0,...,c_{n-1}) -> (alpha*c_{n-1}, c_0,...,c_{n-2})."""
    n = len(v)
    result = np.zeros(n, dtype=int)
    result[0] = int(alpha * GF4(int(v[n - 1])))
    for i in range(1, n):
        result[i] = int(v[i - 1])
    return GF4(result)


def generate_shifts(v, alpha, n):
    """Generate all n constacyclic shifts of v."""
    shifts = [GF4(np.array(v))]
    current = v
    for _ in range(n - 1):
        current = constacyclic_shift(current, alpha)
        shifts.append(current)
    return shifts


# ---------------------------------------------------------------------------
# Polynomial arithmetic over GF(4)
# ---------------------------------------------------------------------------

def poly_to_vec(p, n):
    """Convert polynomial to length-n vector (coeff of x^i at index i)."""
    row = np.zeros(n, dtype=int)
    coeffs = p.coefficients()[::-1]
    for i in range(min(len(coeffs), n)):
        row[i] = int(coeffs[i])
    return row


def build_code_matrix(g, n, alpha):
    """Build generator matrix for alpha-constacyclic code with generator g(x)."""
    modulus_coeffs = [GF4(1)] + [GF4(0)] * (n - 1) + [alpha]
    modulus = galois.Poly(GF4(modulus_coeffs), field=GF4)
    deg_g = g.degree
    k_code = n - deg_g
    if k_code < 1:
        return None

    rows = []
    current = g
    x_poly = galois.Poly([GF4(1), GF4(0)], field=GF4)
    for _ in range(k_code):
        rows.append(poly_to_vec(current, n))
        current = (current * x_poly) % modulus

    return GF4(np.array(rows))


# ---------------------------------------------------------------------------
# Self-orthogonal subcode via constacyclic shifts
# ---------------------------------------------------------------------------

def find_shift_orthogonal_code(seed_vec, alpha, n):
    """Given a seed vector, generate shifts and select maximal mutually
    Hermitian-orthogonal subset."""
    seed = GF4(np.array(seed_vec))

    if hermitian_inner_product(seed, seed) != GF4(0):
        return None

    shifts = generate_shifts(seed, alpha, n)

    selected = [shifts[0]]
    for i in range(1, len(shifts)):
        ok = True
        for s in selected:
            if hermitian_inner_product(shifts[i], s) != GF4(0):
                ok = False
                break
        if ok:
            test_mat = GF4(np.array([list(r) for r in selected + [shifts[i]]]))
            test_sym = gf4_to_symplectic(test_mat)
            prev_mat = GF4(np.array([list(r) for r in selected]))
            prev_sym = gf4_to_symplectic(prev_mat)
            if gf2_rank(test_sym) > gf2_rank(prev_sym):
                selected.append(shifts[i])

    if len(selected) < 1:
        return None

    H = GF4(np.array([list(r) for r in selected]))
    if not is_hermitian_self_orthogonal(H):
        return None
    return H


def find_multi_seed_code(seeds, alpha, n):
    """Combine shifts from multiple seed vectors into a larger stabilizer code."""
    all_shifts = []
    for seed in seeds:
        shifts = generate_shifts(seed, alpha, n)
        all_shifts.extend(shifts)

    extended = list(all_shifts)
    for v in all_shifts:
        extended.append(OMEGA * v)
        extended.append(OMEGA2 * v)

    seen = set()
    unique = []
    for v in extended:
        key = tuple(int(x) for x in v)
        if key not in seen and any(x != 0 for x in key):
            seen.add(key)
            unique.append(v)

    selected = []
    selected_sym = None

    for v in unique:
        if hermitian_inner_product(v, v) != GF4(0):
            continue

        ok = True
        for s in selected:
            if hermitian_inner_product(v, s) != GF4(0):
                ok = False
                break
        if not ok:
            continue

        v_sym = gf4_to_symplectic(GF4(np.array(v).reshape(1, -1)))
        if selected_sym is not None:
            test_sym = np.vstack([selected_sym, v_sym])
            if gf2_rank(test_sym) <= gf2_rank(selected_sym):
                continue
            selected_sym = test_sym
        else:
            selected_sym = v_sym

        selected.append(v)
        if len(selected) >= n - 1:
            break

    if len(selected) < 1:
        return None

    H = GF4(np.array([list(r) for r in selected]))
    if not is_hermitian_self_orthogonal(H):
        return None
    return H


# ---------------------------------------------------------------------------
# Evaluate quantum code
# ---------------------------------------------------------------------------

def evaluate_code(H_gf4, n):
    """Compute [[n,k,d]] and compare against known bounds."""
    sym_H = gf4_to_symplectic(H_gf4)
    code = StabilizerCode(H=sym_H, n=n)

    valid, _ = code.is_valid()
    if not valid:
        return None

    k = code.k
    if k < 1:
        return None

    # Skip codes with high k (distance computation too expensive)
    if n >= 30 and k > 8:
        return None
    if n >= 20 and k > 10:
        return None
    if n >= 15 and k > 12:
        return None

    timeout = min(DIST_TIMEOUT, max(10, 60 - n))

    d = compute_distance_qldpc(code, timeout_sec=timeout)
    if d is None or d < 2:
        return None

    key = (n, k)
    status = ""
    if key in KNOWN_BOUNDS:
        dl, du = KNOWN_BOUNDS[key]
        if d > dl:
            status = f" *** IMPROVEMENT over known d={dl} ***"
        elif d == dl:
            status = f" MATCHES(d={dl})"
        elif d >= dl - 1:
            status = f" close(known={dl})"
        else:
            status = f" (known={dl})"
    else:
        status = " [no bound]"

    return n, k, d, status


# ---------------------------------------------------------------------------
# Phase 1: Factor x^n - omega and x^n - omega^2 for all odd n in [5..50]
# ---------------------------------------------------------------------------

def factor_analysis(n_min=5, n_max=50):
    """Analyze factorizations of x^n - omega and x^n - omega^2 over GF(4)."""
    print("=" * 70)
    print("PHASE 1: FACTORIZATION ANALYSIS (x^n - omega, x^n - omega^2 over GF(4))")
    print("=" * 70)
    print(f"{'n':>4} | {'alpha':>5} | {'#fac':>4} | {'max_deg':>7} | {'richness':>8} | pattern")
    print("-" * 70)

    entries = []  # (richness, n, alpha_val, alpha_name, num_factors, max_deg, pattern)

    for n in range(n_min, n_max + 1):
        if gcd(n, 2) != 1:
            continue

        for alpha, alpha_name in [(OMEGA, "w"), (OMEGA2, "w2")]:
            # Build x^n - alpha = x^n + alpha (char 2)
            modulus_coeffs = [GF4(1)] + [GF4(0)] * (n - 1) + [alpha]
            modulus = galois.Poly(GF4(modulus_coeffs), field=GF4)

            old_handler = signal.signal(signal.SIGALRM, _alarm_handler)
            signal.alarm(FACTOR_TIMEOUT)
            try:
                factors, mults = modulus.factors()
                signal.alarm(0)
            except AlarmTimeout:
                signal.alarm(0)
                signal.signal(signal.SIGALRM, old_handler)
                print(f"{n:>4} | {alpha_name:>5} | TIMEOUT")
                continue
            finally:
                signal.signal(signal.SIGALRM, old_handler)

            degs = []
            for f, m in zip(factors, mults):
                for _ in range(m):
                    degs.append(f.degree)

            num_factors = len(degs)
            max_deg = max(degs) if degs else n
            richness = num_factors / max_deg if max_deg > 0 else 0
            pattern = sorted(degs)

            entries.append((richness, n, alpha, alpha_name, num_factors, max_deg, pattern))
            print(f"{n:>4} | {alpha_name:>5} | {num_factors:>4} | {max_deg:>7} | {richness:>8.2f} | {pattern}")

    # Sort by richness descending
    entries.sort(key=lambda x: -x[0])

    print("\n" + "=" * 70)
    print("TOP 20 BY RICHNESS (num_factors / max_factor_degree):")
    print("=" * 70)
    for i, (richness, n, alpha, alpha_name, num_fac, max_deg, pattern) in enumerate(entries[:20]):
        print(f"  {i+1:>2}. n={n:>3} a={alpha_name:>2}: richness={richness:.2f} "
              f"({num_fac} factors, max_deg={max_deg}) {pattern}")

    return entries


# ---------------------------------------------------------------------------
# Phase 2: Search promising n values
# ---------------------------------------------------------------------------

def search_single_n(n, alphas_to_try, rng=None):
    """Search constacyclic codes for block length n with specified alphas."""
    if rng is None:
        rng = np.random.default_rng(42)

    start = time.time()
    results = []
    best_per_k = {}  # k -> (d, info)

    for alpha, alpha_name in alphas_to_try:
        if time.time() - start > PER_N_TIMEOUT:
            break

        # Factor x^n + alpha
        modulus_coeffs = [GF4(1)] + [GF4(0)] * (n - 1) + [alpha]
        modulus = galois.Poly(GF4(modulus_coeffs), field=GF4)

        old_handler = signal.signal(signal.SIGALRM, _alarm_handler)
        signal.alarm(FACTOR_TIMEOUT)
        try:
            factors, mults = modulus.factors()
            signal.alarm(0)
        except AlarmTimeout:
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old_handler)
            print(f"  a={alpha_name}: factorization timeout (n={n})")
            continue
        finally:
            signal.signal(signal.SIGALRM, old_handler)

        factor_degs = [f.degree for f in factors]
        print(f"  a={alpha_name}: factors={sorted(factor_degs)}", end="")

        nontrivial_degs = [d for d in factor_degs if d > 1]
        if nontrivial_degs and min(nontrivial_degs) > n // 2:
            print(" -- skipping (factors too large)")
            continue

        # Enumerate generator polynomials (products of factor subsets)
        max_combo = min(len(factors), 7)  # slightly more aggressive
        generator_polys = []

        for r in range(1, max_combo + 1):
            for combo in combinations(range(len(factors)), r):
                g = galois.Poly([GF4(1)], field=GF4)
                for idx in combo:
                    for _ in range(mults[idx]):
                        g = g * factors[idx]

                deg_g = g.degree
                if 0 < deg_g < n:
                    generator_polys.append((g, combo))

        print(f", {len(generator_polys)} generators")

        for g, combo in generator_polys:
            if time.time() - start > PER_N_TIMEOUT:
                break

            deg_g = g.degree
            k_classical = n - deg_g

            H_linear = build_code_matrix(g, n, alpha)
            if H_linear is None:
                continue

            # Strategy 1: Use each row as a seed for shift-based construction
            for row_idx in range(H_linear.shape[0]):
                if time.time() - start > PER_N_TIMEOUT:
                    break
                seed = H_linear[row_idx]
                H = find_shift_orthogonal_code(seed, alpha, n)
                if H is None:
                    continue

                result = evaluate_code(H, n)
                if result is not None:
                    _, k, d, status = result
                    if k not in best_per_k or d > best_per_k[k][0]:
                        best_per_k[k] = (d, f"a={alpha_name},g=deg{deg_g},shift")
                        print(f"    [[{n},{k},{d}]] shift(row{row_idx}){status}")
                        results.append((n, k, d, status, alpha_name,
                                        f"shift from g=deg{deg_g}"))

            # Strategy 2: omega/omega^2 scaled seeds
            for row_idx in range(min(H_linear.shape[0], 4)):
                if time.time() - start > PER_N_TIMEOUT:
                    break
                for scalar in [OMEGA, OMEGA2]:
                    seed = scalar * H_linear[row_idx]
                    H = find_shift_orthogonal_code(seed, alpha, n)
                    if H is None:
                        continue

                    result = evaluate_code(H, n)
                    if result is not None:
                        _, k, d, status = result
                        if k not in best_per_k or d > best_per_k[k][0]:
                            best_per_k[k] = (d, f"a={alpha_name},scaled")
                            print(f"    [[{n},{k},{d}]] scaled-shift{status}")
                            results.append((n, k, d, status, alpha_name,
                                            f"scaled shift from g=deg{deg_g}"))

            # Strategy 3: Multi-seed combinations
            if k_classical >= 2:
                for i in range(min(k_classical, 4)):
                    for j in range(i + 1, min(k_classical, 5)):
                        if time.time() - start > PER_N_TIMEOUT:
                            break
                        seeds = [H_linear[i], H_linear[j]]
                        H = find_multi_seed_code(seeds, alpha, n)
                        if H is None:
                            continue
                        result = evaluate_code(H, n)
                        if result is not None:
                            _, k, d, status = result
                            if k not in best_per_k or d > best_per_k[k][0]:
                                best_per_k[k] = (d, f"a={alpha_name},multi")
                                print(f"    [[{n},{k},{d}]] multi-seed{status}")
                                results.append((n, k, d, status, alpha_name,
                                                f"multi-seed from g=deg{deg_g}"))

            # Strategy 4: GF(2) sums of rows as seeds
            for i in range(min(k_classical, 4)):
                for j in range(i + 1, min(k_classical, 5)):
                    if time.time() - start > PER_N_TIMEOUT:
                        break
                    sum_row = H_linear[i] + H_linear[j]
                    if np.all(np.array(sum_row) == 0):
                        continue
                    H = find_shift_orthogonal_code(sum_row, alpha, n)
                    if H is None:
                        continue
                    result = evaluate_code(H, n)
                    if result is not None:
                        _, k, d, status = result
                        if k not in best_per_k or d > best_per_k[k][0]:
                            best_per_k[k] = (d, f"a={alpha_name},sum")
                            print(f"    [[{n},{k},{d}]] sum-seed{status}")
                            results.append((n, k, d, status, alpha_name,
                                            f"sum-seed from g=deg{deg_g}"))

            # Strategy 5: Random GF(4) combinations as seeds
            for _ in range(min(25, k_classical * 3)):
                if time.time() - start > PER_N_TIMEOUT:
                    break
                coeffs = GF4(rng.integers(0, 4, size=k_classical))
                if np.all(np.array(coeffs) == 0):
                    continue
                seed = GF4(np.zeros(n, dtype=int))
                for idx in range(k_classical):
                    seed = seed + coeffs[idx] * H_linear[idx]
                if np.all(np.array(seed) == 0):
                    continue

                H = find_shift_orthogonal_code(seed, alpha, n)
                if H is None:
                    continue
                result = evaluate_code(H, n)
                if result is not None:
                    _, k, d, status = result
                    if k not in best_per_k or d > best_per_k[k][0]:
                        best_per_k[k] = (d, f"a={alpha_name},rand")
                        print(f"    [[{n},{k},{d}]] rand-seed{status}")
                        results.append((n, k, d, status, alpha_name,
                                        f"random seed from g=deg{deg_g}"))

        # Strategy 6: Cross-generator combinations
        if len(generator_polys) >= 4 and time.time() - start < PER_N_TIMEOUT:
            for _ in range(min(40, len(generator_polys))):
                if time.time() - start > PER_N_TIMEOUT:
                    break
                idx1, idx2 = rng.choice(len(generator_polys), 2, replace=False)
                g1, _ = generator_polys[idx1]
                g2, _ = generator_polys[idx2]

                H1 = build_code_matrix(g1, n, alpha)
                H2 = build_code_matrix(g2, n, alpha)
                if H1 is None or H2 is None:
                    continue

                r1 = rng.integers(0, H1.shape[0])
                r2 = rng.integers(0, H2.shape[0])
                seeds = [H1[r1], H2[r2]]
                H = find_multi_seed_code(seeds, alpha, n)
                if H is None:
                    continue
                result = evaluate_code(H, n)
                if result is not None:
                    _, k, d, status = result
                    if k not in best_per_k or d > best_per_k[k][0]:
                        best_per_k[k] = (d, f"a={alpha_name},cross")
                        print(f"    [[{n},{k},{d}]] cross-gen{status}")
                        results.append((n, k, d, status, alpha_name,
                                        f"cross-generator"))

    elapsed = time.time() - start
    print(f"  n={n}: {len(results)} codes ({elapsed:.1f}s)")
    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    overall_start = time.time()

    print("=" * 70)
    print("RICH FACTORIZATION CONSTACYCLIC SEARCH")
    print("=" * 70)
    print(f"Strategy: find n values where x^n - omega factors richly over GF(4),")
    print(f"then search those n values intensively for good quantum codes.")
    print()

    # Phase 1: Factorization analysis
    entries = factor_analysis(n_min=5, n_max=50)

    # Phase 2: Select top 15 unique n values by best richness
    seen_n = set()
    top_n_configs = []  # (n, [(alpha, alpha_name), ...])
    n_to_alphas = {}

    for richness, n, alpha, alpha_name, num_fac, max_deg, pattern in entries:
        if n not in n_to_alphas:
            n_to_alphas[n] = []
        n_to_alphas[n].append((alpha, alpha_name, richness))

    # For each n, get best richness across alphas
    n_best_richness = {}
    for n, alpha_list in n_to_alphas.items():
        n_best_richness[n] = max(r for _, _, r in alpha_list)

    # Sort n by best richness, take top 15
    sorted_n = sorted(n_best_richness.keys(), key=lambda n: -n_best_richness[n])
    top_15 = sorted_n[:15]

    print("\n" + "=" * 70)
    print(f"PHASE 2: CONSTACYCLIC SEARCH ON TOP 15 n VALUES")
    print("=" * 70)
    for n in top_15:
        alphas = n_to_alphas[n]
        alpha_strs = [f"{aname}(r={r:.2f})" for _, aname, r in alphas]
        print(f"  n={n:>3}: richness={n_best_richness[n]:.2f}, alphas={', '.join(alpha_strs)}")
    print()

    all_results = []
    all_matches = []
    all_improvements = []

    rng = np.random.default_rng(42)

    for n in sorted(top_15):
        print(f"\n{'='*60}")
        print(f"n = {n} (richness={n_best_richness[n]:.2f})")
        print(f"{'='*60}")

        # Try all alphas for this n (sorted by richness)
        alphas_for_n = [(alpha, aname) for alpha, aname, _ in
                        sorted(n_to_alphas[n], key=lambda x: -x[2])]

        results = search_single_n(n, alphas_for_n, rng=rng)
        all_results.extend(results)

        for entry in results:
            status = entry[3]
            if "MATCHES" in status:
                all_matches.append(entry)
            elif "IMPROVEMENT" in status:
                all_improvements.append(entry)

    # Final summary
    elapsed = time.time() - overall_start
    print(f"\n{'='*70}")
    print("RICH FACTORIZATION SEARCH COMPLETE")
    print(f"{'='*70}")
    print(f"Time: {elapsed:.1f}s")
    print(f"Total codes found: {len(all_results)}")

    if all_improvements:
        print(f"\n*** IMPROVEMENTS ({len(all_improvements)}) ***")
        for n, k, d, status, alpha_name, info in all_improvements:
            print(f"  [[{n},{k},{d}]] a={alpha_name} {info} {status}")

    if all_matches:
        print(f"\nMATCHES ({len(all_matches)}):")
        for n, k, d, status, alpha_name, info in all_matches:
            print(f"  [[{n},{k},{d}]] a={alpha_name} {info}")

    # Best per (n,k)
    if all_results:
        print(f"\nBEST CODES PER (n,k):")
        best = {}
        for entry in all_results:
            n, k, d, status, alpha_name, info = entry
            key = (n, k)
            if key not in best or d > best[key][2]:
                best[key] = entry
        for key in sorted(best.keys()):
            n, k, d, status, alpha_name, info = best[key]
            print(f"  [[{n},{k},{d}]]{status}")

    print(f"\nDone in {elapsed:.1f}s")
