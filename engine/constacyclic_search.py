"""
Constacyclic code search over GF(4) for quantum error-correcting codes.

Constacyclic codes use x^n - alpha (alpha in GF(4), alpha != 0,1) as the
modulus instead of x^n - 1 (cyclic codes). This gives different polynomial
factorizations and potentially different codes.

For GF(4) = {0, 1, omega, omega^2} with omega^2 + omega + 1 = 0:
  - alpha = omega  gives omega-constacyclic codes
  - alpha = omega^2 gives omega^2-constacyclic codes

The Hermitian dual of an alpha-constacyclic code is alpha^2-constacyclic
(since Frobenius conjugation maps alpha -> alpha^2).

For Hermitian self-orthogonality: we need C subset C^{perp}_H.
We enumerate generator polynomials g(x) | (x^n - alpha) and check
whether the resulting GF(4) code is Hermitian self-orthogonal.

Connection to quasi-twisted codes: an l-quasi-twisted code with twist alpha
decomposes into constacyclic constituents over GF(4)[x]/(x^n - alpha).
The project found [[13,1,5]] via quasi-twisted with twist=omega; constacyclic
codes are the algebraic foundation.

Usage:
    cd ~/workspace/qec-search && venv/bin/python3 engine/constacyclic_search.py
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
from gf4_codes import gf4_to_symplectic, is_hermitian_self_orthogonal
from known_codes import KNOWN_BOUNDS

# ---------------------------------------------------------------------------
# GF(4) field setup
# ---------------------------------------------------------------------------
GF4 = galois.GF(4)
GF2 = galois.GF(2)
OMEGA = GF4.primitive_element       # omega (element 2)
OMEGA2 = OMEGA ** 2                 # omega^2 (element 3)

# Timeout for distance computation per code
DIST_TIMEOUT = 30

# Global timeout for entire search
GLOBAL_TIMEOUT = 3600  # 1 hour

# ---------------------------------------------------------------------------
# Polynomial utilities
# ---------------------------------------------------------------------------

def poly_to_row(p, n, field=GF4):
    """Convert a galois.Poly to a length-n coefficient vector.

    Returns array where index i holds the coefficient of x^i.
    """
    row = np.zeros(n, dtype=int)
    # galois Poly.coefficients() returns highest-degree first
    coeffs = p.coefficients()[::-1]  # now coeffs[i] = coeff of x^i
    for i in range(min(len(coeffs), n)):
        row[i] = int(coeffs[i])
    return row


def hermitian_conjugate_poly(p, field=GF4):
    """Compute the Hermitian conjugate (reciprocal + Frobenius) of polynomial p.

    g^dagger(x) = x^deg(g) * conj(g(1/x))
    where conj applies Frobenius (squaring) to each coefficient.

    Equivalently: reverse coefficients and square each.
    """
    coeffs = list(p.coefficients())  # highest degree first
    # Reverse = reciprocal, square each = Frobenius conjugation
    conj_rev = [c ** 2 for c in coeffs]  # already in reversed order for reciprocal
    # Wait: coefficients() is [a_deg, a_{deg-1}, ..., a_0]
    # Reciprocal means [a_0, a_1, ..., a_deg] = reverse
    # Then conjugate each: [conj(a_0), conj(a_1), ..., conj(a_deg)]
    conj_coeffs = [c ** 2 for c in reversed(coeffs)]
    # Back to highest-degree-first for galois.Poly
    result_coeffs = list(reversed(conj_coeffs))
    return galois.Poly(field(result_coeffs), field=field)


def build_constacyclic_matrix(g, n, alpha, field=GF4):
    """Build generator matrix for alpha-constacyclic code with generator g(x).

    The code C = <g(x)> in GF(4)[x]/(x^n - alpha).
    Rows: g(x), x*g(x), x^2*g(x), ..., x^{k-1}*g(x) mod (x^n - alpha)
    where k = n - deg(g).
    """
    modulus = galois.Poly(field([1] + [0] * (n - 1) + [alpha]), field=field)
    deg_g = g.degree
    k_code = n - deg_g

    if k_code < 1:
        return None

    rows = []
    current = g
    for i in range(k_code):
        row = poly_to_row(current, n, field)
        rows.append(row)
        # Multiply by x and reduce mod (x^n - alpha)
        current = (current * galois.Poly([field(1), field(0)], field=field)) % modulus

    return field(np.array(rows))


# ---------------------------------------------------------------------------
# Factor x^n - alpha and enumerate self-orthogonal generators
# ---------------------------------------------------------------------------

def factor_constacyclic_modulus(n, alpha, field=GF4):
    """Factor x^n - alpha over GF(4).

    In characteristic 2, x^n - alpha = x^n + alpha.
    Returns (factors_list, multiplicities_list).
    """
    modulus = galois.Poly(field([1] + [0] * (n - 1) + [alpha]), field=field)
    factors, mults = modulus.factors()
    return factors, mults


def enumerate_generator_polys(factors, mults):
    """Enumerate all generator polynomials from factor decomposition.

    g(x) = product of some subset of factors (with multiplicity).
    Each factor can appear 0 to mult times.

    Returns list of galois.Poly.
    """
    # Build list of (factor, max_multiplicity) pairs
    factor_choices = []
    for f, m in zip(factors, mults):
        factor_choices.append((f, m))

    # Enumerate all combinations of multiplicities
    generators = []
    _enumerate_recursive(factor_choices, 0, galois.Poly([GF4(1)], field=GF4), generators)
    return generators


def _enumerate_recursive(factor_choices, idx, current_g, results):
    """Recursively enumerate generator polynomials."""
    if idx == len(factor_choices):
        results.append(current_g)
        return

    f, max_m = factor_choices[idx]
    for m in range(max_m + 1):
        g = current_g
        for _ in range(m):
            g = g * f
        _enumerate_recursive(factor_choices, idx + 1, g, results)


def find_self_orthogonal_constacyclic(n, alpha, field=GF4):
    """Find all Hermitian self-orthogonal constacyclic codes for given n, alpha.

    Returns list of (g, H_gf4, n_code, k_quantum) tuples where H_gf4 is
    the GF(4) generator matrix that is Hermitian self-orthogonal.
    """
    factors, mults = factor_constacyclic_modulus(n, alpha, field)

    # Enumerate all possible generator polynomials
    generators = enumerate_generator_polys(factors, mults)

    results = []
    for g in generators:
        deg_g = g.degree
        if deg_g == 0 or deg_g >= n:
            continue  # Skip trivial codes (full space or {0})

        k_classical = n - deg_g  # Classical dimension
        if k_classical < 1 or k_classical >= n:
            continue

        # Build generator matrix
        H = build_constacyclic_matrix(g, n, alpha, field)
        if H is None:
            continue

        # Check Hermitian self-orthogonality
        if is_hermitian_self_orthogonal(H):
            results.append((g, H))

    return results


# ---------------------------------------------------------------------------
# Also try combining rows from BOTH alpha=omega and alpha=omega^2 codes,
# and try extending self-orthogonal codes with additional rows
# ---------------------------------------------------------------------------

def try_extend_with_extra_rows(H_base, n, max_extra=3, rng=None):
    """Try to add random Hermitian-orthogonal rows to enlarge the stabilizer set.

    More stabilizer rows → smaller k → potentially higher d.
    """
    if rng is None:
        rng = np.random.default_rng()

    best_H = H_base
    rows = list(H_base)

    for _ in range(max_extra):
        found = False
        for _ in range(500):
            row = GF4(rng.integers(0, 4, size=n))
            if np.all(row == 0):
                continue

            # Check Hermitian inner product with all existing rows AND self
            from gf4_codes import hermitian_inner_product
            hip_self = hermitian_inner_product(row, row)
            if hip_self != GF4(0):
                continue

            ok = True
            for existing in rows:
                if hermitian_inner_product(row, GF4(existing)) != GF4(0):
                    ok = False
                    break
            if ok:
                rows.append(np.array(row, dtype=int))
                found = True
                break

        if not found:
            break

    if len(rows) > H_base.shape[0]:
        return GF4(np.array(rows))
    return H_base


# ---------------------------------------------------------------------------
# Main search
# ---------------------------------------------------------------------------

def evaluate_code(H_gf4, n, label=""):
    """Convert GF(4) matrix to stabilizer code, compute parameters.

    Returns (n, k, d, status_str) or None if invalid.
    """
    symplectic_H = gf4_to_symplectic(H_gf4)
    code = StabilizerCode(H=symplectic_H, n=n)

    valid, msg = code.is_valid()
    if not valid:
        return None

    k = code.k
    if k < 1:
        return None

    d = compute_distance_qldpc(code, timeout_sec=DIST_TIMEOUT)
    if d is None or d < 2:
        return None

    # Check against known bounds
    key = (n, k)
    status = ""
    if key in KNOWN_BOUNDS:
        dl, du = KNOWN_BOUNDS[key]
        if d > dl:
            status = f" *** IMPROVEMENT over known d={dl} ***"
        elif d == dl:
            status = f" MATCHES(d={dl})"
        elif d == dl - 1:
            status = f" close(known={dl})"
        else:
            status = f" (known={dl})"
    else:
        status = " (no known bound)"

    return n, k, d, status


def search_constacyclic(n_min=5, n_max=25):
    """Main constacyclic code search over GF(4).

    For each odd n in [n_min, n_max], factors x^n - omega and x^n - omega^2,
    enumerates self-orthogonal generator polynomials, evaluates codes.
    """
    start_time = time.time()

    alphas = [
        (OMEGA, "omega"),
        (OMEGA2, "omega^2"),
    ]

    all_results = []   # (n, k, d, status, alpha_name, g_str)
    matches = []       # Codes that match known bounds
    improvements = []  # Codes that beat known bounds

    rng = np.random.default_rng(42)

    for n in range(n_min, n_max + 1):
        if gcd(n, 2) != 1:
            continue  # Skip even n (need gcd(n, char) = 1)

        if time.time() - start_time > GLOBAL_TIMEOUT:
            print(f"\nGlobal timeout reached after {time.time()-start_time:.0f}s")
            break

        print(f"\n{'='*60}")
        print(f"n = {n}")
        print(f"{'='*60}")

        n_codes_found = 0

        for alpha, alpha_name in alphas:
            factors, mults = factor_constacyclic_modulus(n, alpha)
            num_factors = len(factors)
            total_generators = 1
            for m in mults:
                total_generators *= (m + 1)

            print(f"\n  alpha = {alpha_name}: x^{n} + {alpha_name}")
            print(f"  Factors ({num_factors}): ", end="")
            for f, m in zip(factors, mults):
                print(f"({f})^{m} ", end="")
            print(f"\n  Generator polynomials to try: {total_generators}")

            # Find self-orthogonal codes
            so_codes = find_self_orthogonal_constacyclic(n, alpha)

            if not so_codes:
                print(f"  No self-orthogonal codes found")
                continue

            print(f"  Self-orthogonal generators: {len(so_codes)}")

            for g, H_gf4 in so_codes:
                # Evaluate the base constacyclic code
                result = evaluate_code(H_gf4, n, label=f"alpha={alpha_name}")
                if result is not None:
                    _, k, d, status = result
                    g_str = str(g)
                    print(f"    [[{n},{k},{d}]] g={g_str}{status}")

                    entry = (n, k, d, status, alpha_name, g_str)
                    all_results.append(entry)
                    n_codes_found += 1

                    if "MATCHES" in status:
                        matches.append(entry)
                    elif "IMPROVEMENT" in status:
                        improvements.append(entry)

                # Also try extending with additional rows
                H_ext = try_extend_with_extra_rows(H_gf4, n, max_extra=4, rng=rng)
                if H_ext.shape[0] > H_gf4.shape[0]:
                    result_ext = evaluate_code(H_ext, n, label=f"alpha={alpha_name}+ext")
                    if result_ext is not None:
                        _, k_ext, d_ext, status_ext = result_ext
                        # Only report if different from base
                        if k_ext != (result[1] if result else -1) or d_ext != (result[2] if result else -1):
                            print(f"    [[{n},{k_ext},{d_ext}]] extended{status_ext}")

                            entry_ext = (n, k_ext, d_ext, status_ext, alpha_name, f"{g}+ext")
                            all_results.append(entry_ext)
                            n_codes_found += 1

                            if "MATCHES" in status_ext:
                                matches.append(entry_ext)
                            elif "IMPROVEMENT" in status_ext:
                                improvements.append(entry_ext)

        if n_codes_found == 0:
            print(f"  No quantum codes found for n={n}")

    # ---------------------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------------------
    elapsed = time.time() - start_time
    print(f"\n{'='*60}")
    print(f"CONSTACYCLIC SEARCH COMPLETE")
    print(f"{'='*60}")
    print(f"Time: {elapsed:.1f}s")
    print(f"Total codes found: {len(all_results)}")

    if improvements:
        print(f"\n*** IMPROVEMENTS ({len(improvements)}) ***")
        for n, k, d, status, alpha_name, g_str in improvements:
            print(f"  [[{n},{k},{d}]] alpha={alpha_name}, g={g_str}")

    if matches:
        print(f"\nMATCHES ({len(matches)}):")
        for n, k, d, status, alpha_name, g_str in matches:
            print(f"  [[{n},{k},{d}]] alpha={alpha_name}, g={g_str}")

    if all_results:
        print(f"\nALL RESULTS:")
        # Group by n
        by_n = {}
        for entry in all_results:
            n = entry[0]
            if n not in by_n:
                by_n[n] = []
            by_n[n].append(entry)

        for n in sorted(by_n.keys()):
            # Show best d for each k at this n
            best_at_n = {}
            for entry in by_n[n]:
                _, k, d, status, alpha_name, g_str = entry
                if k not in best_at_n or d > best_at_n[k][2]:
                    best_at_n[k] = entry
            for k in sorted(best_at_n.keys()):
                _, _, d, status, alpha_name, g_str = best_at_n[k]
                print(f"  [[{n},{k},{d}]] alpha={alpha_name}{status}")
    else:
        print("\nNo codes found.")

    return all_results, matches, improvements


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("CONSTACYCLIC CODE SEARCH OVER GF(4)")
    print("=" * 60)
    print(f"Searching n = 5..25 (odd only), alpha in {{omega, omega^2}}")
    print(f"Distance timeout per code: {DIST_TIMEOUT}s")
    print(f"Global timeout: {GLOBAL_TIMEOUT}s")
    print()

    # Quick validation: verify [[5,1,3]] can be built from constacyclic codes
    print("--- Validation ---")
    from gf4_codes import hermitian_inner_product

    # The [[5,1,3]] perfect code is actually related to constacyclic structure
    # Let's verify our basic construction works
    n_test = 5
    for alpha, name in [(OMEGA, "omega"), (OMEGA2, "omega^2")]:
        so_codes = find_self_orthogonal_constacyclic(n_test, alpha)
        print(f"n={n_test}, alpha={name}: {len(so_codes)} self-orthogonal generators")
        for g, H in so_codes:
            result = evaluate_code(H, n_test)
            if result:
                _, k, d, status = result
                print(f"  [[{n_test},{k},{d}]]{status} g={g}")

    print("\n--- Main Search ---")
    all_results, matches, improvements = search_constacyclic(n_min=5, n_max=25)
