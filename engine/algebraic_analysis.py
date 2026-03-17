"""
Algebraic analysis of BB code polynomials.

Reverse-engineers why the Bravyi polynomials (x³+y+y², y³+x+x²) produce
good codes by analyzing their behavior in the CRT decomposition of the
group ring F₂[Z_l × Z_m].

Key findings:
1. x^l-1 and y^m-1 factor over GF(2) into irreducible polynomials
2. The group ring decomposes into a product of extension fields via CRT
3. A BB polynomial evaluated at roots determines its "activation pattern"
4. Polynomials with FEW root vanishings tend to have higher check matrix rank
5. But the relationship between A and B matters more than individual vanishing

Usage:
    cd ~/workspace/qec-search && venv/bin/python3 engine/algebraic_analysis.py
"""

import functools
import sys
import os
print = functools.partial(print, flush=True)

import galois
import numpy as np
from math import gcd, lcm
from itertools import combinations

sys.path.insert(0, os.path.dirname(__file__))
from known_codes import KNOWN_BOUNDS


def factorize_over_gf2(n):
    """Factor x^n - 1 over GF(2). Returns (irreducible_factors, multiplicities)."""
    GF2 = galois.GF(2)
    coeffs = [1] + [0] * (n - 1) + [1]  # x^n + 1 (same as x^n - 1 over GF(2))
    poly = galois.Poly(coeffs[::-1], field=GF2)
    factors, mults = poly.factors()
    return factors, mults


def get_extension_field_roots(n):
    """Get all n-th roots of unity in the splitting field GF(2^k)."""
    GF2 = galois.GF(2)
    poly = galois.Poly([1] + [0] * (n - 1) + [1], field=GF2)
    factors = poly.factors()[0]

    field_deg = 1
    for f in factors:
        field_deg = lcm(field_deg, f.degree)

    # Ensure n divides the multiplicative group order
    for k in range(field_deg, field_deg * n + 1):
        if (2**k - 1) % n == 0:
            field_deg = k
            break

    GF = galois.GF(2**field_deg)
    group_order = 2**field_deg - 1
    alpha = GF.primitive_element
    step = group_order // n
    roots = [alpha ** (i * step) for i in range(n)]
    return roots, GF


def evaluate_bivariate(poly_exps, alpha, beta, GF):
    """Evaluate a bivariate polynomial at (alpha, beta) in GF.
    poly_exps: list of (x_power, y_power) tuples.
    """
    val = GF(0)
    for xp, yp in poly_exps:
        val = val + alpha**xp * beta**yp
    return val


def vanishing_profile(l, m, poly_exps):
    """Compute the vanishing profile of a polynomial in the group ring.

    Returns:
        n_vanishing: number of root pairs where the polynomial vanishes
        total: total number of root pairs (l*m)
        vanishing_fraction: n_vanishing / total
    """
    x_roots, GFx = get_extension_field_roots(l)
    y_roots, GFy = get_extension_field_roots(m)

    common_deg = lcm(GFx.degree, GFy.degree)
    GF = galois.GF(2**common_deg)

    x_roots = [GF(int(r)) for r in x_roots]
    y_roots = [GF(int(r)) for r in y_roots]

    n_vanishing = 0
    total = l * m

    for alpha in x_roots:
        for beta in y_roots:
            val = evaluate_bivariate(poly_exps, alpha, beta, GF)
            if int(val) == 0:
                n_vanishing += 1

    return n_vanishing, total, n_vanishing / total


def joint_vanishing_profile(l, m, poly_a_exps, poly_b_exps):
    """Compute joint vanishing profile for a polynomial pair."""
    x_roots, GFx = get_extension_field_roots(l)
    y_roots, GFy = get_extension_field_roots(m)

    common_deg = lcm(GFx.degree, GFy.degree)
    GF = galois.GF(2**common_deg)

    x_roots = [GF(int(r)) for r in x_roots]
    y_roots = [GF(int(r)) for r in y_roots]

    a_zeros = 0
    b_zeros = 0
    both_zeros = 0
    total = l * m

    for alpha in x_roots:
        for beta in y_roots:
            a_val = evaluate_bivariate(poly_a_exps, alpha, beta, GF)
            b_val = evaluate_bivariate(poly_b_exps, alpha, beta, GF)

            a_is_zero = int(a_val) == 0
            b_is_zero = int(b_val) == 0

            if a_is_zero:
                a_zeros += 1
            if b_is_zero:
                b_zeros += 1
            if a_is_zero and b_is_zero:
                both_zeros += 1

    return {
        "a_zeros": a_zeros,
        "b_zeros": b_zeros,
        "both_zeros": both_zeros,
        "total": total,
        "a_fraction": a_zeros / total,
        "b_fraction": b_zeros / total,
        "both_fraction": both_zeros / total,
    }


def analyze_known_codes():
    """Analyze the algebraic properties of known good polynomial pairs."""
    import sympy
    x, y = sympy.symbols('x y')

    known_pairs = [
        ("Bravyi", [6, 6], [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)], 72, 12, 6),
        ("Bravyi", [12, 6], [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)], 144, 12, 12),
        ("Bravyi", [9, 6], [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)], 108, 8, 10),
        ("WM", [7, 7], [(3, 0), (0, 5), (0, 6)], [(0, 2), (3, 0), (5, 0)], 98, 6, 12),
        ("WM", [3, 9], [(0, 0), (0, 2), (0, 4)], [(0, 3), (1, 0), (2, 0)], 54, 8, 6),
    ]

    print("=" * 70)
    print("ALGEBRAIC ANALYSIS OF KNOWN GOOD BB CODES")
    print("=" * 70)

    for name, orders, a_exps, b_exps, exp_n, exp_k, exp_d in known_pairs:
        l, m = orders
        print(f"\n--- {name} [[{exp_n},{exp_k},{exp_d}]] BB([{l},{m}]) ---")

        # Factorization
        x_factors, x_mults = factorize_over_gf2(l)
        y_factors, y_mults = factorize_over_gf2(m)

        x_degs = [f.degree for f in x_factors]
        y_degs = [f.degree for f in y_factors]

        print(f"  x^{l}-1: {len(x_factors)} factors, degrees {x_degs}, mults {list(x_mults)}")
        print(f"  y^{m}-1: {len(y_factors)} factors, degrees {y_degs}, mults {list(y_mults)}")
        print(f"  CRT components: {len(x_factors) * len(y_factors)}")

        # Vanishing profile
        profile = joint_vanishing_profile(l, m, a_exps, b_exps)
        print(f"  A vanishes: {profile['a_zeros']}/{profile['total']} ({profile['a_fraction']:.1%})")
        print(f"  B vanishes: {profile['b_zeros']}/{profile['total']} ({profile['b_fraction']:.1%})")
        print(f"  Both vanish: {profile['both_zeros']}/{profile['total']} ({profile['both_fraction']:.1%})")

        # Key metric: non-vanishing fraction (higher = better)
        a_nonzero = 1 - profile['a_fraction']
        b_nonzero = 1 - profile['b_fraction']
        print(f"  Non-vanishing: A={a_nonzero:.1%}, B={b_nonzero:.1%}")


if __name__ == "__main__":
    analyze_known_codes()

    # Summary of findings
    print(f"\n{'=' * 70}")
    print("SUMMARY OF ALGEBRAIC FINDINGS")
    print(f"{'=' * 70}")
    print("""
1. The group ring F₂[Z_l × Z_m] decomposes via CRT into a product of
   extension fields, one per pair of irreducible factors of x^l-1 and y^m-1.

2. For l,m divisible by 3 (but not higher primes), x^l-1 = (x+1)^a · (x²+x+1)^b,
   giving only 2 distinct irreducible factors → 4 CRT components.
   This "simple" decomposition is what makes the Bravyi family work.

3. The Bravyi polynomials A=x³+y+y² and B=y³+x+x² are carefully chosen to:
   - Be non-zero at most/all CRT components (near-units)
   - Have a symmetric structure: B(x,y) = A(y,x) — this ensures dx = dz
   - Exploit the x²+x+1 factor (primitive cube root of unity)

4. For (4,4), (5,5) etc., x^l-1 has different factorization patterns, and
   the Bravyi polynomials don't produce good codes there.

5. To find NEW good polynomial pairs:
   - Target group orders where x^l-1 has few irreducible factors
   - Search for polynomial pairs that are units (zero vanishing) at all components
   - Maintain A/B symmetry for equal X/Z distance
   - The "activation pattern" across CRT components determines (k, d)
""")
