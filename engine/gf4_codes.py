"""
Quantum code construction via additive self-orthogonal codes over GF(4).

This is how the best known codes in codetables.de are actually built.
A qubit stabilizer code [[n,k,d]] ↔ an additive code C ⊆ GF(4)^n such that
C ⊆ C^⊥_H (Hermitian self-orthogonal), where:
  - GF(4) = {0, 1, ω, ω²}, ω² + ω + 1 = 0
  - Pauli correspondence: I↔0, X↔1, Z↔ω, Y↔ω²
  - Hermitian conjugation: x̄ = x² (Frobenius automorphism)
  - Hermitian inner product: ⟨u,v⟩ = Σ uᵢ·v̄ᵢ
  - k = n - rank_GF2(C) (rank as additive/GF(2)-linear code)

Construction methods:
1. From classical GF(4)-linear codes: if C contains C^⊥_H then C gives quantum code
2. BCH-type: use cyclotomic cosets over GF(4) to build cyclic self-orthogonal codes
3. Quadratic residue codes over GF(4)
4. Random search with self-orthogonality constraint

Usage:
    cd ~/workspace/qec-search && venv/bin/python3 engine/gf4_codes.py
"""

import numpy as np
import galois
import signal
import time
import functools
import json
import os
import sys

print = functools.partial(print, flush=True)
sys.path.insert(0, os.path.dirname(__file__))

from known_codes import KNOWN_BOUNDS

# GF(4) field setup
GF4 = galois.GF(4)
GF2 = galois.GF(2)
OMEGA = GF4.primitive_element       # ω (element 2 in galois representation)
OMEGA2 = OMEGA**2                   # ω² = ω + 1 (element 3)

# Element mapping: int representation → Pauli
# galois GF(4): 0=0, 1=1, 2=ω, 3=ω²
PAULI_MAP = {0: 'I', 1: 'X', 2: 'Z', 3: 'Y'}


# ---------------------------------------------------------------------------
# GF(4) self-orthogonality
# ---------------------------------------------------------------------------

def hermitian_conjugate_vec(v):
    """Hermitian conjugation: x → x² for each element."""
    return v**2


def hermitian_inner_product(u, v):
    """Hermitian inner product: ⟨u,v⟩ = Σ uᵢ · v̄ᵢ."""
    return np.sum(u * hermitian_conjugate_vec(v))


def is_hermitian_self_orthogonal(H):
    """Check if the rows of H are mutually Hermitian-orthogonal."""
    H_conj = H**2
    product = H @ H_conj.T
    return np.all(product == GF4(0))


def gf4_code_params(H):
    """Compute [[n, k]] from a GF(4) generator matrix.

    n = number of columns
    k = n - rank_over_GF2(H) where we treat H as a GF(2)-linear code
    (each GF(4) element maps to 2 bits)
    """
    n = H.shape[1]
    # Convert to binary for GF(2) rank
    binary = gf4_to_binary(H)
    rank_gf2 = int(np.linalg.matrix_rank(GF2(binary)))
    k = n - rank_gf2
    return n, k


def gf4_to_binary(H):
    """Convert GF(4) matrix to binary matrix (each element → 2 bits).

    GF(4) element a + bω maps to (a, b) in GF(2)².
    0 → (0,0), 1 → (1,0), ω → (0,1), ω² → (1,1)
    """
    rows, cols = H.shape
    binary = np.zeros((rows, 2 * cols), dtype=int)
    for i in range(rows):
        for j in range(cols):
            val = int(H[i, j])
            if val == 0:
                binary[i, 2*j] = 0; binary[i, 2*j+1] = 0
            elif val == 1:
                binary[i, 2*j] = 1; binary[i, 2*j+1] = 0
            elif val == 2:  # ω
                binary[i, 2*j] = 0; binary[i, 2*j+1] = 1
            elif val == 3:  # ω²
                binary[i, 2*j] = 1; binary[i, 2*j+1] = 1
    return binary


def gf4_to_symplectic(H):
    """Convert GF(4) matrix to symplectic binary [X|Z] form.

    I(0)→[0,0], X(1)→[1,0], Z(ω)→[0,1], Y(ω²)→[1,1]
    """
    rows, n = H.shape
    Hx = np.zeros((rows, n), dtype=np.uint8)
    Hz = np.zeros((rows, n), dtype=np.uint8)
    for i in range(rows):
        for j in range(n):
            val = int(H[i, j])
            if val == 1:    # X
                Hx[i, j] = 1
            elif val == 2:  # Z (ω)
                Hz[i, j] = 1
            elif val == 3:  # Y (ω²)
                Hx[i, j] = 1
                Hz[i, j] = 1
    return np.hstack([Hx, Hz])


# ---------------------------------------------------------------------------
# Construction 1: Random Hermitian self-orthogonal codes
# ---------------------------------------------------------------------------

def random_gf4_code(n, target_rows=None, rng=None):
    """Build a random Hermitian self-orthogonal code over GF(4).

    Strategy: add rows one at a time, each constrained to be
    Hermitian-orthogonal to all previous rows AND self-orthogonal.
    """
    if rng is None:
        rng = np.random.default_rng()

    if target_rows is None:
        target_rows = n - 1  # Maximize stabilizers

    rows = []
    for _ in range(target_rows):
        for _ in range(500):
            row = GF4(rng.integers(0, 4, size=n))
            if np.all(row == 0):
                continue

            # Self-orthogonal: ⟨row, row⟩ = 0
            if hermitian_inner_product(row, row) != GF4(0):
                continue

            # Orthogonal to all existing rows
            ok = True
            for existing in rows:
                if hermitian_inner_product(row, existing) != GF4(0):
                    ok = False
                    break

            if ok:
                rows.append(row)
                break

    if len(rows) < 1:
        return None

    return GF4(np.array([list(r) for r in rows]))


# ---------------------------------------------------------------------------
# Construction 2: Cyclic codes over GF(4) — use cyclotomic cosets
# ---------------------------------------------------------------------------

def gf4_cyclotomic_cosets(n):
    """Compute the 2-cyclotomic cosets mod n for GF(4) codes.

    For GF(4) = GF(2²), the Frobenius automorphism has order 2,
    so cosets are {i, 2i mod n} (or just {i} if 2i ≡ i mod n).
    """
    visited = set()
    cosets = []
    for i in range(n):
        if i in visited:
            continue
        coset = set()
        j = i
        while j not in coset:
            coset.add(j)
            j = (2 * j) % n
        cosets.append(sorted(coset))
        visited.update(coset)
    return cosets


def gf4_cyclic_code(n, defining_set_cosets):
    """Construct a cyclic code over GF(4) from defining set cosets.

    The code is the set of all polynomials c(x) ∈ GF(4)[x]/(x^n - 1)
    such that c(α^i) = 0 for all i in the defining set, where α is
    a primitive n-th root of unity in some extension of GF(4).

    For the code to be Hermitian self-orthogonal, the defining set T
    must satisfy: T ∩ (-2T mod n) = ∅ (no overlap with -2 times T).
    """
    # Find the splitting field for x^n - 1 over GF(4)
    # The order of 4 mod n gives the extension degree
    if n <= 1:
        return None

    from math import gcd
    if gcd(n, 2) != 1:
        return None  # Need gcd(n, char) = 1 for cyclic codes

    # Find multiplicative order of 4 mod n
    ord4 = 1
    power = 4
    while power % n != 1:
        ord4 += 1
        power = (power * 4) % n
        if ord4 > n:
            return None

    # Extension field GF(4^ord4)
    GFext = galois.GF(4**ord4)

    # Find primitive n-th root of unity
    group_order = 4**ord4 - 1
    if group_order % n != 0:
        return None

    alpha = GFext.primitive_element ** (group_order // n)

    # Build the defining set (union of selected cosets)
    defining_set = set()
    for coset in defining_set_cosets:
        defining_set.update(coset)

    # Generator polynomial: product of minimal polynomials for the defining set
    # Instead, directly build the parity check matrix
    # H[i, j] = α^(i*j) for i in defining_set, j in 0..n-1

    H_ext = []
    for i in sorted(defining_set):
        row = [alpha**(i * j) for j in range(n)]
        H_ext.append(row)

    # Project H back to GF(4) if possible
    # The rows should be expressible over GF(4) (due to cyclotomic coset closure)
    # For now, extract the GF(4) subfield representation

    # This is complex — simplified approach: build the binary check matrix
    # and check self-orthogonality

    # Alternative: use the trace construction
    # For each coset leader i, the minimal polynomial m_i(x) has coefficients in GF(4)
    # if the coset has size ≤ 2

    return None  # TODO: implement full cyclic construction


# ---------------------------------------------------------------------------
# Construction 3: Self-orthogonal from known classical GF(4) codes
# ---------------------------------------------------------------------------

def hermitian_dual_contained(H):
    """Check if a GF(4)-linear code C (with generator H) satisfies C^⊥_H ⊆ C.

    This is equivalent to H being Hermitian self-orthogonal.
    """
    return is_hermitian_self_orthogonal(H)


def extend_to_self_orthogonal(H):
    """Given a GF(4) matrix, find additional rows that are Hermitian-orthogonal
    to all existing rows, enlarging the self-orthogonal code."""
    n = H.shape[1]
    rows = list(H)

    # The Hermitian dual has dimension n - rank(H) over GF(4)
    # We want to find rows in the dual that are also self-orthogonal

    rng = np.random.default_rng()
    for _ in range(n * 2):
        row = GF4(rng.integers(0, 4, size=n))
        if np.all(row == 0):
            continue
        if hermitian_inner_product(row, row) != GF4(0):
            continue
        ok = True
        for existing in rows:
            if hermitian_inner_product(row, existing) != GF4(0):
                ok = False
                break
        if ok:
            rows.append(row)

    return GF4(np.array([list(r) for r in rows]))


# ---------------------------------------------------------------------------
# Construction 4: Weight-constrained search (targeting specific d)
# ---------------------------------------------------------------------------

def search_gf4_codes(n, max_codes=1000, rng=None):
    """Generate many random GF(4) self-orthogonal codes at block length n.

    Returns list of (H_gf4, n, k, d) tuples.
    """
    if rng is None:
        rng = np.random.default_rng()

    results = []

    for _ in range(max_codes):
        # Vary the number of generators
        target_rows = rng.integers(1, n)
        H = random_gf4_code(n, target_rows=target_rows, rng=rng)
        if H is None:
            continue

        if not is_hermitian_self_orthogonal(H):
            continue

        # Convert to symplectic and compute distance
        symplectic_H = gf4_to_symplectic(H)

        from symplectic import StabilizerCode, compute_distance_qldpc
        code = StabilizerCode(H=symplectic_H, n=n)
        valid, _ = code.is_valid()
        if not valid or code.k < 1:
            continue

        d = compute_distance_qldpc(code, timeout_sec=15)
        if d is not None and d >= 2:
            results.append((H, n, code.k, d))

    return results


# ---------------------------------------------------------------------------
# Evolutionary GF(4) search
# ---------------------------------------------------------------------------

def mutate_gf4(H, n, rng):
    """Mutate a GF(4) self-orthogonal matrix.

    Operations (all preserve Hermitian self-orthogonality):
    1. Scalar multiply a row by a nonzero GF(4) element
    2. Add a multiple of one row to another
    3. Replace a row with a random Hermitian-orthogonal vector
    4. Permute columns (qubit relabeling)
    5. Apply field automorphism (conjugation) to a column
    """
    H = H.copy()
    op = rng.choice(5)
    rows = H.shape[0]

    if op == 0 and rows >= 1:
        # Scalar multiply: row *= nonzero element
        i = rng.integers(0, rows)
        scalar = GF4(rng.integers(1, 4))  # 1, ω, or ω²
        H[i] = H[i] * scalar

    elif op == 1 and rows >= 2:
        # Row addition: H[i] += scalar * H[j]
        i, j = rng.choice(rows, 2, replace=False)
        scalar = GF4(rng.integers(1, 4))
        H[i] = H[i] + scalar * H[j]

    elif op == 2:
        # Replace a row with a random Hermitian-orthogonal vector
        i = rng.integers(0, rows)
        for _ in range(200):
            new_row = GF4(rng.integers(0, 4, size=n))
            if np.all(new_row == 0):
                continue
            if hermitian_inner_product(new_row, new_row) != GF4(0):
                continue
            ok = True
            for j in range(rows):
                if j == i:
                    continue
                if hermitian_inner_product(new_row, H[j]) != GF4(0):
                    ok = False
                    break
            if ok:
                H[i] = new_row
                break

    elif op == 3:
        # Column permutation
        perm = rng.permutation(n)
        H = H[:, perm]

    elif op == 4:
        # Apply conjugation to one column: ω ↔ ω²
        j = rng.integers(0, n)
        H[:, j] = H[:, j]**2

    # Remove zero rows
    mask = np.array([not np.all(H[i] == 0) for i in range(H.shape[0])])
    H = H[mask]

    return H


def evolve_gf4(n, pop_size=50, num_gens=100, timeout_sec=15, rng=None):
    """Evolutionary search for GF(4) self-orthogonal codes."""
    if rng is None:
        rng = np.random.default_rng()

    from symplectic import StabilizerCode, compute_distance_qldpc

    print(f"Evolving GF(4) codes: n={n}, pop={pop_size}, gens={num_gens}")

    # Initialize
    pop = []
    fits = []

    for _ in range(pop_size * 5):
        target_rows = rng.integers(max(1, n//3), n)
        H = random_gf4_code(n, target_rows=target_rows, rng=rng)
        if H is None:
            continue

        symplectic_H = gf4_to_symplectic(H)
        code = StabilizerCode(H=symplectic_H, n=n)
        valid, _ = code.is_valid()
        if not valid or code.k < 1:
            continue

        d = compute_distance_qldpc(code, timeout_sec=timeout_sec)
        if d is None:
            d = 1

        pop.append(H)
        fits.append(d)
        if len(pop) >= pop_size:
            break

    if len(pop) < 5:
        print(f"  Failed to initialize (only {len(pop)} codes)")
        return None, 0

    best_d = max(fits)
    best_H = pop[np.argmax(fits)]
    best_k = None
    stagnation = 0

    print(f"  Init: {len(pop)} codes, best d={best_d}")

    for gen in range(num_gens):
        new_pop = []
        new_fits = []

        # Elitism
        elite = np.argsort(fits)[-3:]
        for idx in elite:
            new_pop.append(pop[idx])
            new_fits.append(fits[idx])

        while len(new_pop) < len(pop):
            # Tournament select
            idxs = rng.choice(len(pop), size=min(3, len(pop)), replace=False)
            parent = pop[idxs[np.argmax([fits[i] for i in idxs])]]

            # Mutate (2-4 ops)
            child_H = parent.copy()
            for _ in range(rng.integers(2, 5)):
                child_H = mutate_gf4(child_H, n, rng)

            if child_H.shape[0] < 1:
                continue

            if not is_hermitian_self_orthogonal(child_H):
                continue

            symplectic_H = gf4_to_symplectic(child_H)
            code = StabilizerCode(H=symplectic_H, n=n)
            valid, _ = code.is_valid()
            if not valid or code.k < 1:
                continue

            d = compute_distance_qldpc(code, timeout_sec=timeout_sec)
            if d is None:
                d = 1

            new_pop.append(child_H)
            new_fits.append(d)

        pop = new_pop
        fits = new_fits

        gen_best = max(fits)
        if gen_best > best_d:
            best_d = gen_best
            best_H = pop[np.argmax(fits)]

            symplectic_H = gf4_to_symplectic(best_H)
            code = StabilizerCode(H=symplectic_H, n=n)
            best_k = code.k
            stagnation = 0

            status = ''
            if (n, best_k) in KNOWN_BOUNDS:
                dl, du = KNOWN_BOUNDS[(n, best_k)]
                if best_d > dl:
                    status = f' *** IMPROVEMENT (known: {dl}) ***'
                elif best_d == dl:
                    status = f' MATCHES({dl})'
                else:
                    status = f' (vs {dl})'
            print(f"  Gen {gen}: [[{n},{best_k},{best_d}]]{status}")
        else:
            stagnation += 1

        if gen % 25 == 0 and gen > 0:
            print(f"  Gen {gen}: best d={best_d}, avg d={np.mean(fits):.1f}, stag={stagnation}")

    if best_k is None:
        symplectic_H = gf4_to_symplectic(best_H)
        code = StabilizerCode(H=symplectic_H, n=n)
        best_k = code.k

    print(f"  Final: [[{n},{best_k},{best_d}]]")
    return best_H, best_d


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    print("=" * 60)
    print("GF(4) ADDITIVE CODE SEARCH")
    print("=" * 60)

    # Validate [[5,1,3]]
    print("\nValidation: [[5,1,3]] perfect code")
    H_513 = GF4([
        [1, 2, 2, 1, 0],
        [0, 1, 2, 2, 1],
        [1, 0, 1, 2, 2],
        [2, 1, 0, 1, 2],
    ])
    print(f"  Self-orthogonal: {is_hermitian_self_orthogonal(H_513)}")

    from symplectic import StabilizerCode, compute_distance_qldpc
    sym = gf4_to_symplectic(H_513)
    code = StabilizerCode(H=sym, n=5)
    d = compute_distance_qldpc(code, timeout_sec=10)
    print(f"  [[5, {code.k}, {d}]] {'PASS' if d == 3 else 'FAIL'}")

    # Search n=10-16
    print("\n" + "=" * 60)
    print("EVOLUTIONARY GF(4) SEARCH")
    print("=" * 60)

    all_matches = []
    for n in range(10, 17):
        best_H, best_d = evolve_gf4(n, pop_size=40, num_gens=80, timeout_sec=12)
        if best_H is not None:
            sym = gf4_to_symplectic(best_H)
            code = StabilizerCode(H=sym, n=n)
            k = code.k
            if (n, k) in KNOWN_BOUNDS:
                dl, _ = KNOWN_BOUNDS[(n, k)]
                if best_d >= dl:
                    all_matches.append((n, k, best_d, 'MATCHES' if best_d == dl else 'IMPROVEMENT'))

    print(f"\n{'='*60}")
    print("RESULTS")
    for n, k, d, status in all_matches:
        print(f"  [[{n},{k},{d}]] {status}")
