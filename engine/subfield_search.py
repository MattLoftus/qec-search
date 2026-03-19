"""
Subfield subcode search for quantum error-correcting codes.

Idea: Build Hermitian self-orthogonal codes over GF(16) = GF(4^2),
then restrict to the GF(4) subfield. This can produce GF(4) codes
with parameters not achievable by direct GF(4) construction.

Theory:
  - GF(16) contains GF(4) as a subfield: {x in GF(16) : x^4 = x}
  - Hermitian conjugation over GF(16)/GF(4): x -> x^4 (Frobenius)
  - Hermitian inner product: <u,v> = sum(u_i * v_i^4)
  - If C over GF(16) satisfies C^perp_H <= C, the subfield subcode
    C|_{GF(4)} = {c in C : all entries in GF(4)} inherits self-orthogonality
  - The quantum code parameters [[n, k, d]] come from the GF(4) subfield
    subcode, potentially giving codes unreachable by direct GF(4) search

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

print = functools.partial(print, flush=True)
sys.path.insert(0, os.path.dirname(__file__))

from known_codes import KNOWN_BOUNDS
from gf4_codes import gf4_to_symplectic, GF4
from symplectic import StabilizerCode, compute_distance_qldpc

# Set nice priority
os.nice(15)

# ---------------------------------------------------------------------------
# GF(16) field setup
# ---------------------------------------------------------------------------

GF16 = galois.GF(16)
GF2 = galois.GF(2)

# GF(4) subfield of GF(16): elements x with x^4 = x
# In galois's GF(16), these are: {0, 1, 6, 7}
GF4_SUBFIELD = frozenset(int(x) for x in GF16.elements if x**4 == x)
assert len(GF4_SUBFIELD) == 4, f"Expected 4 subfield elements, got {len(GF4_SUBFIELD)}"

# Mapping from GF(16) subfield elements to GF(4) elements
# GF16(0) -> GF4(0), GF16(1) -> GF4(1), GF16(6) -> GF4(2)=omega, GF16(7) -> GF4(3)=omega^2
GF16_TO_GF4 = {0: 0, 1: 1, 6: 2, 7: 3}

# All nonzero GF(16) elements for random generation
GF16_NONZERO = [int(x) for x in GF16.elements if x != 0]

print(f"GF(4) subfield in GF(16): {sorted(GF4_SUBFIELD)}")
print(f"Mapping GF16->GF4: {GF16_TO_GF4}")


# ---------------------------------------------------------------------------
# Hermitian operations over GF(16)/GF(4)
# ---------------------------------------------------------------------------

def hermitian_conjugate_gf16(v):
    """Hermitian conjugation over GF(16)/GF(4): x -> x^4."""
    return v ** 4


def hermitian_inner_product_gf16(u, v):
    """Hermitian inner product over GF(16): <u,v> = sum(u_i * v_i^4)."""
    return np.sum(u * hermitian_conjugate_gf16(v))


def is_hermitian_self_orthogonal_gf16(H):
    """Check if rows of H over GF(16) are mutually Hermitian self-orthogonal."""
    H_conj = H ** 4
    product = H @ H_conj.T
    return np.all(product == GF16(0))


# ---------------------------------------------------------------------------
# Random Hermitian self-orthogonal code over GF(16)
# ---------------------------------------------------------------------------

def random_gf16_self_orthogonal(n, target_rows, rng):
    """Build a random Hermitian self-orthogonal code over GF(16).

    Adds rows one at a time, each constrained to be Hermitian-orthogonal
    to all existing rows and self-orthogonal.
    """
    rows = []
    for _ in range(target_rows):
        for _ in range(800):
            row = GF16(rng.integers(0, 16, size=n))
            if np.all(row == 0):
                continue
            # Self-orthogonal: <row, row>_H = 0
            if hermitian_inner_product_gf16(row, row) != GF16(0):
                continue
            # Orthogonal to all existing
            ok = True
            for existing in rows:
                if hermitian_inner_product_gf16(row, existing) != GF16(0):
                    ok = False
                    break
            if ok:
                rows.append(row)
                break

    if len(rows) < 1:
        return None
    return GF16(np.array([list(r) for r in rows]))


# ---------------------------------------------------------------------------
# Subfield subcode extraction
# ---------------------------------------------------------------------------

def extract_gf4_subfield_subcode(H_gf16):
    """Extract the GF(4) subfield subcode from a GF(16) code.

    The subfield subcode consists of all GF(4)-linear combinations
    of rows of H_gf16 whose entries all lie in GF(4) subset of GF(16).

    Strategy: express H_gf16 rows over a GF(4)-basis of GF(16),
    then find which GF(4)-linear combinations yield all-GF(4) entries.

    GF(16) is a 2-dimensional vector space over GF(4) with basis {1, beta}
    where beta is a root of x^2 + x + omega (irreducible over GF(4)).
    A codeword is in the subfield subcode iff its beta-component is zero.
    """
    rows, n = H_gf16.shape

    # Approach: enumerate GF(4)-linear combinations of rows
    # Each row over GF(16) can be scaled by GF(4) elements {0,1,omega,omega^2}
    # which live in GF(16) as {0, 1, 6, 7}
    gf4_in_gf16 = sorted(GF4_SUBFIELD)  # [0, 1, 6, 7]

    # For small number of rows, we can enumerate all GF(4)-linear combos
    # With r rows, there are 4^r combinations
    if rows > 12:
        # Too many to enumerate; use systematic approach
        return _extract_subfield_systematic(H_gf16)

    subfield_rows = []
    num_combos = 4 ** rows

    for combo_idx in range(1, num_combos):  # skip all-zero
        coeffs_idx = combo_idx
        coeffs = []
        for _ in range(rows):
            coeffs.append(gf4_in_gf16[coeffs_idx % 4])
            coeffs_idx //= 4

        # Compute linear combination
        vec = GF16(np.zeros(n, dtype=int))
        for j, c in enumerate(coeffs):
            if c != 0:
                vec = vec + GF16(c) * H_gf16[j]

        # Check if all entries are in GF(4) subfield
        if all(int(vec[i]) in GF4_SUBFIELD for i in range(n)):
            # Convert to GF(4)
            gf4_row = np.array([GF16_TO_GF4[int(vec[i])] for i in range(n)], dtype=int)
            subfield_rows.append(gf4_row)

    if not subfield_rows:
        return None

    # Remove duplicates and find a basis
    mat = GF4(np.array(subfield_rows))
    # Row reduce to find linearly independent rows
    mat = _gf4_row_reduce(mat)

    if mat.shape[0] < 1:
        return None
    return mat


def _extract_subfield_systematic(H_gf16):
    """Systematic subfield subcode extraction for large matrices.

    Decompose each GF(16) element as a + b*beta where a,b in GF(4).
    Find beta: an element of GF(16) not in GF(4).
    Express each row as (row_a, row_b) with row_a + row_b*beta.
    The subfield subcode is the kernel of the beta-component map.
    """
    rows, n = H_gf16.shape

    # Find a GF(4)-basis element for GF(16)/GF(4)
    # beta = primitive element of GF(16) works if it's not in GF(4)
    beta = GF16.primitive_element
    assert int(beta) not in GF4_SUBFIELD

    # Build the decomposition: for each element x in GF(16),
    # find a, b in GF(4) such that x = a + b*beta
    # We precompute this lookup table
    decomp_table = {}
    gf4_list = sorted(GF4_SUBFIELD)
    for a_val in gf4_list:
        for b_val in gf4_list:
            x = GF16(a_val) + GF16(b_val) * beta
            decomp_table[int(x)] = (a_val, b_val)

    assert len(decomp_table) == 16, "Decomposition should cover all 16 elements"

    # For each row, extract the beta-component
    # Subfield subcode = rows with all beta-components zero
    # For GF(4)-linear combinations: if row = (a_1,...,a_n) + beta*(b_1,...,b_n)
    # then lambda*row for lambda in GF(4) has beta-component = lambda*(b_1,...,b_n)
    # So we need combinations where the beta-components cancel

    # Build the beta-component matrix B where B[i,j] = b-component of H[i,j]
    B = np.zeros((rows, n), dtype=int)
    A = np.zeros((rows, n), dtype=int)
    for i in range(rows):
        for j in range(n):
            a_val, b_val = decomp_table[int(H_gf16[i, j])]
            A[i, j] = a_val
            B[i, j] = b_val

    # The subfield subcode is: find GF(4)-coefficients c_1,...,c_r such that
    # sum(c_i * B[i,:]) = 0 over GF(4)
    # This is a system of n equations over GF(4) in r unknowns

    # Convert B to GF(4) and find its kernel
    B_gf4 = GF4(B)
    # We want nullspace of B^T (treating rows as vectors)
    # Actually: for coefficients c (length r), we need c @ B = 0 over GF(4)
    # That means c is in left nullspace of B, i.e., nullspace of B.T
    B_T = B_gf4.T  # n x r matrix
    # Nullspace of B_T: find c such that B_T @ c = 0, i.e., B.T @ c = 0
    # Equivalently: c in nullspace of B (as row operation: c @ B = 0)

    # Use row reduction on B to find its left nullspace
    kernel = _gf4_left_nullspace(B_gf4)

    if kernel is None or kernel.shape[0] == 0:
        return None

    # Each kernel vector gives coefficients for a GF(4)-linear combination
    # that lies entirely in the GF(4) subfield
    subfield_rows = []
    for k_idx in range(kernel.shape[0]):
        coeffs = kernel[k_idx]  # GF(4) coefficients
        # Map to GF(16): GF4->GF16 embedding
        gf4_to_gf16 = {0: 0, 1: 1, 2: 6, 3: 7}
        vec = GF16(np.zeros(n, dtype=int))
        for i in range(rows):
            c_gf16 = GF16(gf4_to_gf16[int(coeffs[i])])
            if c_gf16 != GF16(0):
                vec = vec + c_gf16 * H_gf16[i]

        # Convert to GF(4)
        if all(int(vec[j]) in GF4_SUBFIELD for j in range(n)):
            gf4_row = np.array([GF16_TO_GF4[int(vec[j])] for j in range(n)], dtype=int)
            subfield_rows.append(gf4_row)

    if not subfield_rows:
        return None

    mat = GF4(np.array(subfield_rows))
    mat = _gf4_row_reduce(mat)
    if mat.shape[0] < 1:
        return None
    return mat


def _gf4_row_reduce(M):
    """Row reduce a GF(4) matrix and return only nonzero rows."""
    M = M.copy()
    rows, cols = M.shape
    pivot_row = 0
    for col in range(cols):
        # Find pivot
        found = -1
        for r in range(pivot_row, rows):
            if M[r, col] != GF4(0):
                found = r
                break
        if found == -1:
            continue
        # Swap
        if found != pivot_row:
            M[[pivot_row, found]] = M[[found, pivot_row]]
        # Scale pivot row
        inv = M[pivot_row, col] ** (-1)
        M[pivot_row] = M[pivot_row] * inv
        # Eliminate
        for r in range(rows):
            if r != pivot_row and M[r, col] != GF4(0):
                M[r] = M[r] + M[r, col] * M[pivot_row]  # GF(4) subtraction = addition
        pivot_row += 1

    # Remove zero rows
    nonzero = [i for i in range(rows) if np.any(M[i] != GF4(0))]
    if not nonzero:
        return GF4(np.zeros((0, cols), dtype=int))
    return M[nonzero]


def _gf4_left_nullspace(M):
    """Find the left nullspace of a GF(4) matrix M.

    Left nullspace: {c : c @ M = 0}, equivalently nullspace of M.T.
    """
    Mt = M.T.copy()  # cols_M x rows_M
    rows, cols = Mt.shape  # rows = original cols, cols = original rows

    # Augment with identity
    aug = GF4(np.hstack([Mt, GF4(np.eye(cols, dtype=int))]))

    pivot_row = 0
    for col in range(rows):
        found = -1
        for r in range(pivot_row, cols):
            if aug[r, col] != GF4(0):
                found = r
                break
        if found == -1:
            continue
        if found != pivot_row:
            aug[[pivot_row, found]] = aug[[found, pivot_row]]
        inv = aug[pivot_row, col] ** (-1)
        aug[pivot_row] = aug[pivot_row] * inv
        for r in range(cols):
            if r != pivot_row and aug[r, col] != GF4(0):
                aug[r] = aug[r] + aug[r, col] * aug[pivot_row]
        pivot_row += 1

    rank = pivot_row
    # Nullspace rows are those with zero in the first `rows` columns
    kernel_rows = []
    for r in range(cols):
        if np.all(aug[r, :rows] == GF4(0)):
            kernel_rows.append(aug[r, rows:])

    if not kernel_rows:
        return None
    return GF4(np.array([list(r) for r in kernel_rows]))


# ---------------------------------------------------------------------------
# Mutation operators for GF(16) codes
# ---------------------------------------------------------------------------

def mutate_gf16(H, n, rng):
    """Mutate a GF(16) Hermitian self-orthogonal matrix."""
    H = H.copy()
    op = rng.choice(6)
    rows = H.shape[0]

    if op == 0 and rows >= 1:
        # Scalar multiply row by nonzero GF(16) element
        i = rng.integers(0, rows)
        scalar = GF16(rng.choice(GF16_NONZERO))
        H[i] = H[i] * scalar

    elif op == 1 and rows >= 2:
        # Row addition: H[i] += scalar * H[j]
        i, j = rng.choice(rows, 2, replace=False)
        scalar = GF16(rng.choice(GF16_NONZERO))
        H[i] = H[i] + scalar * H[j]

    elif op == 2:
        # Replace a row with a random self-orthogonal vector
        i = rng.integers(0, rows)
        for _ in range(400):
            new_row = GF16(rng.integers(0, 16, size=n))
            if np.all(new_row == 0):
                continue
            if hermitian_inner_product_gf16(new_row, new_row) != GF16(0):
                continue
            ok = True
            for j in range(rows):
                if j == i:
                    continue
                if hermitian_inner_product_gf16(new_row, H[j]) != GF16(0):
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
        # Apply Frobenius x -> x^4 to one column
        j = rng.integers(0, n)
        H[:, j] = H[:, j] ** 4

    elif op == 5:
        # Apply x -> x^2 to one column (another automorphism)
        j = rng.integers(0, n)
        H[:, j] = H[:, j] ** 2

    # Remove zero rows
    mask = np.array([not np.all(H[i] == 0) for i in range(H.shape[0])])
    if np.any(mask):
        H = H[mask]
    else:
        return H  # Don't lose everything

    return H


# ---------------------------------------------------------------------------
# Evaluate a GF(16) code via its GF(4) subfield subcode
# ---------------------------------------------------------------------------

def evaluate_gf16_code(H_gf16, n, timeout_sec=15):
    """Extract GF(4) subfield subcode and compute quantum parameters.

    Returns (n, k, d) or None if the code is trivial.
    """
    # Extract subfield subcode
    H_gf4 = extract_gf4_subfield_subcode(H_gf16)
    if H_gf4 is None or H_gf4.shape[0] < 1:
        return None

    # Convert to symplectic binary
    symplectic_H = gf4_to_symplectic(H_gf4)
    code = StabilizerCode(H=symplectic_H, n=n)
    valid, _ = code.is_valid()
    if not valid:
        return None
    k = code.k
    if k < 1:
        return None

    d = compute_distance_qldpc(code, timeout_sec=timeout_sec)
    if d is None or d < 2:
        return None

    return n, k, d


# ---------------------------------------------------------------------------
# Main search: evolutionary + random
# ---------------------------------------------------------------------------

def search_subfield_subcodes(n, num_trials=200, evolve_gens=60,
                              pop_size=30, timeout_sec=12, rng=None):
    """Search for good quantum codes via GF(16) subfield subcodes at length n."""
    if rng is None:
        rng = np.random.default_rng()

    best_codes = {}  # (n, k) -> (d, H_gf16)
    codes_found = 0

    # Phase 1: Random sampling
    print(f"\n  [n={n}] Phase 1: Random sampling ({num_trials} trials)")
    for trial in range(num_trials):
        target_rows = rng.integers(max(1, n // 4), n)
        H_gf16 = random_gf16_self_orthogonal(n, target_rows, rng)
        if H_gf16 is None:
            continue
        if not is_hermitian_self_orthogonal_gf16(H_gf16):
            continue

        result = evaluate_gf16_code(H_gf16, n, timeout_sec=timeout_sec)
        if result is None:
            continue

        _, k, d = result
        codes_found += 1
        key = (n, k)

        if key not in best_codes or d > best_codes[key][0]:
            best_codes[key] = (d, H_gf16)

            status = ""
            if key in KNOWN_BOUNDS:
                dl, du = KNOWN_BOUNDS[key]
                if d > dl:
                    status = f" *** IMPROVEMENT (known d={dl}) ***"
                elif d == dl:
                    status = f" MATCHES(d={dl})"
                else:
                    status = f" (vs known d={dl})"
            print(f"    Trial {trial}: [[{n},{k},{d}]]{status}")

    if not best_codes:
        print(f"  [n={n}] No valid subfield subcodes found in random phase")
        return best_codes

    print(f"  [n={n}] Phase 1 found {codes_found} codes, {len(best_codes)} distinct (n,k)")

    # Phase 2: Evolutionary refinement on best seeds
    print(f"  [n={n}] Phase 2: Evolutionary refinement ({evolve_gens} gens, pop={pop_size})")

    for key in list(best_codes.keys()):
        seed_d, seed_H = best_codes[key]
        k = key[1]

        # Initialize population from seed + mutations
        population = [(seed_H, seed_d)]
        for _ in range(pop_size * 3):
            child = seed_H.copy()
            for _ in range(rng.integers(1, 4)):
                child = mutate_gf16(child, n, rng)
            if not is_hermitian_self_orthogonal_gf16(child):
                continue
            result = evaluate_gf16_code(child, n, timeout_sec=timeout_sec)
            if result is None:
                continue
            _, ck, cd = result
            population.append((child, cd))
            if len(population) >= pop_size:
                break

        if len(population) < 3:
            continue

        local_best_d = max(d for _, d in population)

        for gen in range(evolve_gens):
            # Sort by fitness
            population.sort(key=lambda x: x[1], reverse=True)
            # Keep top half
            survivors = population[:max(3, len(population) // 2)]

            new_pop = list(survivors)
            while len(new_pop) < pop_size:
                # Tournament select
                idx = rng.integers(0, len(survivors))
                parent_H = survivors[idx][0]

                child = parent_H.copy()
                for _ in range(rng.integers(1, 5)):
                    child = mutate_gf16(child, n, rng)

                if not is_hermitian_self_orthogonal_gf16(child):
                    continue
                result = evaluate_gf16_code(child, n, timeout_sec=timeout_sec)
                if result is None:
                    continue
                _, ck, cd = result
                new_pop.append((child, cd))

            population = new_pop
            gen_best = max(d for _, d in population)
            if gen_best > local_best_d:
                local_best_d = gen_best
                best_H = max(population, key=lambda x: x[1])[0]

                # Re-evaluate to get correct k
                result = evaluate_gf16_code(best_H, n, timeout_sec=timeout_sec)
                if result:
                    _, ek, ed = result
                    ekey = (n, ek)
                    if ekey not in best_codes or ed > best_codes[ekey][0]:
                        best_codes[ekey] = (ed, best_H)
                        status = ""
                        if ekey in KNOWN_BOUNDS:
                            dl, du = KNOWN_BOUNDS[ekey]
                            if ed > dl:
                                status = f" *** IMPROVEMENT ***"
                            elif ed == dl:
                                status = f" MATCHES(d={dl})"
                            else:
                                status = f" (vs d={dl})"
                        print(f"    Evolve gen {gen}: [[{n},{ek},{ed}]]{status}")

    return best_codes


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate():
    """Validate with known constructions."""
    print("=" * 60)
    print("VALIDATION: Subfield subcode machinery")
    print("=" * 60)

    rng = np.random.default_rng(42)

    # Test 1: Hermitian self-orthogonality over GF(16)
    print("\nTest 1: Random self-orthogonal codes over GF(16)")
    for n in [5, 7, 10]:
        H = random_gf16_self_orthogonal(n, target_rows=n - 1, rng=rng)
        if H is not None:
            ok = is_hermitian_self_orthogonal_gf16(H)
            print(f"  n={n}: {H.shape[0]} rows, self-orthogonal={ok}")
        else:
            print(f"  n={n}: failed to build")

    # Test 2: Subfield subcode extraction
    print("\nTest 2: Subfield subcode extraction")
    for n in [5, 7]:
        H = random_gf16_self_orthogonal(n, target_rows=n - 1, rng=rng)
        if H is None:
            print(f"  n={n}: no GF(16) code")
            continue
        sub = extract_gf4_subfield_subcode(H)
        if sub is not None:
            print(f"  n={n}: GF(16) code {H.shape[0]}x{n} -> GF(4) subfield subcode {sub.shape[0]}x{n}")
            result = evaluate_gf16_code(H, n, timeout_sec=10)
            if result:
                print(f"    Quantum code: [[{result[0]},{result[1]},{result[2]}]]")
            else:
                print(f"    No valid quantum code")
        else:
            print(f"  n={n}: empty subfield subcode")

    # Test 3: Quick search at n=6
    print("\nTest 3: Quick search at n=6")
    results = search_subfield_subcodes(6, num_trials=50, evolve_gens=10,
                                        pop_size=10, timeout_sec=8, rng=rng)
    for (sn, sk), (sd, _) in sorted(results.items()):
        status = ""
        if (sn, sk) in KNOWN_BOUNDS:
            dl, _ = KNOWN_BOUNDS[(sn, sk)]
            if sd >= dl:
                status = " MATCHES" if sd == dl else " IMPROVEMENT!"
        print(f"  [[{sn},{sk},{sd}]]{status}")

    print("\nValidation complete.\n")


# ---------------------------------------------------------------------------
# Main search
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Global timeout
    GLOBAL_TIMEOUT = 1800  # 30 minutes

    def alarm_handler(signum, frame):
        print("\n" + "=" * 60)
        print("GLOBAL TIMEOUT REACHED")
        print("=" * 60)
        # Print summary before exiting
        print_summary()
        sys.exit(0)

    all_results = {}

    def print_summary():
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

    signal.signal(signal.SIGALRM, alarm_handler)
    signal.alarm(GLOBAL_TIMEOUT)

    t_start = time.time()

    # Run validation first
    validate()

    # Main search
    print("=" * 60)
    print("SUBFIELD SUBCODE SEARCH: n=5..20")
    print("=" * 60)

    rng = np.random.default_rng()

    for n in range(5, 21):
        t_n = time.time()
        print(f"\n{'='*50}")
        print(f"SEARCHING n={n}")
        print(f"{'='*50}")

        # Scale effort with n
        if n <= 8:
            num_trials = 300
            evolve_gens = 40
            pop_size = 25
        elif n <= 12:
            num_trials = 200
            evolve_gens = 50
            pop_size = 25
        elif n <= 16:
            num_trials = 150
            evolve_gens = 40
            pop_size = 20
        else:
            num_trials = 100
            evolve_gens = 30
            pop_size = 15

        results = search_subfield_subcodes(
            n, num_trials=num_trials, evolve_gens=evolve_gens,
            pop_size=pop_size, timeout_sec=12, rng=rng
        )

        for key, val in results.items():
            if key not in all_results or val[0] > all_results[key][0]:
                all_results[key] = val

        elapsed = time.time() - t_n
        print(f"  n={n} complete in {elapsed:.1f}s")

    signal.alarm(0)
    total = time.time() - t_start
    print(f"\nTotal time: {total:.1f}s")

    print_summary()
