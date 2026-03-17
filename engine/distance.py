"""
Code distance computation for stabilizer codes.

Three tiers:
1. Exact: enumerate all logical operators, find minimum weight (exponential, n < ~25)
2. Approximate: use stim's shortest_graphlike_error (fast, may overestimate)
3. Decoder-based: BP+OSD upper bound (medium speed, good estimate)
"""

import numpy as np
from itertools import product as cartesian_product
from codes import CSSCode, gf2_rank, gf2_nullspace
import time


def exact_distance_css(code: CSSCode, timeout_sec=300):
    """Compute exact code distance for a CSS code.

    For CSS codes, dx and dz can be computed independently:
    - dx = min weight of vector in ker(Hz) \\ rowspace(Hx)
    - dz = min weight of vector in ker(Hx) \\ rowspace(Hz)
    - d = min(dx, dz)

    This is exponential in the nullspace dimension, feasible for n < ~25.
    """
    start = time.time()
    n = code.n

    dx = _min_weight_coset(code.Hz, code.Hx, n, timeout_sec, start)
    if dx is None:
        return None  # timeout

    elapsed = time.time() - start
    remaining = timeout_sec - elapsed
    if remaining <= 0:
        return None

    dz = _min_weight_coset(code.Hx, code.Hz, n, remaining, start)
    if dz is None:
        return None

    return min(dx, dz)


def _min_weight_coset(H_kernel, H_stabilizer, n, timeout_sec, start_time):
    """Find minimum weight vector in ker(H_kernel) that is NOT in rowspace(H_stabilizer).

    This is the distance for one side (X or Z) of a CSS code.

    Optimized: enumerate all nullspace elements once, track minimum weight
    non-stabilizer vector. Uses batch computation for speed.
    """
    # Compute nullspace of H_kernel (vectors that commute with all checks)
    if H_kernel.size == 0:
        null_basis = np.eye(n, dtype=np.uint8)
    else:
        null_basis = gf2_nullspace(H_kernel)

    if null_basis.shape[0] == 0:
        return n  # No non-trivial logical operators

    null_dim = null_basis.shape[0]

    if null_dim > 22:
        # Too large for exact enumeration (~4M elements)
        return None

    # Pre-compute the stabilizer augmented matrix for fast rowspace check
    if H_stabilizer.size > 0:
        stab_rank = gf2_rank(H_stabilizer)
        # Row-reduce H_stabilizer for fast rowspace membership
        stab_rref = _row_reduce(H_stabilizer.copy())
    else:
        stab_rank = 0
        stab_rref = None

    min_weight = n + 1
    total = 2**null_dim

    # Process in chunks to allow timeout checks
    chunk_size = min(total, 50000)

    for chunk_start in range(1, total, chunk_size):
        if time.time() - start_time > timeout_sec:
            return None

        chunk_end = min(chunk_start + chunk_size, total)

        for bits in range(chunk_start, chunk_end):
            # Build the vector from nullspace coefficients
            coeffs = np.zeros(null_dim, dtype=np.uint8)
            b = bits
            for j in range(null_dim):
                coeffs[j] = b & 1
                b >>= 1

            vec = (coeffs @ null_basis) % 2
            w = int(np.sum(vec))

            if w == 0 or w >= min_weight:
                continue

            # Check if NOT in stabilizer rowspace
            if not _in_rowspace_fast(vec, stab_rref, stab_rank):
                min_weight = w
                if w == 1:
                    return 1  # Can't do better than 1

    return min_weight if min_weight <= n else n


def _row_reduce(M):
    """Row-reduce binary matrix to reduced row echelon form over GF(2)."""
    if M.size == 0:
        return M
    M = M.astype(np.uint8).copy()
    rows, cols = M.shape
    pivot_row = 0
    for col in range(cols):
        pivot = None
        for row in range(pivot_row, rows):
            if M[row, col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        M[[pivot_row, pivot]] = M[[pivot, pivot_row]]
        for row in range(rows):
            if row != pivot_row and M[row, col] == 1:
                M[row] = (M[row] + M[pivot_row]) % 2
        pivot_row += 1
    return M[:pivot_row]  # Return only non-zero rows


def _in_rowspace_fast(vec, rref, rank):
    """Fast check if vector is in rowspace using pre-computed RREF."""
    if rref is None or rank == 0:
        return False  # Empty rowspace
    aug = np.vstack([rref, vec.reshape(1, -1)])
    return gf2_rank(aug) == rank


def _in_rowspace(vec, H):
    """Check if binary vector is in the rowspace of H over GF(2)."""
    if H.size == 0:
        return False

    # Augment H with vec and check if rank changes
    aug = np.vstack([H, vec.reshape(1, -1)])
    return gf2_rank(aug) == gf2_rank(H)


def exact_distance_bruteforce(code: CSSCode, timeout_sec=300):
    """Alternative brute-force distance computation.

    Enumerate candidate logical operators by weight.
    Slower but simpler — useful as cross-check.
    """
    start = time.time()
    n = code.n

    # For CSS codes, check X-logicals and Z-logicals separately
    dx = _bruteforce_one_side(code.Hz, code.Hx, n, timeout_sec, start)
    if dx is None:
        return None

    elapsed = time.time() - start
    remaining = timeout_sec - elapsed
    if remaining <= 0:
        return None

    dz = _bruteforce_one_side(code.Hx, code.Hz, n, remaining, start)
    if dz is None:
        return None

    return min(dx, dz)


def _bruteforce_one_side(H_check, H_stab, n, timeout_sec, start_time):
    """Find min-weight vector in ker(H_check) \\ rowspace(H_stab) by weight enumeration."""
    from itertools import combinations

    for w in range(1, n + 1):
        if time.time() - start_time > timeout_sec:
            return None

        for positions in combinations(range(n), w):
            vec = np.zeros(n, dtype=np.uint8)
            vec[list(positions)] = 1

            # Check: in kernel of H_check?
            if H_check.size > 0:
                syndrome = (H_check @ vec) % 2
                if np.any(syndrome):
                    continue

            # Check: NOT in rowspace of H_stab?
            if not _in_rowspace(vec, H_stab):
                return w

    return n


def estimate_distance_random(code: CSSCode, num_samples=10000, rng=None):
    """Estimate distance by random sampling of logical operators.

    Returns an upper bound on the distance (may overestimate).
    Fast but imprecise — useful for filtering during search.
    """
    if rng is None:
        rng = np.random.default_rng()

    n = code.n

    # Compute nullspace bases
    null_hz = gf2_nullspace(code.Hz) if code.Hz.size > 0 else np.eye(n, dtype=np.uint8)
    null_hx = gf2_nullspace(code.Hx) if code.Hx.size > 0 else np.eye(n, dtype=np.uint8)

    min_weight = n + 1

    for null_basis, H_stab in [(null_hz, code.Hx), (null_hx, code.Hz)]:
        if null_basis.shape[0] == 0:
            continue

        for _ in range(num_samples):
            # Random combination of nullspace vectors
            coeffs = rng.integers(0, 2, size=null_basis.shape[0])
            if not np.any(coeffs):
                continue
            vec = (coeffs @ null_basis) % 2
            w = int(np.sum(vec))
            if w == 0:
                continue

            # Check: NOT in stabilizer rowspace?
            if not _in_rowspace(vec, H_stab):
                min_weight = min(min_weight, w)

    return min_weight if min_weight <= n else None
