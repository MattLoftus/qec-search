"""
General (non-CSS) stabilizer code representation and search.

Non-CSS stabilizer codes use the full symplectic representation:
- An [[n, k, d]] code is defined by (n-k) independent commuting Pauli operators
- Each operator is a row in a (n-k) x 2n binary matrix H = [X | Z]
  where the left half encodes X components and right half encodes Z components
- Row i represents the Pauli: X^{H[i,0:n]} Z^{H[i,n:2n]}
- Commutativity: rows i,j commute iff H[i] * Omega * H[j]^T = 0 (mod 2)
  where Omega = [[0, I], [I, 0]] is the symplectic form

This is a strictly larger space than CSS codes (which have H = [Hx 0; 0 Hz]).
Many of the best known codes in codetables.de are non-CSS.

Usage:
    cd ~/workspace/qec-search && venv/bin/python3 engine/symplectic.py
"""

import numpy as np
import time
import signal
import functools
from dataclasses import dataclass, field
from typing import Optional, List, Tuple

print = functools.partial(print, flush=True)

import sys
import os
sys.path.insert(0, os.path.dirname(__file__))
from codes import gf2_rank, gf2_nullspace


# ---------------------------------------------------------------------------
# Symplectic inner product and commutativity
# ---------------------------------------------------------------------------

def symplectic_product(a, b, n):
    """Compute the symplectic inner product of two 2n-bit vectors.

    For Pauli operators represented as [x1..xn | z1..zn]:
    <a, b>_s = sum(a_x * b_z + a_z * b_x) mod 2

    Two Paulis commute iff their symplectic inner product is 0.
    """
    ax, az = a[:n], a[n:]
    bx, bz = b[:n], b[n:]
    return int(np.sum(ax * bz + az * bx)) % 2


def check_commutation(H, n):
    """Check that all rows of the stabilizer matrix H commute.

    H is (n-k) x 2n binary matrix.
    Returns True if all pairs of rows have symplectic inner product 0.
    """
    rows = H.shape[0]
    for i in range(rows):
        for j in range(i + 1, rows):
            if symplectic_product(H[i], H[j], n) != 0:
                return False
    return True


def symplectic_matrix_commutation(H, n):
    """Compute the full commutation matrix.

    Returns a binary matrix C where C[i,j] = <H[i], H[j]>_s.
    The code is valid iff C = 0.
    """
    Hx = H[:, :n]
    Hz = H[:, n:]
    # C = Hx @ Hz.T + Hz @ Hx.T (mod 2)
    C = (Hx @ Hz.T + Hz @ Hx.T) % 2
    return C


# ---------------------------------------------------------------------------
# StabilizerCode class
# ---------------------------------------------------------------------------

@dataclass
class StabilizerCode:
    """A general [[n, k, d]] stabilizer code.

    Defined by a (n-k) x 2n binary matrix H in symplectic form [X | Z].
    Each row is a stabilizer generator.
    """
    H: np.ndarray   # (n-k, 2n) binary matrix
    n: int           # number of physical qubits
    _distance: Optional[int] = field(default=None, repr=False)

    def __post_init__(self):
        self.H = np.asarray(self.H, dtype=np.uint8)
        assert self.H.shape[1] == 2 * self.n

    @property
    def num_generators(self):
        return self.H.shape[0]

    @property
    def k(self):
        """Number of logical qubits."""
        return self.n - gf2_rank(self.H)

    @property
    def distance(self):
        return self._distance

    @distance.setter
    def distance(self, d):
        self._distance = d

    @property
    def is_css(self):
        """Check if this code is CSS (no mixed X/Z generators)."""
        Hx = self.H[:, :self.n]
        Hz = self.H[:, self.n:]
        for i in range(self.H.shape[0]):
            has_x = np.any(Hx[i])
            has_z = np.any(Hz[i])
            if has_x and has_z:
                return False
        return True

    def is_valid(self):
        """Check that all stabilizers commute."""
        C = symplectic_matrix_commutation(self.H, self.n)
        if np.any(C):
            return False, "Non-commuting stabilizers"
        return True, "Valid stabilizer code"

    @property
    def weight_profile(self):
        """Weight of each stabilizer generator."""
        Hx = self.H[:, :self.n]
        Hz = self.H[:, self.n:]
        # Weight of a Pauli = number of qubits it acts on nontrivially
        weights = []
        for i in range(self.H.shape[0]):
            w = int(np.sum(np.logical_or(Hx[i], Hz[i])))
            weights.append(w)
        return weights

    def __repr__(self):
        d_str = str(self._distance) if self._distance is not None else "?"
        css_str = " CSS" if self.is_css else ""
        return f"StabilizerCode([[{self.n}, {self.k}, {d_str}]]{css_str})"


# ---------------------------------------------------------------------------
# Distance computation for non-CSS codes
# ---------------------------------------------------------------------------

def compute_distance_qldpc(code: StabilizerCode, timeout_sec=60):
    """Compute distance using qLDPC's optimized QuditCode (preferred).

    Uses bit-packed vectorized distance computation. Fast for n <= ~30.
    Falls back to our own implementation for larger codes or if qLDPC unavailable.
    """
    class TO(Exception):
        pass

    def handler(s, f):
        raise TO()

    try:
        import galois
        from qldpc.codes.quantum import QuditCode
        GF2 = galois.GF(2)

        old = signal.signal(signal.SIGALRM, handler)
        signal.alarm(timeout_sec)
        try:
            qcode = QuditCode(GF2(code.H.astype(int)))
            _, _, d = qcode.get_code_params()
            signal.alarm(0)
            return int(d)
        except TO:
            signal.alarm(0)
            return None
        except Exception:
            signal.alarm(0)
            return compute_distance(code, timeout_sec)
        finally:
            signal.signal(signal.SIGALRM, old)
    except ImportError:
        return compute_distance(code, timeout_sec)


def compute_distance(code: StabilizerCode, timeout_sec=60):
    """Compute the distance of a general stabilizer code (fallback).

    The distance is the minimum weight of a Pauli operator that:
    1. Commutes with all stabilizers (is in the centralizer)
    2. Is NOT itself a stabilizer (is a nontrivial logical operator)

    For non-CSS codes, we can't separate X and Z — we must search
    the full 2n-dimensional symplectic space.
    """
    start = time.time()
    n = code.n
    H = code.H

    # Find the centralizer: all 2n-bit vectors v such that
    # symplectic_product(v, H[i]) = 0 for all i
    # This is the nullspace of H @ Omega^T (or equivalently, the
    # symplectic complement of the stabilizer group)

    # Build the symplectic form matrix
    # For row v = [vx | vz], the condition <v, H[i]>_s = 0 means:
    # vx @ Hz[i]^T + vz @ Hx[i]^T = 0 (mod 2)
    # This is equivalent to: v @ [Hz^T; Hx^T] = 0, i.e., v @ Omega @ H^T = 0

    Hx = H[:, :n]
    Hz = H[:, n:]
    # Omega @ H^T has rows [Hz^T | Hx^T]^T → we need the matrix M such that
    # v @ M = 0 for v in centralizer
    # M = [Hz | Hx]^T  (swap X and Z parts, then transpose)
    M = np.hstack([Hz, Hx]).T  # 2n x (n-k)
    # But we want vectors v such that v @ M = 0 mod 2
    # Equivalently, M^T @ v^T = 0, so v is in nullspace of M^T = [Hz | Hx]
    M_check = np.hstack([Hz, Hx])  # (n-k) x 2n

    centralizer_basis = gf2_nullspace(M_check)

    if centralizer_basis.shape[0] == 0:
        return n  # No centralizer elements beyond stabilizers

    # The stabilizer group itself is spanned by H
    # A nontrivial logical is in centralizer but NOT in rowspace of H
    # Distance = min weight of such an operator

    cent_dim = centralizer_basis.shape[0]
    if cent_dim > 22:
        return None  # Too large for exact enumeration

    min_weight = 2 * n + 1

    total = 2**cent_dim
    for bits in range(1, total):
        if time.time() - start > timeout_sec:
            return None

        # Build vector from centralizer basis
        coeffs = np.zeros(cent_dim, dtype=np.uint8)
        b = bits
        for j in range(cent_dim):
            coeffs[j] = b & 1
            b >>= 1

        vec = (coeffs @ centralizer_basis) % 2

        # Compute Pauli weight: number of qubits where X or Z is nonzero
        vx = vec[:n]
        vz = vec[n:]
        w = int(np.sum(np.logical_or(vx, vz)))

        if w == 0 or w >= min_weight:
            continue

        # Check: NOT in rowspace of H (not a stabilizer)
        aug = np.vstack([H, vec.reshape(1, -1)])
        if gf2_rank(aug) == gf2_rank(H):
            continue  # It's a stabilizer — skip

        min_weight = w
        if w == 1:
            return 1

    return min_weight if min_weight <= 2 * n else n


# ---------------------------------------------------------------------------
# Known non-CSS codes
# ---------------------------------------------------------------------------

def perfect_5_1_3():
    """The [[5, 1, 3]] perfect code — smallest QEC code.

    Stabilizer generators:
    XZZXI
    IXZZX
    XIXZZ
    ZXIXZ
    """
    n = 5
    H = np.zeros((4, 10), dtype=np.uint8)
    # Row 0: XZZXI → X=[1,0,0,0,0] Z=[0,1,1,0,0] → but wait, Y=XZ
    # Actually: X component marks X or Y, Z component marks Z or Y
    # XZZXI: qubit 0=X, qubit 1=Z, qubit 2=Z, qubit 3=X, qubit 4=I
    H[0] = [1, 0, 0, 1, 0,  0, 1, 1, 0, 0]  # XZZXI
    H[1] = [0, 1, 0, 0, 1,  0, 0, 1, 1, 0]  # IXZZX
    H[2] = [1, 0, 1, 0, 0,  0, 0, 0, 1, 1]  # XIXZZ
    H[3] = [0, 1, 0, 1, 0,  1, 0, 0, 0, 1]  # ZXIXZ
    return StabilizerCode(H=H, n=n)


# ---------------------------------------------------------------------------
# Random non-CSS code generation
# ---------------------------------------------------------------------------

def random_stabilizer_code(n, target_k=1, rng=None):
    """Generate a random stabilizer code on n qubits.

    Strategy: build rows one at a time, each constrained to
    symplectically commute with all previous rows.
    """
    if rng is None:
        rng = np.random.default_rng()

    num_generators = n - target_k
    rows = []

    for _ in range(num_generators):
        # Try to find a row that commutes with all existing rows
        for _ in range(200):
            row = rng.integers(0, 2, size=2 * n).astype(np.uint8)
            if not np.any(row):
                continue

            # Check commutativity with all existing rows
            ok = True
            for existing in rows:
                if symplectic_product(row, existing, n) != 0:
                    ok = False
                    break
            if ok:
                rows.append(row)
                break

    if len(rows) < max(1, num_generators - 1):
        return None

    H = np.array(rows, dtype=np.uint8)
    code = StabilizerCode(H=H, n=n)
    valid, _ = code.is_valid()
    if not valid:
        return None
    if code.k < 1:
        return None

    return code


def random_non_css_code(n, target_k=1, rng=None):
    """Generate a random non-CSS stabilizer code.

    Like random_stabilizer_code but ensures at least one generator
    has both X and Z components (genuinely non-CSS).
    """
    if rng is None:
        rng = np.random.default_rng()

    for _ in range(100):
        code = random_stabilizer_code(n, target_k=target_k, rng=rng)
        if code is None:
            continue
        if not code.is_css:
            return code

    return None


# ---------------------------------------------------------------------------
# Validation and testing
# ---------------------------------------------------------------------------

def validate_known_codes():
    """Validate the [[5,1,3]] perfect code."""
    print("=" * 60)
    print("NON-CSS CODE VALIDATION")
    print("=" * 60)

    code = perfect_5_1_3()
    print(f"\n[[5,1,3]] perfect code:")
    print(f"  Repr: {code}")
    valid, msg = code.is_valid()
    print(f"  Valid: {valid} ({msg})")
    print(f"  n={code.n}, k={code.k}")
    print(f"  CSS: {code.is_css}")
    print(f"  Weights: {code.weight_profile}")

    t0 = time.time()
    d = compute_distance(code, timeout_sec=30)
    elapsed = time.time() - t0
    print(f"  Distance: {d} (expected 3) [{'PASS' if d == 3 else 'FAIL'}] [{elapsed:.2f}s]")

    # Test random generation
    print(f"\nRandom non-CSS code generation (n=8, k=1):")
    rng = np.random.default_rng(42)
    d_counts = {}
    generated = 0
    for _ in range(500):
        code = random_non_css_code(8, target_k=1, rng=rng)
        if code is None:
            continue
        generated += 1
        d = compute_distance(code, timeout_sec=5)
        if d is not None:
            d_counts[d] = d_counts.get(d, 0) + 1

    print(f"  Generated {generated} valid non-CSS codes")
    print(f"  Distance distribution: {dict(sorted(d_counts.items()))}")

    # Compare with CSS codes at same parameters
    print(f"\nRandom CSS codes (n=8, k=1) for comparison:")
    from codes import random_css_code
    from distance import exact_distance_css
    css_d_counts = {}
    for _ in range(500):
        css = random_css_code(8, target_k=1, rng=rng)
        if css is None or css.k < 1:
            continue
        d = exact_distance_css(css, timeout_sec=5)
        if d is not None:
            css_d_counts[d] = css_d_counts.get(d, 0) + 1

    print(f"  Distance distribution: {dict(sorted(css_d_counts.items()))}")


# ---------------------------------------------------------------------------
# Clifford mutation operators (commutativity-preserving)
# ---------------------------------------------------------------------------

def mutate_symplectic(H, n, rng):
    """Mutate a stabilizer matrix using Clifford-group operations.

    All operations preserve commutativity by construction.
    """
    H = H.copy()
    op = rng.choice(5)

    if op == 0 and H.shape[0] >= 2:
        # Row addition: S[i] += S[j] (bilinearity preserves commutativity)
        i, j = rng.choice(H.shape[0], 2, replace=False)
        H[i] = (H[i] + H[j]) % 2

    elif op == 1:
        # Row replacement from nullspace of other rows' symplectic complement
        i = rng.integers(0, H.shape[0])
        others = np.delete(H, i, axis=0)
        if others.shape[0] > 0:
            constraint = np.hstack([others[:, n:], others[:, :n]])
            null = gf2_nullspace(constraint)
            if null.shape[0] > 0:
                coeffs = rng.integers(0, 2, size=null.shape[0])
                if not np.any(coeffs):
                    coeffs[rng.integers(0, null.shape[0])] = 1
                new_row = (coeffs @ null) % 2
                if np.any(new_row):
                    H[i] = new_row

    elif op == 2:
        # Qubit permutation (relabeling)
        perm = rng.permutation(n)
        H[:, :n] = H[:, perm]
        H[:, n:] = H[:, n + perm]

    elif op == 3:
        # Local Hadamard on qubit j: swap X_j <-> Z_j
        j = rng.integers(0, n)
        H[:, [j, n + j]] = H[:, [n + j, j]]

    elif op == 4:
        # Local S gate on qubit j: X_j -> Y_j (toggle Z where X=1)
        j = rng.integers(0, n)
        mask = H[:, j] == 1
        H[mask, n + j] = (H[mask, n + j] + 1) % 2

    return H


def crossover_symplectic(H1, H2, n, rng):
    """Crossover two stabilizer codes.

    Take rows from both parents, keeping only rows that commute.
    """
    all_rows = np.vstack([H1, H2])
    rng.shuffle(all_rows)

    selected = []
    target = max(H1.shape[0], H2.shape[0])

    for row in all_rows:
        if not np.any(row):
            continue
        ok = True
        for existing in selected:
            if symplectic_product(row, existing, n) != 0:
                ok = False
                break
        if ok:
            selected.append(row)
            if len(selected) >= target:
                break

    if len(selected) < 2:
        return None
    return np.array(selected, dtype=np.uint8)


# ---------------------------------------------------------------------------
# Evolutionary non-CSS search
# ---------------------------------------------------------------------------

def evolve_non_css(n, target_k=1, pop_size=50, num_gens=100,
                   timeout_sec=15, rng=None):
    """Evolutionary search for non-CSS stabilizer codes.

    Uses Clifford mutations (Hadamard, S gate, row ops, qubit permutation)
    which preserve commutativity by construction.
    """
    if rng is None:
        rng = np.random.default_rng()

    from known_codes import KNOWN_BOUNDS

    # Initialize population
    population = []
    fitnesses = []

    print(f"Evolving non-CSS codes: n={n}, k={target_k}, pop={pop_size}, gens={num_gens}")

    for _ in range(pop_size * 5):
        code = random_stabilizer_code(n, target_k=target_k, rng=rng)
        if code is None or code.k < 1:
            continue
        d = compute_distance_qldpc(code, timeout_sec=timeout_sec)
        if d is None:
            d = 1
        code.distance = d
        population.append(code)
        fitnesses.append(d)
        if len(population) >= pop_size:
            break

    if len(population) < 2:
        print("Failed to initialize population")
        return None

    best_d = max(fitnesses)
    best_code = population[np.argmax(fitnesses)]
    stagnation = 0

    print(f"Init: {len(population)} codes, best d={best_d}")

    for gen in range(num_gens):
        # Selection + reproduction
        new_pop = []
        new_fit = []

        # Elitism
        elite_idx = np.argsort(fitnesses)[-3:]
        for idx in elite_idx:
            new_pop.append(population[idx])
            new_fit.append(fitnesses[idx])

        while len(new_pop) < pop_size:
            # Tournament select parent
            idxs = rng.choice(len(population), size=3, replace=False)
            parent_idx = idxs[np.argmax([fitnesses[i] for i in idxs])]
            parent = population[parent_idx]

            # Mutate
            new_H = mutate_symplectic(parent.H, n, rng)

            # Sometimes double-mutate
            if rng.random() < 0.5:
                new_H = mutate_symplectic(new_H, n, rng)

            # Remove zero rows
            new_H = new_H[np.any(new_H, axis=1)]
            if new_H.shape[0] < 1:
                continue

            child = StabilizerCode(H=new_H, n=n)
            valid, _ = child.is_valid()
            if not valid or child.k < 1:
                continue

            d = compute_distance_qldpc(child, timeout_sec=timeout_sec)
            if d is None:
                d = 1
            child.distance = d
            new_pop.append(child)
            new_fit.append(d)

        population = new_pop
        fitnesses = new_fit

        gen_best_d = max(fitnesses)
        if gen_best_d > best_d:
            best_d = gen_best_d
            best_code = population[np.argmax(fitnesses)]
            stagnation = 0

            key = (n, target_k)
            status = ''
            if key in KNOWN_BOUNDS:
                dl, du = KNOWN_BOUNDS[key]
                if best_d > dl:
                    status = ' *** IMPROVEMENT ***'
                elif best_d == dl:
                    status = f' MATCHES({dl})'
                else:
                    status = f' (vs {dl})'
            print(f"  Gen {gen}: NEW BEST [[{n},{best_code.k},{best_d}]]{status}")
        else:
            stagnation += 1

        if gen % 20 == 0:
            print(f"  Gen {gen}: best d={best_d}, avg d={np.mean(fitnesses):.1f}, stag={stagnation}")

    print(f"Final: [[{n}, {best_code.k}, {best_d}]]")
    return best_code


if __name__ == "__main__":
    validate_known_codes()

    # Quick evolutionary test
    print(f"\n{'='*60}")
    print("EVOLUTIONARY NON-CSS SEARCH (n=10, k=1)")
    print(f"{'='*60}")
    best = evolve_non_css(10, target_k=1, pop_size=30, num_gens=50, timeout_sec=10)
