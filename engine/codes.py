"""
Stabilizer code representation and construction.

Core abstractions:
- CSSCode: CSS stabilizer code from Hx, Hz parity check matrices
- StabilizerCode: General stabilizer code from symplectic matrix (future)
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional
import hashlib
import json


def gf2_rank(matrix):
    """Compute rank of a binary matrix over GF(2) via row echelon form."""
    if matrix.size == 0:
        return 0
    M = matrix.astype(np.uint8).copy()
    rows, cols = M.shape
    rank = 0
    for col in range(cols):
        # Find pivot row
        pivot = None
        for row in range(rank, rows):
            if M[row, col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        # Swap pivot row to current rank position
        M[[rank, pivot]] = M[[pivot, rank]]
        # Eliminate below
        for row in range(rows):
            if row != rank and M[row, col] == 1:
                M[row] = (M[row] + M[rank]) % 2
        rank += 1
    return rank


def gf2_nullspace(matrix):
    """Compute nullspace of a binary matrix over GF(2).

    Returns a matrix whose rows span the nullspace.
    """
    if matrix.size == 0:
        return np.zeros((matrix.shape[1], matrix.shape[1]), dtype=np.uint8)

    M = matrix.astype(np.uint8).copy()
    rows, cols = M.shape

    # Augment with identity
    aug = np.hstack([M, np.eye(rows, dtype=np.uint8)])

    # Row echelon
    rank = 0
    pivot_cols = []
    for col in range(cols):
        pivot = None
        for row in range(rank, rows):
            if aug[row, col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        aug[[rank, pivot]] = aug[[pivot, rank]]
        for row in range(rows):
            if row != rank and aug[row, col] == 1:
                aug[row] = (aug[row] + aug[rank]) % 2
        pivot_cols.append(col)
        rank += 1

    # Free columns (not pivot columns) give nullspace vectors
    free_cols = [c for c in range(cols) if c not in pivot_cols]
    if not free_cols:
        return np.zeros((0, cols), dtype=np.uint8)

    # Build nullspace basis
    null_basis = []
    for fc in free_cols:
        vec = np.zeros(cols, dtype=np.uint8)
        vec[fc] = 1
        for i, pc in enumerate(pivot_cols):
            vec[pc] = aug[i, fc]
        null_basis.append(vec)

    return np.array(null_basis, dtype=np.uint8)


@dataclass
class CSSCode:
    """A CSS (Calderbank-Shor-Steane) stabilizer code.

    Defined by two binary parity check matrices Hx and Hz where:
    - n = number of physical qubits (columns)
    - Hx rows are Z-type stabilizer generators (detect X errors)
    - Hz rows are X-type stabilizer generators (detect Z errors)
    - Commutativity: Hx @ Hz.T = 0 (mod 2)
    - k = n - rank(Hx) - rank(Hz)
    """
    Hx: np.ndarray  # (rx, n) binary matrix
    Hz: np.ndarray  # (rz, n) binary matrix
    _distance: Optional[int] = field(default=None, repr=False)
    _hash: Optional[str] = field(default=None, repr=False)

    def __post_init__(self):
        self.Hx = np.asarray(self.Hx, dtype=np.uint8)
        self.Hz = np.asarray(self.Hz, dtype=np.uint8)

    @property
    def n(self):
        """Number of physical qubits."""
        return self.Hx.shape[1]

    @property
    def k(self):
        """Number of logical qubits."""
        return self.n - gf2_rank(self.Hx) - gf2_rank(self.Hz)

    @property
    def rx(self):
        """Number of X-check rows."""
        return self.Hx.shape[0]

    @property
    def rz(self):
        """Number of Z-check rows."""
        return self.Hz.shape[0]

    @property
    def distance(self):
        """Code distance (cached). Must be computed externally and set."""
        return self._distance

    @distance.setter
    def distance(self, d):
        self._distance = d

    @property
    def params(self):
        """[[n, k, d]] parameters as tuple."""
        return (self.n, self.k, self._distance)

    @property
    def rate(self):
        """Code rate k/n."""
        return self.k / self.n if self.n > 0 else 0

    @property
    def max_stabilizer_weight(self):
        """Maximum weight (number of non-zero entries) across all stabilizers."""
        weights = []
        if self.Hx.size > 0:
            weights.extend(np.sum(self.Hx, axis=1).tolist())
        if self.Hz.size > 0:
            weights.extend(np.sum(self.Hz, axis=1).tolist())
        return max(weights) if weights else 0

    @property
    def avg_stabilizer_weight(self):
        """Average stabilizer weight."""
        weights = []
        if self.Hx.size > 0:
            weights.extend(np.sum(self.Hx, axis=1).tolist())
        if self.Hz.size > 0:
            weights.extend(np.sum(self.Hz, axis=1).tolist())
        return float(np.mean(weights)) if weights else 0

    def code_hash(self):
        """Deterministic hash for deduplication."""
        if self._hash is None:
            data = {
                "Hx": self.Hx.tolist(),
                "Hz": self.Hz.tolist(),
            }
            self._hash = hashlib.sha256(
                json.dumps(data, sort_keys=True).encode()
            ).hexdigest()[:16]
        return self._hash

    def is_valid(self):
        """Check CSS commutativity constraint: Hx @ Hz.T = 0 (mod 2)."""
        if self.Hx.shape[1] != self.Hz.shape[1]:
            return False, "Hx and Hz have different number of columns"
        if self.Hx.size == 0 or self.Hz.size == 0:
            return True, "Trivial (empty check matrix)"
        product = (self.Hx @ self.Hz.T) % 2
        if np.any(product):
            return False, "Commutativity violated: Hx @ Hz.T != 0 (mod 2)"
        return True, "Valid CSS code"

    def to_dict(self):
        """Serialize to JSON-safe dict."""
        return {
            "type": "css",
            "n": int(self.n),
            "k": int(self.k),
            "d": int(self._distance) if self._distance is not None else None,
            "Hx": self.Hx.tolist(),
            "Hz": self.Hz.tolist(),
            "hash": self.code_hash(),
            "max_weight": int(self.max_stabilizer_weight),
            "avg_weight": float(self.avg_stabilizer_weight),
            "rate": float(self.rate),
        }

    @classmethod
    def from_dict(cls, d):
        """Deserialize from dict."""
        code = cls(
            Hx=np.array(d["Hx"], dtype=np.uint8),
            Hz=np.array(d["Hz"], dtype=np.uint8),
        )
        if d.get("d") is not None:
            code._distance = d["d"]
        return code

    def __repr__(self):
        d_str = str(self._distance) if self._distance is not None else "?"
        return f"CSSCode([[{self.n}, {self.k}, {d_str}]])"


# ---------------------------------------------------------------------------
# CSS code generation utilities
# ---------------------------------------------------------------------------

def random_css_code(n, target_k=None, density=0.3, rng=None):
    """Generate a random CSS code on n qubits.

    Strategy: generate random Hx, then construct Hz to satisfy commutativity.

    Args:
        n: number of physical qubits
        target_k: desired number of logical qubits (approximate)
        density: probability of 1 in random matrices
        rng: numpy random generator

    Returns:
        CSSCode or None if construction fails
    """
    if rng is None:
        rng = np.random.default_rng()

    if target_k is None:
        target_k = 1

    # Target ranks: rx + rz = n - k
    total_checks = n - target_k
    rx = total_checks // 2
    rz = total_checks - rx

    if rx <= 0 or rz <= 0:
        return None

    # Generate random Hx
    Hx = rng.random((rx, n)) < density
    Hx = Hx.astype(np.uint8)

    # Remove zero rows
    Hx = Hx[np.any(Hx, axis=1)]
    if Hx.shape[0] == 0:
        return None

    # Hz must satisfy: Hz @ Hx.T = 0 (mod 2)
    # This means each row of Hz is in the nullspace of Hx.T
    # i.e., Hz rows are in the row nullspace of Hx over GF(2)
    null_hx = gf2_nullspace(Hx)
    if null_hx.shape[0] == 0:
        return None

    # Select rz random rows from the nullspace
    if null_hx.shape[0] < rz:
        rz = null_hx.shape[0]

    # Take random combinations of nullspace vectors
    Hz = np.zeros((rz, n), dtype=np.uint8)
    for i in range(rz):
        # Random combination of nullspace vectors
        coeffs = rng.integers(0, 2, size=null_hx.shape[0])
        Hz[i] = (coeffs @ null_hx) % 2

    # Remove zero rows
    Hz = Hz[np.any(Hz, axis=1)]
    if Hz.shape[0] == 0:
        return None

    code = CSSCode(Hx=Hx, Hz=Hz)
    valid, _ = code.is_valid()
    if not valid:
        return None

    # Check we got at least 1 logical qubit
    if code.k < 1:
        return None

    return code


def random_css_from_classical(n, rng=None):
    """Generate a CSS code from a random classical linear code.

    Uses a random parity check matrix H for both Hx and Hz.
    This automatically satisfies commutativity since H @ H.T = 0 (mod 2)
    for self-orthogonal codes.

    For non-self-orthogonal H, uses H for Hx and builds Hz from
    the intersection of rowspace(H) and nullspace(H).
    """
    if rng is None:
        rng = np.random.default_rng()

    # Generate random parity check matrix
    r = rng.integers(2, n - 1)
    H = rng.integers(0, 2, size=(r, n), dtype=np.uint8)
    H = H[np.any(H, axis=1)]  # Remove zero rows
    if H.shape[0] < 2:
        return None

    # Check if H is self-orthogonal: H @ H.T = 0 (mod 2)
    product = (H @ H.T) % 2
    if not np.any(product):
        # Self-orthogonal — use H for both
        code = CSSCode(Hx=H.copy(), Hz=H.copy())
        if code.k >= 1:
            return code

    # Not self-orthogonal — use H for Hx, build Hz from nullspace
    null_H = gf2_nullspace(H)
    if null_H.shape[0] == 0:
        return None

    # Pick subset of nullspace rows
    num_hz = min(null_H.shape[0], n // 3)
    if num_hz == 0:
        return None
    idx = rng.choice(null_H.shape[0], size=num_hz, replace=False)
    Hz = null_H[idx]
    Hz = Hz[np.any(Hz, axis=1)]
    if Hz.shape[0] == 0:
        return None

    code = CSSCode(Hx=H, Hz=Hz)
    valid, _ = code.is_valid()
    if not valid:
        return None
    if code.k < 1:
        return None
    return code


def random_self_orthogonal(n, target_k=1, rng=None):
    """Generate a random self-orthogonal code over GF(2).

    A self-orthogonal code has H @ H.T = 0 (mod 2), meaning every pair
    of rows has even overlap. Using H for both Hx and Hz automatically
    satisfies CSS commutativity and tends to produce higher-distance codes.

    Strategy: build rows one at a time, each constrained to be orthogonal
    to all previous rows.
    """
    if rng is None:
        rng = np.random.default_rng()

    # Number of check rows: (n - k) / 2 since same H used for both
    num_rows = (n - target_k) // 2
    if num_rows < 1:
        return None

    rows = []
    for _ in range(num_rows):
        # Try to find a row orthogonal to all existing rows
        for _ in range(100):  # attempts
            row = rng.integers(0, 2, size=n).astype(np.uint8)
            if not np.any(row):
                continue

            # Check self-orthogonality: even weight
            if np.sum(row) % 2 != 0:
                continue

            # Check orthogonality with existing rows
            ok = True
            for existing in rows:
                if np.sum(row & existing) % 2 != 0:
                    ok = False
                    break
            if ok:
                rows.append(row)
                break

    if len(rows) < 1:
        return None

    H = np.array(rows, dtype=np.uint8)

    # Verify self-orthogonality
    product = (H @ H.T) % 2
    if np.any(product):
        return None

    code = CSSCode(Hx=H.copy(), Hz=H.copy())
    if code.k < 1:
        return None
    return code


def random_bicycle_code(n, rng=None):
    """Generate a random bicycle (circulant) CSS code.

    Bicycle codes use circulant matrices — defined by a single row that
    is cyclically shifted. They have algebraic structure that often
    produces good distance.

    For a bicycle code on n qubits (n must be even):
    - m = n/2
    - Hx = [A | B] where A, B are m x m circulant matrices
    - Hz = [B.T | A.T]
    - Commutativity: AB.T + BA.T = 0 (mod 2), which holds when
      A and B commute (always true for circulants over GF(2))
    """
    if rng is None:
        rng = np.random.default_rng()

    if n % 2 != 0 or n < 6:
        return None

    m = n // 2

    # Generate random first rows for two circulant matrices
    # Use moderate weight (not too sparse, not too dense)
    weight = max(2, m // 3)

    a_row = np.zeros(m, dtype=np.uint8)
    a_positions = rng.choice(m, size=weight, replace=False)
    a_row[a_positions] = 1

    b_row = np.zeros(m, dtype=np.uint8)
    b_positions = rng.choice(m, size=weight, replace=False)
    b_row[b_positions] = 1

    # Build circulant matrices
    A = np.zeros((m, m), dtype=np.uint8)
    B = np.zeros((m, m), dtype=np.uint8)
    for i in range(m):
        A[i] = np.roll(a_row, i)
        B[i] = np.roll(b_row, i)

    # Hx = [A | B], Hz = [B.T | A.T]
    Hx = np.hstack([A, B]).astype(np.uint8)
    Hz = np.hstack([B.T, A.T]).astype(np.uint8)

    # Remove zero rows
    Hx = Hx[np.any(Hx, axis=1)]
    Hz = Hz[np.any(Hz, axis=1)]

    if Hx.shape[0] == 0 or Hz.shape[0] == 0:
        return None

    code = CSSCode(Hx=Hx, Hz=Hz)
    valid, _ = code.is_valid()
    if not valid:
        return None
    if code.k < 1:
        return None
    return code


def punctured_code(code, positions=None, rng=None):
    """Create a new code by removing qubits (columns) from an existing code.

    Puncturing can increase code rate at the cost of distance.
    Useful for exploring the neighborhood of known good codes.
    """
    if rng is None:
        rng = np.random.default_rng()

    n = code.n
    if n <= 5:
        return None

    if positions is None:
        # Remove 1-2 random qubits
        num_remove = rng.integers(1, min(3, n - 4))
        positions = rng.choice(n, size=num_remove, replace=False)

    keep = np.array([i for i in range(n) if i not in positions])
    if len(keep) < 4:
        return None

    Hx = code.Hx[:, keep]
    Hz = code.Hz[:, keep]

    # Remove zero rows
    Hx = Hx[np.any(Hx, axis=1)]
    Hz = Hz[np.any(Hz, axis=1)]

    if Hx.shape[0] == 0 or Hz.shape[0] == 0:
        return None

    new_code = CSSCode(Hx=Hx, Hz=Hz)
    valid, _ = new_code.is_valid()
    if not valid:
        return None
    if new_code.k < 1:
        return None
    return new_code


def hypergraph_product(H1, H2):
    """Hypergraph product of two classical codes.

    Given classical parity check matrices H1 (r1 x n1) and H2 (r2 x n2),
    produces a CSS code with:
        n = n1*n2 + r1*r2
        Hx = [H1 ⊗ I_n2  |  I_r1 ⊗ H2.T]
        Hz = [I_n1 ⊗ H2   |  H1.T ⊗ I_r2]

    The commutativity Hx @ Hz.T = 0 follows from H1 @ H1.T over GF(2)... no,
    it follows algebraically from the product structure.
    """
    H1 = np.asarray(H1, dtype=np.uint8)
    H2 = np.asarray(H2, dtype=np.uint8)

    r1, n1 = H1.shape
    r2, n2 = H2.shape

    # Hx = [H1 ⊗ I_n2 | I_r1 ⊗ H2.T]
    Hx_left = np.kron(H1, np.eye(n2, dtype=np.uint8))  # (r1*n2, n1*n2)
    Hx_right = np.kron(np.eye(r1, dtype=np.uint8), H2.T)  # (r1*n2, r1*r2)
    Hx = np.hstack([Hx_left, Hx_right]) % 2

    # Hz = [I_n1 ⊗ H2 | H1.T ⊗ I_r2]
    Hz_left = np.kron(np.eye(n1, dtype=np.uint8), H2)  # (n1*r2, n1*n2)
    Hz_right = np.kron(H1.T, np.eye(r2, dtype=np.uint8))  # (n1*r2, r1*r2)
    Hz = np.hstack([Hz_left, Hz_right]) % 2

    Hx = Hx.astype(np.uint8)
    Hz = Hz.astype(np.uint8)

    # Remove zero rows
    Hx = Hx[np.any(Hx, axis=1)]
    Hz = Hz[np.any(Hz, axis=1)]

    return CSSCode(Hx=Hx, Hz=Hz)
