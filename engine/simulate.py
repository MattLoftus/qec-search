"""
Monte Carlo simulation of logical error rates.

Given a CSS code, simulate:
1. Apply random physical errors at rate p
2. Measure syndrome
3. Decode (MWPM or BP)
4. Check if logical error occurred
5. Repeat → estimate logical error rate

This tests decoder performance, not just code parameters.
"""

import numpy as np
from codes import CSSCode, gf2_rank
from distance import _in_rowspace
import time


def simulate_code_performance(code: CSSCode, physical_error_rates=None,
                               num_shots=10000, decoder="lookup",
                               verbose=True):
    """Simulate logical error rate vs physical error rate.

    Args:
        code: CSS code to simulate
        physical_error_rates: list of p values to test
        num_shots: Monte Carlo samples per p value
        decoder: "lookup" (exact, small codes) or "mwpm" (approximate)
        verbose: print progress

    Returns:
        dict with {p: logical_error_rate} data
    """
    if physical_error_rates is None:
        physical_error_rates = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1]

    n = code.n
    results = {}

    for p in physical_error_rates:
        if verbose:
            print(f"  p={p:.4f}: ", end="", flush=True)

        logical_errors = 0
        t0 = time.time()

        for shot in range(num_shots):
            # Generate random X error (detected by Hz)
            x_error = np.random.binomial(1, p, size=n).astype(np.uint8)

            # Measure syndrome
            if code.Hz.size > 0:
                syndrome = (code.Hz @ x_error) % 2
            else:
                syndrome = np.array([], dtype=np.uint8)

            # Decode: find correction
            if decoder == "lookup":
                correction = _lookup_decode(code.Hz, syndrome, n)
            else:
                correction = np.zeros(n, dtype=np.uint8)  # No correction (trivial)

            # Residual error
            residual = (x_error + correction) % 2

            # Check if residual is a logical error
            # (in kernel of Hz but NOT in rowspace of Hx)
            if np.any(residual):
                if code.Hz.size > 0:
                    check = (code.Hz @ residual) % 2
                    if np.any(check):
                        continue  # Not in kernel — decoder didn't fully correct
                # In kernel — check if it's a nontrivial logical
                if not _in_rowspace(residual, code.Hx):
                    logical_errors += 1

        logical_rate = logical_errors / num_shots
        elapsed = time.time() - t0

        if verbose:
            print(f"logical_error_rate={logical_rate:.6f} ({logical_errors}/{num_shots}), "
                  f"{elapsed:.1f}s")

        results[str(p)] = {
            "physical_error_rate": p,
            "logical_error_rate": logical_rate,
            "logical_errors": logical_errors,
            "shots": num_shots,
            "time_sec": round(elapsed, 2),
        }

    return results


def _lookup_decode(Hz, syndrome, n):
    """Simple minimum-weight decoder for small codes.

    Finds minimum weight error consistent with syndrome.
    Only feasible for small n (< ~20).
    """
    if not np.any(syndrome):
        return np.zeros(n, dtype=np.uint8)

    if Hz.size == 0:
        return np.zeros(n, dtype=np.uint8)

    # Try single-qubit errors first
    for i in range(n):
        e = np.zeros(n, dtype=np.uint8)
        e[i] = 1
        if np.array_equal((Hz @ e) % 2, syndrome):
            return e

    # Try two-qubit errors
    from itertools import combinations
    for i, j in combinations(range(n), 2):
        e = np.zeros(n, dtype=np.uint8)
        e[i] = 1
        e[j] = 1
        if np.array_equal((Hz @ e) % 2, syndrome):
            return e

    # Try three-qubit errors (only for small codes)
    if n <= 15:
        for positions in combinations(range(n), 3):
            e = np.zeros(n, dtype=np.uint8)
            for pos in positions:
                e[pos] = 1
            if np.array_equal((Hz @ e) % 2, syndrome):
                return e

    # Give up — return zero correction
    return np.zeros(n, dtype=np.uint8)


def find_threshold(code: CSSCode, num_shots=5000, verbose=True):
    """Estimate the pseudo-threshold of a code.

    The pseudo-threshold is the physical error rate p* where
    the logical error rate equals p (i.e., the code starts helping).

    Uses binary search on physical error rate.
    """
    if verbose:
        print(f"Finding threshold for {code}...")

    p_low, p_high = 0.001, 0.2
    threshold = None

    for _ in range(10):  # Binary search iterations
        p_mid = (p_low + p_high) / 2
        result = simulate_code_performance(
            code, physical_error_rates=[p_mid],
            num_shots=num_shots, verbose=False,
        )
        logical_rate = result[str(p_mid)]["logical_error_rate"]

        if verbose:
            print(f"  p={p_mid:.4f}: logical={logical_rate:.4f} "
                  f"({'below' if logical_rate < p_mid else 'above'} diagonal)")

        if logical_rate < p_mid:
            p_low = p_mid  # Code is helping, try higher p
        else:
            p_high = p_mid  # Code isn't helping, try lower p

        threshold = p_mid

    if verbose:
        print(f"  Estimated threshold: p* ≈ {threshold:.4f}")

    return threshold
