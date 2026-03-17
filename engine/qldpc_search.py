"""
Search for codes using the qLDPC library's algebraic constructions.

qLDPC provides:
- Bivariate bicycle codes (BBCode) — Google/IBM's breakthrough construction
- Hypergraph product codes (HGPCode) — from pairs of classical codes
- Lifted product codes (LPCode) — generalization of HGP
- Exact distance computation via get_distance()

This module wraps qLDPC to systematically explore code families and
benchmark against known bounds.

Requires: Python 3.12 venv with qLDPC installed
  Run with: ~/workspace/qec-search/venv/bin/python3 qldpc_search.py
"""

import sys
import os
import json
import time
import itertools
import signal
import numpy as np

# Flush prints immediately for real-time output
import functools
print = functools.partial(print, flush=True)

# qLDPC imports
from qldpc.codes import (
    ClassicalCode, CSSCode as QCSSCode,
    HGPCode, BBCode,
)
import sympy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from known_codes import KNOWN_BOUNDS


def get_distance_with_timeout(code, timeout_sec=60):
    """Compute code distance with a timeout.

    Uses SIGALRM on macOS/Linux. Falls back to get_distance_bound() if
    exact computation times out.
    """
    class TimeoutError(Exception):
        pass

    def handler(signum, frame):
        raise TimeoutError()

    # Try exact distance with timeout
    old_handler = signal.signal(signal.SIGALRM, handler)
    signal.alarm(timeout_sec)
    try:
        d = code.get_distance()
        signal.alarm(0)
        return d, "exact"
    except TimeoutError:
        signal.alarm(0)
        # Try upper bound (faster)
        try:
            d = code.get_distance_bound()
            return d, "bound"
        except Exception:
            return None, "timeout"
    except Exception:
        signal.alarm(0)
        return None, "error"
    finally:
        signal.signal(signal.SIGALRM, old_handler)

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
BEST_CODES_PATH = os.path.join(RESULTS_DIR, "best_codes.json")


def load_best_codes():
    if os.path.exists(BEST_CODES_PATH):
        with open(BEST_CODES_PATH) as f:
            return json.load(f)
    return []


def save_best_codes(codes):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(BEST_CODES_PATH, 'w') as f:
        json.dump(codes, f, indent=2)


def update_best(existing, n, k, d, strategy):
    """Update best codes list if this code improves on existing."""
    by_nk = {}
    for e in existing:
        key = (e['n'], e['k'])
        if key not in by_nk or (e.get('d') or 0) > (by_nk[key].get('d') or 0):
            by_nk[key] = e

    key = (n, k)
    if key not in by_nk or d > (by_nk[key].get('d') or 0):
        by_nk[key] = {
            'type': 'css', 'n': n, 'k': k, 'd': d,
            'strategy': strategy,
            'Hx': [], 'Hz': [],  # Omit matrices for now (large)
            'hash': f'{strategy}_{n}_{k}_{d}',
            'max_weight': 0, 'avg_weight': 0,
            'rate': k / n,
        }
        return True, list(sorted(by_nk.values(), key=lambda x: (x['n'], x['k'])))
    return False, existing


def check_bounds(n, k, d):
    """Check code against known bounds."""
    key = (n, k)
    if key in KNOWN_BOUNDS:
        dl, du = KNOWN_BOUNDS[key]
        if d > dl:
            return "IMPROVEMENT", dl, du
        elif d == dl:
            return "MATCHES", dl, du
        else:
            return "BELOW", dl, du
    return "NO_REF", None, None


# ---------------------------------------------------------------------------
# HGP search
# ---------------------------------------------------------------------------

def build_classical_library():
    """Build a library of classical codes for HGP construction."""
    codes = {}

    # Repetition codes
    for n in range(3, 12):
        H = np.zeros((n - 1, n), dtype=int)
        for i in range(n - 1):
            H[i, i] = 1
            H[i, i + 1] = 1
        codes[f'rep({n})'] = ClassicalCode(H)

    # Hamming codes
    for r in [2, 3, 4]:
        n = 2**r - 1
        H = np.zeros((r, n), dtype=int)
        for j in range(n):
            col = j + 1
            for i in range(r):
                H[i, j] = (col >> i) & 1
        codes[f'ham(r={r})'] = ClassicalCode(H)

    # Extended Hamming (even-weight subcode)
    for r in [3, 4]:
        n = 2**r
        H = np.zeros((r + 1, n), dtype=int)
        for j in range(n - 1):
            col = j + 1
            for i in range(r):
                H[i, j] = (col >> i) & 1
        # Overall parity check row
        H[r, :] = 1
        codes[f'ext_ham(r={r})'] = ClassicalCode(H)

    # Random codes with varied densities
    rng = np.random.default_rng(12345)
    for trial in range(30):
        n = rng.integers(4, 12)
        r = rng.integers(2, max(3, n - 1))
        H = rng.integers(0, 2, size=(r, n)).astype(int)
        H = H[np.any(H, axis=1)]
        if H.shape[0] >= 2:
            codes[f'rand{trial}({H.shape[0]}x{n})'] = ClassicalCode(H)

    return codes


def search_hgp(max_n=80, verbose=True):
    """Search for codes via hypergraph products using qLDPC."""
    classical = build_classical_library()
    names = list(classical.keys())

    if verbose:
        print(f"HGP search: {len(names)} classical codes, max n={max_n}")

    best = load_best_codes()
    results = []
    improvements = 0
    matches = 0

    for i, name1 in enumerate(names):
        for j, name2 in enumerate(names):
            if j < i:
                continue

            try:
                code = HGPCode(classical[name1], classical[name2])
                n, k = code.num_qubits, code.dimension
                if k < 1 or n > max_n:
                    continue

                t0 = time.time()
                timeout = 30 if n <= 40 else 10
                d, method = get_distance_with_timeout(code, timeout_sec=timeout)
                elapsed = time.time() - t0

                if d is None or d < 2:
                    continue

                status, dl, du = check_bounds(n, k, d)
                updated, best = update_best(best, n, k, d, f"hgp({name1},{name2})")

                label = ""
                if status == "IMPROVEMENT":
                    label = f" *** IMPROVEMENT (known: {dl}) ***"
                    improvements += 1
                elif status == "MATCHES":
                    label = f" (matches best: {dl})"
                    matches += 1

                method_tag = f" [{method}]" if method != "exact" else ""
                if updated or status in ("IMPROVEMENT", "MATCHES"):
                    if verbose:
                        print(f"  [[{n}, {k}, {d}]] HGP({name1},{name2}){label}{method_tag} [{elapsed:.1f}s]")

                results.append((n, k, d, f"hgp({name1},{name2})"))

            except Exception:
                pass

    save_best_codes(best)

    if verbose:
        print(f"\nHGP: {len(results)} codes, {matches} matches, {improvements} improvements")

    return results


# ---------------------------------------------------------------------------
# Bivariate bicycle code search
# ---------------------------------------------------------------------------

def search_bb(max_n=80, verbose=True):
    """Search for bivariate bicycle codes.

    BB codes are defined by:
    - Group orders [l, m] (code acts on Z_l x Z_m)
    - Two polynomials a(x,y) and b(x,y) in the group ring

    Systematically tries combinations of orders and polynomial forms.
    """
    x, y = sympy.symbols('x y')

    # Polynomial templates to try
    def poly_templates(l, m):
        """Generate polynomial candidates for given group orders."""
        polys = [
            1 + x,
            1 + y,
            1 + x + y,
            1 + x * y,
            x + y,
            1 + x + x**2,
            1 + y + y**2,
            x + y + x * y,
            1 + x**2,
            1 + y**2,
        ]
        # Add higher powers that make sense for the group order
        if l >= 4:
            polys.extend([1 + x + x**3, x + x**2])
        if m >= 4:
            polys.extend([1 + y + y**3, y + y**2])
        if l >= 3 and m >= 3:
            polys.extend([
                1 + x + y + x * y,
                x + y**2,
                x**2 + y,
            ])
        return polys

    if verbose:
        print(f"BB code search: max n={max_n}")

    best = load_best_codes()
    results = []
    improvements = 0
    matches = 0

    for l in range(3, 13):
        for m in range(l, 13):
            if 2 * l * m > max_n:
                continue

            polys = poly_templates(l, m)

            for pa in polys:
                for pb in polys:
                    if pa == pb:
                        continue

                    try:
                        code = BBCode([l, m], pa, pb)
                        n, k = code.num_qubits, code.dimension
                        if k < 1 or n > max_n:
                            continue

                        t0 = time.time()
                        timeout = 60 if n <= 30 else (30 if n <= 50 else 10)
                        d, method = get_distance_with_timeout(code, timeout_sec=timeout)
                        elapsed = time.time() - t0

                        if d is None or d < 2:
                            continue

                        status, dl, du = check_bounds(n, k, d)
                        updated, best = update_best(best, n, k, d,
                                                     f"bb([{l},{m}],{pa},{pb})")

                        label = ""
                        if status == "IMPROVEMENT":
                            label = f" *** IMPROVEMENT (known: {dl}) ***"
                            improvements += 1
                        elif status == "MATCHES":
                            label = f" (matches best: {dl})"
                            matches += 1

                        method_tag = f" [{method}]" if method != "exact" else ""
                        if updated or status in ("IMPROVEMENT", "MATCHES"):
                            if verbose:
                                print(f"  [[{n}, {k}, {d}]] BB([{l},{m}], {pa}, {pb}){label}{method_tag} [{elapsed:.1f}s]")

                        results.append((n, k, d, f"bb([{l},{m}],{pa},{pb})"))

                    except Exception:
                        pass

    save_best_codes(best)

    if verbose:
        print(f"\nBB: {len(results)} codes, {matches} matches, {improvements} improvements")

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    import argparse
    parser = argparse.ArgumentParser(description="qLDPC-powered code search")
    parser.add_argument("--mode", choices=["hgp", "bb", "all"], default="all")
    parser.add_argument("--max-n", type=int, default=80)
    args = parser.parse_args()

    print(f"QEC Code Search — qLDPC mode ({args.mode})")
    print(f"Max block length: {args.max_n}")
    print("=" * 60)

    if args.mode in ("hgp", "all"):
        print("\n--- Hypergraph Product Search ---")
        hgp_results = search_hgp(max_n=args.max_n)

    if args.mode in ("bb", "all"):
        print("\n--- Bivariate Bicycle Code Search ---")
        bb_results = search_bb(max_n=args.max_n)

    # Final summary
    best = load_best_codes()
    print(f"\n{'=' * 60}")
    print(f"Total best codes: {len(best)}")

    match_count = 0
    improve_count = 0
    for entry in best:
        status, dl, du = check_bounds(entry['n'], entry['k'], entry.get('d', 0))
        if status == "MATCHES":
            match_count += 1
        elif status == "IMPROVEMENT":
            improve_count += 1

    print(f"Matches known bounds: {match_count}")
    print(f"Improvements: {improve_count}")


if __name__ == "__main__":
    main()
