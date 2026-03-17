"""
Validation of stabilizer codes and reference code recovery.

Validates:
1. CSS commutativity constraint
2. Code parameters match expected values
3. Known reference codes produce correct distances
"""

import numpy as np
from codes import CSSCode, gf2_rank
from distance import exact_distance_css
from known_codes import get_reference_codes
import time
import sys


def validate_css_properties(code: CSSCode):
    """Run all property checks on a CSS code.

    Returns dict with check results.
    """
    results = {}

    # 1. Commutativity
    valid, reason = code.is_valid()
    results["commutativity"] = {"pass": valid, "reason": reason}

    # 2. Dimensions
    n = code.n
    k = code.k
    rx = gf2_rank(code.Hx)
    rz = gf2_rank(code.Hz)
    results["dimensions"] = {
        "n": n,
        "k": k,
        "rx": rx,
        "rz": rz,
        "check": f"k = n - rx - rz = {n} - {rx} - {rz} = {n - rx - rz}",
        "pass": k >= 0,
    }

    # 3. Stabilizer independence (no redundant rows)
    results["hx_independence"] = {
        "rows": int(code.Hx.shape[0]),
        "rank": rx,
        "pass": code.Hx.shape[0] == rx,
        "reason": "All rows independent" if code.Hx.shape[0] == rx else f"{code.Hx.shape[0] - rx} redundant rows",
    }
    results["hz_independence"] = {
        "rows": int(code.Hz.shape[0]),
        "rank": rz,
        "pass": code.Hz.shape[0] == rz,
        "reason": "All rows independent" if code.Hz.shape[0] == rz else f"{code.Hz.shape[0] - rz} redundant rows",
    }

    # 4. Stabilizer weights
    results["weights"] = {
        "max": int(code.max_stabilizer_weight),
        "avg": round(code.avg_stabilizer_weight, 2),
    }

    return results


def validate_reference_codes(verbose=True):
    """Validate all reference codes — compute distance, check against expected.

    This is the equivalent of exoplanet's known_planet_recovery.
    """
    refs = get_reference_codes()
    results = []

    if verbose:
        print("=" * 70)
        print("REFERENCE CODE VALIDATION")
        print("=" * 70)

    total = len(refs)
    passed = 0

    for name, Hx, Hz, expected_n, expected_k, expected_d in refs:
        if verbose:
            print(f"\n--- {name} ---")

        code = CSSCode(Hx=Hx, Hz=Hz)
        t0 = time.time()

        # Check basic properties
        props = validate_css_properties(code)
        valid = props["commutativity"]["pass"]
        actual_n = code.n
        actual_k = code.k

        if verbose:
            status = "PASS" if valid else "FAIL"
            print(f"  Commutativity: {status}")
            print(f"  n={actual_n} (expected {expected_n}): {'PASS' if actual_n == expected_n else 'FAIL'}")
            print(f"  k={actual_k} (expected {expected_k}): {'PASS' if actual_k == expected_k else 'FAIL'}")

        # Compute distance
        d = exact_distance_css(code, timeout_sec=60)
        elapsed = time.time() - t0

        if d is not None:
            code.distance = d
            d_pass = (d == expected_d)
            if verbose:
                status = "PASS" if d_pass else "FAIL"
                print(f"  d={d} (expected {expected_d}): {status}")
                print(f"  Time: {elapsed:.2f}s")
        else:
            d_pass = False
            if verbose:
                print(f"  d=TIMEOUT (expected {expected_d}): FAIL")

        all_pass = (valid and actual_n == expected_n and actual_k == expected_k and d_pass)
        if all_pass:
            passed += 1

        results.append({
            "name": name,
            "expected": [expected_n, expected_k, expected_d],
            "actual": [actual_n, actual_k, d],
            "pass": all_pass,
            "time_sec": round(elapsed, 3),
        })

        if verbose:
            print(f"  Overall: {'PASS' if all_pass else 'FAIL'}")

    if verbose:
        print(f"\n{'=' * 70}")
        print(f"RESULTS: {passed}/{total} passed")
        print(f"{'=' * 70}")

    return results


if __name__ == "__main__":
    results = validate_reference_codes(verbose=True)
    all_passed = all(r["pass"] for r in results)
    sys.exit(0 if all_passed else 1)
