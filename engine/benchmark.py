"""
Benchmark search results against known code tables (codetables.de).
"""

import json
import os
from known_codes import KNOWN_BOUNDS


def benchmark_results(best_codes, verbose=True):
    """Compare found codes against known bounds.

    Args:
        best_codes: list of SearchResult or dicts with n, k, d fields

    Returns:
        dict with summary and per-code comparisons
    """
    comparisons = []
    matches = 0
    improvements = 0
    below = 0
    no_reference = 0

    for item in best_codes:
        if hasattr(item, 'n'):
            n, k, d = item.n, item.k, item.d
        else:
            n, k, d = item['n'], item['k'], item['d']

        key = (n, k)
        if key in KNOWN_BOUNDS:
            known_lower, known_upper = KNOWN_BOUNDS[key]

            if d > known_lower:
                status = "IMPROVEMENT"
                improvements += 1
            elif d == known_lower:
                status = "MATCHES_BEST"
                matches += 1
            else:
                status = "BELOW_BEST"
                below += 1

            comp = {
                "n": n, "k": k, "d_found": d,
                "d_known_lower": known_lower,
                "d_known_upper": known_upper,
                "status": status,
                "gap": d - known_lower,
            }
        else:
            status = "NO_REFERENCE"
            no_reference += 1
            comp = {
                "n": n, "k": k, "d_found": d,
                "d_known_lower": None,
                "d_known_upper": None,
                "status": status,
                "gap": None,
            }

        comparisons.append(comp)

        if verbose:
            if comp["d_known_lower"] is not None:
                print(f"  [[{n}, {k}, {d}]] vs known [{comp['d_known_lower']}, {comp['d_known_upper']}]: "
                      f"{status}" + (f" (+{comp['gap']})" if comp['gap'] and comp['gap'] > 0 else ""))
            else:
                print(f"  [[{n}, {k}, {d}]]: no reference data")

    summary = {
        "total": len(best_codes),
        "matches_best": matches,
        "improvements": improvements,
        "below_best": below,
        "no_reference": no_reference,
    }

    if verbose:
        print(f"\nSummary: {matches} match, {improvements} improve, "
              f"{below} below, {no_reference} no reference")

    return {"summary": summary, "comparisons": comparisons}


def save_benchmark(benchmark_result, path="results/benchmark.json"):
    """Save benchmark results to JSON."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        json.dump(benchmark_result, f, indent=2)


def load_benchmark(path="results/benchmark.json"):
    """Load benchmark results from JSON."""
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)
