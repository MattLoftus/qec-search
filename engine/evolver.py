#!/usr/bin/env python3
"""
QEC Code Search — Main CLI entry point.

Usage:
    python3 evolver.py --validate                     # Validate reference codes
    python3 evolver.py --mode random --n 15 --k 1     # Random search
    python3 evolver.py --mode genetic --n 15 --k 1    # Genetic algorithm
    python3 evolver.py --mode algebraic               # Hypergraph product search
    python3 evolver.py --benchmark                    # Compare results vs known
    python3 evolver.py --simulate --n 7 --k 1         # Simulate error rates
    python3 evolver.py --status                       # Show search status
"""

import argparse
import json
import os
import sys
import time

# Add engine dir to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from codes import CSSCode
from validate import validate_reference_codes
from search import random_search, genetic_search, algebraic_search, SearchResult
from benchmark import benchmark_results, save_benchmark
from simulate import simulate_code_performance, find_threshold
from known_codes import get_reference_codes, steane_code, KNOWN_BOUNDS


RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
CAMPAIGNS_DIR = os.path.join(RESULTS_DIR, "campaigns")
BEST_CODES_PATH = os.path.join(RESULTS_DIR, "best_codes.json")


def ensure_dirs():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(CAMPAIGNS_DIR, exist_ok=True)


def save_campaign(stats, name=None):
    """Save search campaign results."""
    ensure_dirs()

    if name is None:
        name = f"{stats.strategy}_{int(time.time())}"

    campaign = {
        "name": name,
        "strategy": stats.strategy,
        "total_tested": stats.total_codes_tested,
        "valid_codes": stats.valid_codes,
        "codes_with_k": stats.codes_with_k_ge_1,
        "generations": stats.generations_completed,
        "elapsed_sec": round(stats.elapsed_sec, 2),
        "best_distance": {f"{n},{k}": d for (n, k), d in stats.best_distance.items()},
        "best_codes": [r.to_dict() for r in stats.best_codes],
        "timestamp": time.time(),
    }

    path = os.path.join(CAMPAIGNS_DIR, f"{name}.json")
    with open(path, 'w') as f:
        json.dump(campaign, f, indent=2)
    print(f"\nCampaign saved: {path}")

    # Update all-time best codes
    _update_best_codes(stats.best_codes)

    return path


def _update_best_codes(new_codes):
    """Merge new codes into all-time best codes file."""
    ensure_dirs()

    existing = []
    if os.path.exists(BEST_CODES_PATH):
        with open(BEST_CODES_PATH) as f:
            existing = json.load(f)

    # Index existing by (n, k)
    best_by_nk = {}
    for entry in existing:
        key = (entry["n"], entry["k"])
        if key not in best_by_nk or entry["d"] > best_by_nk[key]["d"]:
            best_by_nk[key] = entry

    # Merge new
    updated = False
    for result in new_codes:
        entry = result.to_dict()
        key = (entry["n"], entry["k"])
        if key not in best_by_nk or entry["d"] > best_by_nk[key]["d"]:
            best_by_nk[key] = entry
            updated = True

    if updated:
        all_best = sorted(best_by_nk.values(), key=lambda x: (x["n"], x["k"]))
        with open(BEST_CODES_PATH, 'w') as f:
            json.dump(all_best, f, indent=2)
        print(f"Updated best codes: {BEST_CODES_PATH}")


def cmd_validate(args):
    """Validate reference codes — the first thing to run."""
    results = validate_reference_codes(verbose=True)
    all_passed = all(r["pass"] for r in results)
    return 0 if all_passed else 1


def cmd_search(args):
    """Run a search campaign."""
    n = args.n
    k = args.k

    print(f"Starting {args.mode} search: n={n}, k={k}")
    print(f"Known bounds for [[{n}, {k}, ?]]: ", end="")

    key = (n, k)
    if key in KNOWN_BOUNDS:
        lower, upper = KNOWN_BOUNDS[key]
        print(f"d in [{lower}, {upper}]")
    else:
        print("no reference data")

    def progress_callback(*cb_args):
        pass  # Search functions print their own progress

    if args.mode == "random":
        stats = random_search(
            n=n, target_k=k,
            num_samples=args.samples or 10000,
            distance_timeout=args.timeout or 30,
            callback=progress_callback,
        )
    elif args.mode == "genetic":
        stats = genetic_search(
            n=n, target_k=k,
            pop_size=args.pop or 100,
            num_generations=args.gens or 200,
            mutation_rate=args.mutation_rate or 0.05,
            distance_timeout=args.timeout or 30,
            callback=progress_callback,
        )
    elif args.mode == "algebraic":
        stats = algebraic_search(
            max_n=args.n_max or 30,
            callback=progress_callback,
        )
    else:
        print(f"Unknown mode: {args.mode}")
        return 1

    # Print summary
    print(f"\n{'=' * 60}")
    print(f"SEARCH COMPLETE: {args.mode}")
    print(f"{'=' * 60}")
    print(f"  Codes tested: {stats.total_codes_tested}")
    print(f"  Valid CSS: {stats.valid_codes}")
    print(f"  With k >= 1: {stats.codes_with_k_ge_1}")
    print(f"  Time: {stats.elapsed_sec:.1f}s")

    if stats.best_distance:
        print(f"\nBest distances found:")
        for (bn, bk), bd in sorted(stats.best_distance.items()):
            status = ""
            if (bn, bk) in KNOWN_BOUNDS:
                kl, ku = KNOWN_BOUNDS[(bn, bk)]
                if bd > kl:
                    status = " *** BEATS KNOWN ***"
                elif bd == kl:
                    status = " (matches best known)"
                else:
                    status = f" (known best: {kl})"
            print(f"  [[{bn}, {bk}, {bd}]]{status}")

    # Benchmark
    if stats.best_codes:
        print(f"\nBenchmark vs codetables.de:")
        benchmark = benchmark_results(stats.best_codes, verbose=True)
        save_benchmark(benchmark, os.path.join(RESULTS_DIR, "benchmark.json"))

    # Save campaign
    campaign_name = f"{args.mode}_n{n}_k{k}"
    save_campaign(stats, name=campaign_name)

    return 0


def cmd_benchmark(args):
    """Benchmark all-time best codes against known tables."""
    if not os.path.exists(BEST_CODES_PATH):
        print("No best codes found. Run a search first.")
        return 1

    with open(BEST_CODES_PATH) as f:
        best = json.load(f)

    print(f"Benchmarking {len(best)} best codes against codetables.de:")
    benchmark = benchmark_results(best, verbose=True)
    save_benchmark(benchmark, os.path.join(RESULTS_DIR, "benchmark.json"))
    return 0


def cmd_simulate(args):
    """Simulate error rates for a code."""
    if args.code:
        with open(args.code) as f:
            data = json.load(f)
        if isinstance(data, list):
            entry = data[args.index or 0]
        else:
            entry = data
        code = CSSCode.from_dict(entry)
    else:
        # Default: Steane code
        Hx, Hz = steane_code()
        code = CSSCode(Hx=Hx, Hz=Hz)

    print(f"Simulating {code}...")
    results = simulate_code_performance(
        code,
        num_shots=args.shots or 10000,
        verbose=True,
    )

    # Save results
    ensure_dirs()
    path = os.path.join(RESULTS_DIR, f"simulation_{code.code_hash()}.json")
    with open(path, 'w') as f:
        json.dump({"code": code.to_dict(), "results": results}, f, indent=2)
    print(f"\nResults saved: {path}")

    if args.threshold:
        threshold = find_threshold(code, verbose=True)

    return 0


def cmd_status(args):
    """Show search status."""
    ensure_dirs()

    print("QEC Code Search Status")
    print("=" * 50)

    # Campaigns
    campaigns = []
    if os.path.exists(CAMPAIGNS_DIR):
        for f in sorted(os.listdir(CAMPAIGNS_DIR)):
            if f.endswith('.json'):
                with open(os.path.join(CAMPAIGNS_DIR, f)) as fh:
                    campaigns.append(json.load(fh))

    if campaigns:
        print(f"\nCampaigns: {len(campaigns)}")
        for c in campaigns:
            print(f"  {c['name']}: {c['total_tested']} tested, "
                  f"{len(c.get('best_codes', []))} best codes, "
                  f"{c['elapsed_sec']:.0f}s")
    else:
        print("\nNo campaigns yet. Run: python3 evolver.py --mode random --n 10 --k 1")

    # Best codes
    if os.path.exists(BEST_CODES_PATH):
        with open(BEST_CODES_PATH) as f:
            best = json.load(f)
        print(f"\nAll-time best codes: {len(best)}")
        for entry in best:
            d = entry.get('d', '?')
            beats = " ***" if entry.get('beats_known') else ""
            print(f"  [[{entry['n']}, {entry['k']}, {d}]]{beats}")

    # Known bounds coverage
    print(f"\nKnown bounds loaded: {len(KNOWN_BOUNDS)} parameter sets")
    gaps = [(n, k, dl, du) for (n, k), (dl, du) in KNOWN_BOUNDS.items() if dl < du]
    print(f"Open gaps (d_lower < d_upper): {len(gaps)}")
    if gaps:
        for n, k, dl, du in sorted(gaps)[:10]:
            print(f"  [[{n}, {k}, ?]]: d in [{dl}, {du}] — gap of {du - dl}")

    return 0


def main():
    parser = argparse.ArgumentParser(description="QEC Code Search")
    parser.add_argument("--validate", action="store_true", help="Validate reference codes")
    parser.add_argument("--mode", choices=["random", "genetic", "algebraic"],
                        help="Search mode")
    parser.add_argument("--n", type=int, default=10, help="Number of physical qubits")
    parser.add_argument("--k", type=int, default=1, help="Number of logical qubits")
    parser.add_argument("--n-max", type=int, default=30, help="Max n for algebraic search")
    parser.add_argument("--samples", type=int, help="Samples for random search")
    parser.add_argument("--pop", type=int, help="Population for genetic search")
    parser.add_argument("--gens", type=int, help="Generations for genetic search")
    parser.add_argument("--mutation-rate", type=float, help="Mutation rate")
    parser.add_argument("--timeout", type=int, help="Distance computation timeout (sec)")
    parser.add_argument("--benchmark", action="store_true", help="Benchmark results")
    parser.add_argument("--simulate", action="store_true", help="Simulate error rates")
    parser.add_argument("--code", type=str, help="Code JSON file for simulation")
    parser.add_argument("--index", type=int, help="Index in code JSON array")
    parser.add_argument("--shots", type=int, help="Monte Carlo shots")
    parser.add_argument("--threshold", action="store_true", help="Find pseudo-threshold")
    parser.add_argument("--status", action="store_true", help="Show status")

    args = parser.parse_args()

    if args.validate:
        return cmd_validate(args)
    elif args.mode:
        return cmd_search(args)
    elif args.benchmark:
        return cmd_benchmark(args)
    elif args.simulate:
        return cmd_simulate(args)
    elif args.status:
        return cmd_status(args)
    else:
        parser.print_help()
        return 0


if __name__ == "__main__":
    sys.exit(main())
