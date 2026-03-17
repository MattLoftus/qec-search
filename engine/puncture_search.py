#!/usr/bin/env python3
"""
Puncture and shorten search: systematically reduce published large codes
to find good smaller codes that beat known bounds.

Two reduction operations:
  Puncturing: remove columns from Hx and Hz. Restore commutativity via
              symplectic projection (kernel of Hx*Hz^T).
  Shortening: fix qubits to 0 -- keep rows with 0 in that column, then
              remove the column. Always preserves commutativity.

Usage:
    ../venv/bin/python3 puncture_search.py --code bravyi72 --trials 200
    ../venv/bin/python3 puncture_search.py --all --trials 200
"""

import argparse
import json
import os
import signal
import sys
import time
from datetime import datetime, timezone

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

from codes import CSSCode, gf2_rank, gf2_nullspace
from known_codes import KNOWN_BOUNDS


# ---------------------------------------------------------------------------
# Nice priority & SIGALRM timeout
# ---------------------------------------------------------------------------
try:
    os.nice(15)
except OSError:
    pass


class AlarmTimeout(Exception):
    pass


def _alarm_handler(signum, frame):
    raise AlarmTimeout("Timed out")


signal.signal(signal.SIGALRM, _alarm_handler)


# ---------------------------------------------------------------------------
# Published codes
# ---------------------------------------------------------------------------

def build_published_codes():
    """Build published BB codes using qLDPC."""
    from qldpc.codes import BBCode
    import sympy
    x, y = sympy.symbols('x y')

    specs = [
        ("bravyi72", [6, 6], x**3 + y + y**2, y**3 + x + x**2, 72, 12, 6),
        ("bravyi108", [9, 6], x**3 + y + y**2, y**3 + x + x**2, 108, 8, 10),
        ("bravyi144", [12, 6], x**3 + y + y**2, y**3 + x + x**2, 144, 12, 12),
        ("wangmueller98", [7, 7], x**3 + y**5 + y**6, y**2 + x**3 + x**5, 98, 6, 12),
    ]

    codes = {}
    for name, orders, poly_a, poly_b, exp_n, exp_k, exp_d in specs:
        print(f"  Building {name}...", end=" ", flush=True)
        t0 = time.time()
        bb = BBCode(orders, poly_a, poly_b)
        Hx = np.array(bb.matrix_x, dtype=np.uint8) % 2
        Hz = np.array(bb.matrix_z, dtype=np.uint8) % 2
        n = Hx.shape[1]
        k = n - gf2_rank(Hx) - gf2_rank(Hz)
        print(f"[[{n},{k},{exp_d}]] in {time.time()-t0:.1f}s", flush=True)
        assert n == exp_n and k == exp_k
        codes[name] = (Hx, Hz, exp_n, exp_k, exp_d)

    return codes


# ---------------------------------------------------------------------------
# Puncturing
# ---------------------------------------------------------------------------

def puncture(Hx, Hz, positions_to_remove):
    """Remove columns, project onto maximal commuting subspace if needed."""
    n = Hx.shape[1]
    keep = np.array(sorted(set(range(n)) - set(positions_to_remove)))
    if len(keep) < 4:
        return None

    Hx_r = Hx[:, keep]
    Hz_r = Hz[:, keep]
    Hx_r = Hx_r[np.any(Hx_r, axis=1)]
    Hz_r = Hz_r[np.any(Hz_r, axis=1)]

    if Hx_r.shape[0] == 0 or Hz_r.shape[0] == 0:
        return None

    S = (Hx_r @ Hz_r.T) % 2
    if not np.any(S):
        return Hx_r, Hz_r

    ker_ST = gf2_nullspace(S.T)
    ker_S = gf2_nullspace(S)
    if ker_ST.shape[0] == 0 or ker_S.shape[0] == 0:
        return None

    Hx_p = (ker_ST @ Hx_r) % 2
    Hz_p = (ker_S @ Hz_r) % 2
    Hx_p = Hx_p[np.any(Hx_p, axis=1)]
    Hz_p = Hz_p[np.any(Hz_p, axis=1)]

    if Hx_p.shape[0] == 0 or Hz_p.shape[0] == 0:
        return None
    if np.any((Hx_p @ Hz_p.T) % 2):
        return None

    return Hx_p, Hz_p


# ---------------------------------------------------------------------------
# Shortening
# ---------------------------------------------------------------------------

def shorten(Hx, Hz, positions_to_fix):
    """Fix qubits to 0. Always preserves commutativity."""
    positions = sorted(positions_to_fix)
    Hx_r = Hx.copy()
    Hz_r = Hz.copy()

    for col in positions:
        Hx_r = Hx_r[Hx_r[:, col] == 0]
        Hz_r = Hz_r[Hz_r[:, col] == 0]
        if Hx_r.shape[0] == 0 or Hz_r.shape[0] == 0:
            return None

    keep = np.array(sorted(set(range(Hx_r.shape[1])) - set(positions)))
    if len(keep) < 4:
        return None

    Hx_r = Hx_r[:, keep]
    Hz_r = Hz_r[:, keep]
    Hx_r = Hx_r[np.any(Hx_r, axis=1)]
    Hz_r = Hz_r[np.any(Hz_r, axis=1)]

    if Hx_r.shape[0] == 0 or Hz_r.shape[0] == 0:
        return None
    if np.any((Hx_r @ Hz_r.T) % 2):
        return None

    return Hx_r, Hz_r


# ---------------------------------------------------------------------------
# Distance computation
# ---------------------------------------------------------------------------

def compute_distance(Hx, Hz, timeout_sec=60):
    """Compute distance. Tries exact (small codes) then graphlike then MC."""
    n = Hx.shape[1]
    k = n - gf2_rank(Hx) - gf2_rank(Hz)
    if k < 1:
        return None, "k<1"

    # Exact distance for small codes
    null_dim = max(n - gf2_rank(Hx), n - gf2_rank(Hz))
    if null_dim <= 20:
        try:
            from distance import exact_distance_css
            code_obj = CSSCode(Hx=Hx, Hz=Hz)
            t = min(timeout_sec, 30)
            signal.alarm(t)
            try:
                d = exact_distance_css(code_obj, timeout_sec=t)
                signal.alarm(0)
                if d is not None:
                    return d, "exact"
            except (AlarmTimeout, Exception):
                pass
            finally:
                signal.alarm(0)
        except Exception:
            signal.alarm(0)

    # Graphlike distance
    try:
        from stim_distance import graphlike_distance
        signal.alarm(20)
        try:
            d, dx, dz = graphlike_distance(Hx, Hz)
            signal.alarm(0)
            if d is not None:
                return d, "graphlike"
        except (AlarmTimeout, Exception):
            pass
        finally:
            signal.alarm(0)
    except Exception:
        signal.alarm(0)

    # MC fallback (only for codes worth investigating)
    if k <= n * 0.4:
        try:
            from stim_distance import monte_carlo_distance_estimate
            signal.alarm(min(timeout_sec, 60))
            try:
                result = monte_carlo_distance_estimate(
                    Hx, Hz, num_shots=10000, verbose=False
                )
                signal.alarm(0)
                d = result.get("estimated_distance")
                if d is not None:
                    return d, "mc"
            except (AlarmTimeout, Exception):
                pass
            finally:
                signal.alarm(0)
        except Exception:
            signal.alarm(0)

    return None, None


# ---------------------------------------------------------------------------
# Bounds check
# ---------------------------------------------------------------------------

def check_against_bounds(n, k, d):
    if d is None:
        return {"status": "unknown"}
    key = (n, k)
    if key in KNOWN_BOUNDS:
        d_lo, d_hi = KNOWN_BOUNDS[key]
        if d > d_lo:
            return {"status": "BEATS_KNOWN", "known_lower": d_lo,
                    "known_upper": d_hi, "improvement": d - d_lo}
        elif d == d_lo:
            return {"status": "matches_best", "known_lower": d_lo,
                    "known_upper": d_hi}
        else:
            return {"status": "below_best", "known_lower": d_lo,
                    "known_upper": d_hi, "gap": d_lo - d}
    else:
        nearby = {}
        for dk in range(-2, 3):
            nk = (n, k + dk)
            if nk in KNOWN_BOUNDS:
                nearby[f"k={k+dk}"] = list(KNOWN_BOUNDS[nk])
        return {"status": "no_bound",
                "note": f"No entry for [[{n},{k},*]]",
                "nearby": nearby}


# ---------------------------------------------------------------------------
# Pattern generators
# ---------------------------------------------------------------------------

def random_puncture_sets(n, num_remove, num_trials, rng):
    seen = set()
    for _ in range(num_trials * 5):
        if len(seen) >= num_trials:
            break
        idx = tuple(sorted(rng.choice(n, size=num_remove, replace=False)))
        if idx not in seen:
            seen.add(idx)
            yield list(idx)


def structured_puncture_sets(n, num_remove):
    patterns = []
    # Contiguous blocks
    for frac in [0, 0.25, 0.5, 0.75, 1.0]:
        s = int(frac * max(n - num_remove, 0))
        s = max(0, min(s, n - num_remove))
        patterns.append(list(range(s, s + num_remove)))
    # Strided
    for stride in range(2, min(n // max(num_remove, 1) + 1, 8)):
        for off in range(min(stride, 2)):
            pat = list(range(off, n, stride))[:num_remove]
            if len(pat) == num_remove:
                patterns.append(pat)
    # Even/odd
    for s in (0, 1):
        pat = list(range(s, n, 2))[:num_remove]
        if len(pat) == num_remove:
            patterns.append(pat)
    # Deduplicate
    seen = set()
    out = []
    for p in patterns:
        k = tuple(p)
        if k not in seen and len(p) == num_remove:
            seen.add(k)
            out.append(p)
    return out


# ---------------------------------------------------------------------------
# Search
# ---------------------------------------------------------------------------

def search_one_code(name, Hx, Hz, exp_n, exp_k, exp_d,
                    num_punctures_range, num_trials, rng,
                    do_shorten=True):
    n_orig = exp_n
    results = []
    best_by_nk = {}
    stats = {"tried": 0, "valid": 0, "computed": 0, "interesting": 0}

    for operation in ["puncture", "shorten"]:
        if operation == "shorten" and not do_shorten:
            continue

        print(f"\n  === {operation.upper()} ===", flush=True)

        for num_remove in num_punctures_range:
            target_n = n_orig - num_remove
            if target_n < 10:
                break

            struct_pats = structured_puncture_sets(n_orig, num_remove)
            rand_pats = list(random_puncture_sets(n_orig, num_remove,
                                                   num_trials, rng))
            all_pats = struct_pats + rand_pats

            batch_t0 = time.time()
            batch_new = 0
            batch_skip_highk = 0

            for pi, positions in enumerate(all_pats):
                stats["tried"] += 1
                ptype = "struct" if pi < len(struct_pats) else "random"

                if operation == "puncture":
                    res = puncture(Hx, Hz, positions)
                else:
                    res = shorten(Hx, Hz, positions)
                if res is None:
                    continue

                Hx_r, Hz_r = res
                n_r = Hx_r.shape[1]
                k_r = n_r - gf2_rank(Hx_r) - gf2_rank(Hz_r)
                if k_r < 1:
                    continue

                stats["valid"] += 1
                nk = (n_r, k_r)

                # Filter: skip very high k/n codes not in bounds table
                if k_r > n_r * 0.5 and nk not in KNOWN_BOUNDS:
                    batch_skip_highk += 1
                    continue

                # Skip if we already have a good result for this (n,k)
                if nk in best_by_nk:
                    if nk not in KNOWN_BOUNDS:
                        continue
                    if best_by_nk[nk] >= KNOWN_BOUNDS[nk][0]:
                        continue

                d, method = compute_distance(Hx_r, Hz_r, timeout_sec=60)
                if d is None:
                    continue

                stats["computed"] += 1

                if nk not in best_by_nk or d > best_by_nk[nk]:
                    best_by_nk[nk] = d
                    comparison = check_against_bounds(n_r, k_r, d)

                    entry = {
                        "source": name,
                        "operation": operation,
                        "positions": positions,
                        "pattern_type": ptype,
                        "n": int(n_r), "k": int(k_r), "d": int(d),
                        "d_method": method,
                        "comparison": comparison,
                        "Hx_shape": list(Hx_r.shape),
                        "Hz_shape": list(Hz_r.shape),
                    }
                    results.append(entry)
                    batch_new += 1

                    status = comparison.get("status", "?")
                    flag = ""
                    if status == "BEATS_KNOWN":
                        flag = " ***"
                        stats["interesting"] += 1
                        entry["Hx"] = Hx_r.tolist()
                        entry["Hz"] = Hz_r.tolist()
                    elif status == "matches_best":
                        flag = " (=best)"

                    print(f"    [{operation[:5]:5s} {num_remove:2d} {ptype:6s}] "
                          f"[[{n_r},{k_r},{d}]] ({method}) "
                          f"{status}{flag}", flush=True)

            elapsed = time.time() - batch_t0
            print(f"    {operation} {num_remove:2d}: {batch_new} new, "
                  f"{batch_skip_highk} skipped(high-k), "
                  f"{elapsed:.1f}s", flush=True)

    print(f"\n  {name}: {stats['tried']} tried, {stats['valid']} valid, "
          f"{stats['computed']} distances, {stats['interesting']} beating bounds",
          flush=True)
    if best_by_nk:
        print(f"  Best codes:", flush=True)
        for (nr, kr), d in sorted(best_by_nk.items()):
            c = check_against_bounds(nr, kr, d)
            s = c.get("status", "?")
            kn = ""
            if "known_lower" in c:
                kn = f" (known: {c['known_lower']}-{c['known_upper']})"
            print(f"    [[{nr},{kr},{d}]] {s}{kn}", flush=True)

    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Puncture/shorten published codes to find good smaller codes"
    )
    parser.add_argument("--code", type=str, default=None)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--num-punctures", type=int, default=10)
    parser.add_argument("--min-punctures", type=int, default=5)
    parser.add_argument("--trials", type=int, default=200)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--no-shorten", action="store_true")
    args = parser.parse_args()

    if not args.code and not args.all:
        parser.error("Specify --code NAME or --all")

    rng = np.random.default_rng(args.seed)

    print("=" * 70, flush=True)
    print("Puncture & Shorten Search", flush=True)
    print(f"  {datetime.now(timezone.utc).isoformat()}", flush=True)
    print(f"  trials={args.trials}, range={args.min_punctures}-{args.num_punctures}",
          flush=True)
    print("=" * 70, flush=True)

    print("\nBuilding codes...", flush=True)
    all_codes = build_published_codes()

    if args.code:
        cn = args.code.lower().replace("-", "").replace("_", "")
        if cn not in all_codes:
            print(f"Unknown: {args.code}. Available: {', '.join(all_codes.keys())}",
                  flush=True)
            sys.exit(1)
        to_search = {cn: all_codes[cn]}
    else:
        to_search = all_codes

    punct_range = range(args.min_punctures, args.num_punctures + 1)
    all_results = []

    for name, (Hx, Hz, en, ek, ed) in to_search.items():
        print(f"\n{'='*70}", flush=True)
        print(f"{name} [[{en},{ek},{ed}]]", flush=True)
        print(f"{'='*70}", flush=True)
        results = search_one_code(name, Hx, Hz, en, ek, ed,
                                   punct_range, args.trials, rng,
                                   do_shorten=not args.no_shorten)
        all_results.extend(results)

    # Save
    out_dir = os.path.join(os.path.dirname(__file__), "..", "results")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "puncture_results.json")

    beating = [r for r in all_results
               if r.get("comparison", {}).get("status") == "BEATS_KNOWN"]
    matching = [r for r in all_results
                if r.get("comparison", {}).get("status") == "matches_best"]

    output = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "config": {"codes": list(to_search.keys()),
                    "puncture_range": [args.min_punctures, args.num_punctures],
                    "trials": args.trials, "seed": args.seed},
        "results": all_results,
        "summary": {"total": len(all_results),
                     "beating_known": len(beating),
                     "matching_best": len(matching)},
    }

    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.integer,)):
                return int(obj)
            if isinstance(obj, (np.floating,)):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return super().default(obj)

    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, cls=NumpyEncoder)
    print(f"\nSaved to {out_path}", flush=True)

    print(f"\n{'='*70}", flush=True)
    print(f"SUMMARY: {len(all_results)} codes, "
          f"{len(beating)} beating bounds, {len(matching)} matching best",
          flush=True)
    for r in beating + matching:
        c = r["comparison"]
        kn = f" (known: d={c.get('known_lower','?')})" if "known_lower" in c else ""
        print(f"  [[{r['n']},{r['k']},{r['d']}]] {r['source']} "
              f"{r['operation']} -- {c['status']}{kn}", flush=True)
    print("=" * 70, flush=True)


if __name__ == "__main__":
    main()
