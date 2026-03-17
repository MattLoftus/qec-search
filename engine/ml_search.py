"""
ML-guided search for quantum error correcting codes.

Trains a RandomForest model on (polynomial features -> distance) data from
QCCode constructions, then uses the model to prioritize which polynomial
pairs to evaluate with exact (expensive) distance computation.

Usage:
  cd ~/workspace/qec-search/engine
  ../venv/bin/python3 ml_search.py --all
  ../venv/bin/python3 ml_search.py --collect   # generate training data only
  ../venv/bin/python3 ml_search.py --train     # train model only
  ../venv/bin/python3 ml_search.py --search    # ML-guided search only
"""

import sys
import os
import json
import time
import signal
import random
import argparse
import functools
import math
import pickle
from collections import defaultdict

print = functools.partial(print, flush=True)

# Be gentle on CPU
os.nice(15)

import numpy as np
import sympy
from qldpc.codes import QCCode

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from known_codes import KNOWN_BOUNDS

try:
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import cross_val_score
except ImportError:
    print("Installing scikit-learn...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "scikit-learn"])
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import cross_val_score

x = sympy.Symbol('x')

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
BEST_CODES_PATH = os.path.join(RESULTS_DIR, "best_codes.json")
TRAINING_DATA_PATH = os.path.join(RESULTS_DIR, "ml_training_data.json")
ML_RESULTS_PATH = os.path.join(RESULTS_DIR, "ml_search_results.json")
MODEL_PATH = os.path.join(RESULTS_DIR, "ml_model.pkl")


# ---------------------------------------------------------------------------
# Timeout-protected distance computation
# ---------------------------------------------------------------------------

class DistanceTimeout(Exception):
    pass


def _alarm_handler(signum, frame):
    raise DistanceTimeout()


def get_distance_safe(code, timeout_sec=15):
    """Compute code distance with SIGALRM timeout. Returns (d, method) or (None, 'error')."""
    old_handler = signal.signal(signal.SIGALRM, _alarm_handler)
    signal.alarm(timeout_sec)
    try:
        d = code.get_distance()
        signal.alarm(0)
        return d, "exact"
    except DistanceTimeout:
        signal.alarm(0)
        return None, "timeout"
    except Exception:
        signal.alarm(0)
        return None, "error"
    finally:
        signal.signal(signal.SIGALRM, old_handler)


# ---------------------------------------------------------------------------
# Best codes persistence (shared with other modules)
# ---------------------------------------------------------------------------

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
            'Hx': [], 'Hz': [],
            'hash': f'ml_search_{n}_{k}_{d}',
            'max_weight': 0, 'avg_weight': 0,
            'rate': k / n,
        }
        return True, list(sorted(by_nk.values(), key=lambda x: (x['n'], x['k'])))
    return False, existing


def check_bounds(n, k, d):
    """Check against known bounds."""
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
# Feature engineering
# ---------------------------------------------------------------------------

def _gcd(a, b):
    """GCD for non-negative integers."""
    while b:
        a, b = b, a % b
    return a


def _is_prime(n):
    """Simple primality check."""
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


def _is_prime_power(n):
    """Check if n is a prime power (p^k for some prime p, k >= 1)."""
    if n < 2:
        return False
    for p in range(2, int(math.sqrt(n)) + 1):
        if n % p == 0:
            while n % p == 0:
                n //= p
            return n == 1
    return True  # n itself is prime


def _num_divisors(n):
    """Count number of divisors."""
    if n < 1:
        return 0
    count = 0
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            count += 2 if i * i != n else 1
    return count


def extract_features(order, exp_a, exp_b):
    """Extract feature vector from (order, exponents_a, exponents_b).

    Supports weight-3 (2 exponents each) and weight-4 (3 exponents each).
    Always returns a fixed-size feature vector by padding shorter exponent
    lists with zeros.
    """
    # Normalize exponents to sorted lists, pad to length 3
    a = sorted(exp_a)
    b = sorted(exp_b)
    while len(a) < 3:
        a.append(0)
    while len(b) < 3:
        b.append(0)

    a1, a2, a3 = a[0], a[1], a[2]
    b1, b2, b3 = b[0], b[1], b[2]
    n = 2 * order

    features = []

    # Basic parameters
    features.append(order)                          # 0: order
    features.append(order / 30.0)                   # 1: order normalized
    features.append(n)                              # 2: n = 2*order

    # Exponent values (normalized by order)
    features.append(a1 / order)                     # 3
    features.append(a2 / order)                     # 4
    features.append(a3 / order)                     # 5: 0 for weight-3
    features.append(b1 / order)                     # 6
    features.append(b2 / order)                     # 7
    features.append(b3 / order)                     # 8: 0 for weight-3

    # Raw exponent values
    features.append(a1)                             # 9
    features.append(a2)                             # 10
    features.append(a3)                             # 11
    features.append(b1)                             # 12
    features.append(b2)                             # 13
    features.append(b3)                             # 14

    # Exponent gaps (within each polynomial)
    features.append(abs(a2 - a1))                   # 15
    features.append(abs(a3 - a2) if a3 > 0 else 0) # 16
    features.append(abs(b2 - b1))                   # 17
    features.append(abs(b3 - b2) if b3 > 0 else 0) # 18

    # Cross-polynomial gaps
    features.append(abs(a1 - b1))                   # 19
    features.append(abs(a2 - b2))                   # 20
    features.append(abs(a1 - b2))                   # 21
    features.append(abs(a2 - b1))                   # 22

    # GCDs within polynomials
    features.append(_gcd(a1, a2) if a2 > 0 else a1)  # 23
    features.append(_gcd(b1, b2) if b2 > 0 else b1)  # 24

    # GCDs with order
    features.append(_gcd(order, a1))                # 25
    features.append(_gcd(order, a2))                # 26
    features.append(_gcd(order, b1))                # 27
    features.append(_gcd(order, b2))                # 28

    # GCD of all exponents with order
    all_exp = [e for e in [a1, a2, a3, b1, b2, b3] if e > 0]
    g = all_exp[0] if all_exp else 0
    for e in all_exp[1:]:
        g = _gcd(g, e)
    features.append(_gcd(order, g))                 # 29

    # Number-theoretic properties of order
    features.append(1.0 if _is_prime(order) else 0.0)         # 30
    features.append(1.0 if _is_prime_power(order) else 0.0)   # 31
    features.append(_num_divisors(order))                       # 32
    features.append(order % 2)                                  # 33: parity

    # Polynomial symmetry
    a_sorted = tuple(sorted(e for e in [a1, a2, a3] if e > 0))
    b_sorted = tuple(sorted(e for e in [b1, b2, b3] if e > 0))
    features.append(1.0 if a_sorted == b_sorted else 0.0)     # 34: identical
    # Similarity: fraction of shared exponents
    a_set = set(e for e in [a1, a2, a3] if e > 0)
    b_set = set(e for e in [b1, b2, b3] if e > 0)
    if a_set or b_set:
        features.append(len(a_set & b_set) / len(a_set | b_set))  # 35: Jaccard
    else:
        features.append(0.0)

    # Weight info
    weight_a = len([e for e in [a1, a2, a3] if e > 0]) + 1  # +1 for constant term
    weight_b = len([e for e in [b1, b2, b3] if e > 0]) + 1
    features.append(weight_a)                                  # 36
    features.append(weight_b)                                  # 37

    # Sum and product of exponents (mod order)
    sum_a = sum(e for e in [a1, a2, a3] if e > 0) % order
    sum_b = sum(e for e in [b1, b2, b3] if e > 0) % order
    features.append(sum_a / order)                             # 38
    features.append(sum_b / order)                             # 39

    # Min/max exponent spread
    all_nonzero = [e for e in [a1, a2, a3, b1, b2, b3] if e > 0]
    if all_nonzero:
        features.append(min(all_nonzero) / order)             # 40
        features.append(max(all_nonzero) / order)             # 41
        features.append((max(all_nonzero) - min(all_nonzero)) / order)  # 42
    else:
        features.extend([0.0, 0.0, 0.0])

    return features


FEATURE_DIM = len(extract_features(10, [1, 3], [2, 5]))


# ---------------------------------------------------------------------------
# Phase 1: Collect training data
# ---------------------------------------------------------------------------

def collect_training_data(orders=range(9, 26), samples_per_order=500, timeout_sec=15, seed=42):
    """Generate training data by evaluating random QCCode polynomial pairs.

    For each order, generate random weight-3 polynomial pairs, compute exact
    distance, and store (features, distance, n, k, order, exp_a, exp_b).

    Many orders (especially primes) yield k=0 for most polynomial pairs, so
    we try up to 10x samples_per_order attempts and also try all exhaustive
    pairs for small orders.
    """
    rng = random.Random(seed)
    training_data = []

    print(f"Collecting training data for orders {orders[0]}-{orders[-1]}, "
          f"target ~{samples_per_order} valid samples per order")
    print(f"Timeout: {timeout_sec}s per distance computation")
    print("-" * 60)

    total_tried = 0
    total_valid = 0
    total_timeout = 0
    t_start = time.time()

    for order in orders:
        order_start = time.time()
        valid = 0
        tried = 0
        timeouts = 0
        seen = set()

        # For small orders, exhaustive search is feasible
        max_attempts = max(samples_per_order * 15, (order - 1) ** 4)
        # Cap at a reasonable number
        max_attempts = min(max_attempts, samples_per_order * 20)

        while tried < max_attempts:
            tried += 1
            total_tried += 1

            # Generate random weight-3 polynomial: 1 + x^a1 + x^a2
            exps = sorted(rng.sample(range(1, order), 2))
            a1, a2 = exps
            exps = sorted(rng.sample(range(1, order), 2))
            b1, b2 = exps

            key = (order, a1, a2, b1, b2)
            if key in seen:
                continue
            seen.add(key)

            try:
                poly_a = 1 + x**a1 + x**a2
                poly_b = 1 + x**b1 + x**b2
                code = QCCode([order], poly_a, poly_b)
                n = code.num_qubits
                k = code.dimension

                if k < 1:
                    continue

                d, method = get_distance_safe(code, timeout_sec=timeout_sec)

                if d is None:
                    if method == "timeout":
                        timeouts += 1
                        total_timeout += 1
                    continue
                if d < 1:
                    continue

                features = extract_features(order, [a1, a2], [b1, b2])
                training_data.append({
                    'features': features,
                    'distance': d,
                    'n': n,
                    'k': k,
                    'order': order,
                    'exp_a': [a1, a2],
                    'exp_b': [b1, b2],
                    'method': method,
                })
                valid += 1
                total_valid += 1

                # Stop early if we have enough for this order
                if valid >= samples_per_order:
                    break

            except Exception:
                continue

        order_elapsed = time.time() - order_start
        total_elapsed = time.time() - t_start
        print(f"  order={order:3d} (n={2*order:3d}) | "
              f"{valid:4d} valid / {tried:4d} tried ({100*valid/max(tried,1):.0f}%) | "
              f"{timeouts:3d} timeouts | "
              f"{order_elapsed:.1f}s | "
              f"total: {total_valid} samples in {total_elapsed:.0f}s")

    # Save training data
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(TRAINING_DATA_PATH, 'w') as f:
        json.dump(training_data, f)

    elapsed = time.time() - t_start
    print(f"\nCollection complete: {total_valid} valid samples from "
          f"{total_tried} attempts in {elapsed:.0f}s")
    print(f"Timeouts: {total_timeout}")
    print(f"Saved to {TRAINING_DATA_PATH}")

    # Print distance distribution
    distances = [d['distance'] for d in training_data]
    if distances:
        dist_counts = defaultdict(int)
        for d in distances:
            dist_counts[d] += 1
        print(f"\nDistance distribution:")
        for d in sorted(dist_counts.keys()):
            bar = "#" * min(dist_counts[d] // 5, 50)
            print(f"  d={d:3d}: {dist_counts[d]:5d} {bar}")

    return training_data


# ---------------------------------------------------------------------------
# Phase 2: Train model
# ---------------------------------------------------------------------------

def train_model(training_data=None):
    """Train a RandomForestRegressor to predict distance from polynomial features."""
    # Load training data if not provided
    if training_data is None:
        if not os.path.exists(TRAINING_DATA_PATH):
            print(f"ERROR: No training data at {TRAINING_DATA_PATH}")
            print("Run with --collect first.")
            return None
        with open(TRAINING_DATA_PATH) as f:
            training_data = json.load(f)

    if len(training_data) < 20:
        print(f"ERROR: Only {len(training_data)} training samples (need >= 20)")
        return None

    print(f"Training on {len(training_data)} samples, {FEATURE_DIM} features")

    # Build feature matrix and target vector
    X = np.array([d['features'] for d in training_data], dtype=np.float64)
    y = np.array([d['distance'] for d in training_data], dtype=np.float64)

    print(f"Target stats: min={y.min():.0f}, max={y.max():.0f}, "
          f"mean={y.mean():.2f}, std={y.std():.2f}")

    # Train RandomForest
    model = RandomForestRegressor(
        n_estimators=200,
        max_depth=20,
        min_samples_leaf=3,
        max_features='sqrt',
        n_jobs=-1,
        random_state=42,
    )

    # Cross-validation
    print("Running 5-fold cross-validation...")
    cv_scores = cross_val_score(model, X, y, cv=5, scoring='r2')
    print(f"  CV R^2: {cv_scores.mean():.4f} (+/- {cv_scores.std():.4f})")
    cv_mae = cross_val_score(model, X, y, cv=5, scoring='neg_mean_absolute_error')
    print(f"  CV MAE: {-cv_mae.mean():.4f} (+/- {cv_mae.std():.4f})")

    # Train on full dataset
    print("Training final model on all data...")
    model.fit(X, y)

    # Feature importances
    importances = model.feature_importances_
    feature_names = [
        'order', 'order_norm', 'n',
        'a1_norm', 'a2_norm', 'a3_norm', 'b1_norm', 'b2_norm', 'b3_norm',
        'a1', 'a2', 'a3', 'b1', 'b2', 'b3',
        'gap_a12', 'gap_a23', 'gap_b12', 'gap_b23',
        'cross_a1b1', 'cross_a2b2', 'cross_a1b2', 'cross_a2b1',
        'gcd_a12', 'gcd_b12',
        'gcd_ord_a1', 'gcd_ord_a2', 'gcd_ord_b1', 'gcd_ord_b2',
        'gcd_all',
        'is_prime', 'is_prime_power', 'num_divisors', 'parity',
        'poly_identical', 'poly_jaccard',
        'weight_a', 'weight_b',
        'sum_a_norm', 'sum_b_norm',
        'min_exp_norm', 'max_exp_norm', 'exp_spread',
    ]
    top_indices = np.argsort(importances)[::-1][:15]
    print("\nTop 15 feature importances:")
    for i, idx in enumerate(top_indices):
        name = feature_names[idx] if idx < len(feature_names) else f"feat_{idx}"
        print(f"  {i+1:2d}. {name:20s} {importances[idx]:.4f}")

    # Save model
    with open(MODEL_PATH, 'wb') as f:
        pickle.dump(model, f)
    print(f"\nModel saved to {MODEL_PATH}")

    return model


# ---------------------------------------------------------------------------
# Phase 3: ML-guided search
# ---------------------------------------------------------------------------

def ml_guided_search(model=None, order_range=range(9, 31), candidates_per_order=100_000,
                     top_k=200, timeout_sec=15, seed=123):
    """Use trained model to guide code search.

    For each order:
    1. Generate many random polynomial pairs
    2. Predict distance with model
    3. Evaluate only the top predictions with exact distance
    """
    # Load model if not provided
    if model is None:
        if not os.path.exists(MODEL_PATH):
            print(f"ERROR: No trained model at {MODEL_PATH}")
            print("Run with --train first.")
            return []
        with open(MODEL_PATH, 'rb') as f:
            model = pickle.load(f)

    rng = random.Random(seed)
    best_codes = load_best_codes()
    all_results = []
    improvements = 0
    matches = 0

    print(f"ML-guided search: orders {order_range[0]}-{order_range[-1]}")
    print(f"  {candidates_per_order:,} candidates per order, top {top_k} evaluated exactly")
    print(f"  Timeout: {timeout_sec}s per exact distance computation")
    print("=" * 70)

    t_start = time.time()

    for order in order_range:
        order_start = time.time()

        # Generate random polynomial pairs and extract features
        candidates = []
        seen = set()
        max_exps = order - 1

        if max_exps < 2:
            continue

        for _ in range(candidates_per_order):
            a1, a2 = sorted(rng.sample(range(1, order), 2))
            b1, b2 = sorted(rng.sample(range(1, order), 2))
            key = (a1, a2, b1, b2)
            if key in seen:
                continue
            seen.add(key)
            candidates.append((a1, a2, b1, b2))

        if not candidates:
            continue

        # Extract features for all candidates (vectorized-ish)
        features = np.array(
            [extract_features(order, [a1, a2], [b1, b2])
             for a1, a2, b1, b2 in candidates],
            dtype=np.float64
        )

        # Predict distances
        predictions = model.predict(features)

        # Take top-k by predicted distance
        top_indices = np.argsort(predictions)[::-1][:top_k]

        # Evaluate exact distance on top candidates
        order_results = []
        order_improvements = 0
        order_matches = 0
        exact_count = 0
        best_d_found = 0

        for idx in top_indices:
            a1, a2, b1, b2 = candidates[idx]
            pred_d = predictions[idx]

            try:
                poly_a = 1 + x**a1 + x**a2
                poly_b = 1 + x**b1 + x**b2
                code = QCCode([order], poly_a, poly_b)
                n = code.num_qubits
                k = code.dimension

                if k < 1:
                    continue

                d, method = get_distance_safe(code, timeout_sec=timeout_sec)
                if d is None or d < 1:
                    continue

                exact_count += 1
                if d > best_d_found:
                    best_d_found = d

                status, dl, du = check_bounds(n, k, d)
                updated, best_codes = update_best(
                    best_codes, n, k, d,
                    f"ml_qc({order}, 1+x^{a1}+x^{a2}, 1+x^{b1}+x^{b2})"
                )

                result = {
                    'order': order,
                    'n': n,
                    'k': k,
                    'd': d,
                    'exp_a': [a1, a2],
                    'exp_b': [b1, b2],
                    'predicted_d': float(pred_d),
                    'method': method,
                    'status': status,
                }

                if status == "IMPROVEMENT":
                    print(f"  ******* IMPROVEMENT: [[{n},{k},{d}]] "
                          f"QC({order}, 1+x^{a1}+x^{a2}, 1+x^{b1}+x^{b2}) "
                          f"(known d={dl}, upper={du}, pred={pred_d:.1f}) *******")
                    order_improvements += 1
                    improvements += 1
                elif status == "MATCHES":
                    if d >= 4:  # Only log interesting matches
                        print(f"  MATCHES: [[{n},{k},{d}]] "
                              f"QC({order}, 1+x^{a1}+x^{a2}, 1+x^{b1}+x^{b2}) "
                              f"(pred={pred_d:.1f})")
                    order_matches += 1
                    matches += 1

                order_results.append(result)
                all_results.append(result)

                if updated:
                    save_best_codes(best_codes)

            except Exception:
                continue

        order_elapsed = time.time() - order_start
        total_elapsed = time.time() - t_start
        print(f"order={order:3d} (n={2*order:3d}) | "
              f"{len(seen):6d} unique candidates | "
              f"{exact_count:3d} evaluated | "
              f"best_d={best_d_found:2d} | "
              f"{order_matches:2d} match, {order_improvements:2d} improve | "
              f"{order_elapsed:.1f}s | total {total_elapsed:.0f}s")

    # Save results
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(ML_RESULTS_PATH, 'w') as f:
        json.dump(all_results, f, indent=2)

    elapsed = time.time() - t_start
    print(f"\n{'=' * 70}")
    print(f"ML-guided search complete in {elapsed:.0f}s")
    print(f"Total results: {len(all_results)}")
    print(f"Matches: {matches}, Improvements: {improvements}")
    print(f"Results saved to {ML_RESULTS_PATH}")

    # Summary of best codes found per (n,k)
    best_by_nk = {}
    for r in all_results:
        key = (r['n'], r['k'])
        if key not in best_by_nk or r['d'] > best_by_nk[key]['d']:
            best_by_nk[key] = r

    if best_by_nk:
        print(f"\nBest codes found:")
        for key in sorted(best_by_nk.keys()):
            r = best_by_nk[key]
            status_tag = ""
            if r['status'] == "IMPROVEMENT":
                status_tag = " *** IMPROVEMENT ***"
            elif r['status'] == "MATCHES":
                status_tag = " (matches known)"
            print(f"  [[{r['n']},{r['k']},{r['d']}]] "
                  f"QC({r['order']}, [1,{r['exp_a']}], [1,{r['exp_b']}]) "
                  f"pred={r['predicted_d']:.1f}{status_tag}")

    return all_results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="ML-guided QEC code search"
    )
    parser.add_argument("--collect", action="store_true",
                        help="Collect training data from random QCCode evaluations")
    parser.add_argument("--train", action="store_true",
                        help="Train RandomForest model on collected data")
    parser.add_argument("--search", action="store_true",
                        help="Run ML-guided search using trained model")
    parser.add_argument("--all", action="store_true",
                        help="Run all phases: collect, train, search")
    parser.add_argument("--samples", type=int, default=500,
                        help="Samples per order during collection (default: 500)")
    parser.add_argument("--candidates", type=int, default=100_000,
                        help="Candidates per order during search (default: 100000)")
    parser.add_argument("--top-k", type=int, default=200,
                        help="Top-k candidates to evaluate exactly (default: 200)")
    parser.add_argument("--timeout", type=int, default=15,
                        help="Distance computation timeout in seconds (default: 15)")
    args = parser.parse_args()

    if not (args.collect or args.train or args.search or args.all):
        parser.print_help()
        sys.exit(1)

    print("ML-Guided QEC Code Search")
    print(f"Python: {sys.version.split()[0]}")
    print(f"Feature dimension: {FEATURE_DIM}")
    print("=" * 70)

    training_data = None
    model = None

    if args.collect or args.all:
        print("\n--- Phase 1: Collecting Training Data ---")
        training_data = collect_training_data(
            orders=range(9, 26),
            samples_per_order=args.samples,
            timeout_sec=args.timeout,
        )

    if args.train or args.all:
        print("\n--- Phase 2: Training Model ---")
        model = train_model(training_data)

    if args.search or args.all:
        print("\n--- Phase 3: ML-Guided Search ---")
        ml_guided_search(
            model=model,
            order_range=range(9, 31),
            candidates_per_order=args.candidates,
            top_k=args.top_k,
            timeout_sec=args.timeout,
        )


if __name__ == "__main__":
    main()
