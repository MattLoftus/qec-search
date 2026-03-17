"""
Genetic algorithm for evolving quasi-cyclic (QC) code polynomials.

Evolves polynomial exponents directly for QCCode construction via qLDPC.
Genome: (order, [a1,a2,...], [b1,b2,...]) where polynomials are
  1 + x^a1 + x^a2 + ...  and  1 + x^b1 + x^b2 + ...

Usage:
  cd ~/workspace/qec-search/engine
  ../venv/bin/python3 evolve_polys.py --pop 50 --gens 200 --weight 3
"""

import sys
import os
import json
import time
import signal
import random
import argparse
import functools
import copy

print = functools.partial(print, flush=True)

# Be gentle on CPU
os.nice(15)

import sympy
from qldpc.codes import QCCode

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from known_codes import KNOWN_BOUNDS

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
BEST_CODES_PATH = os.path.join(RESULTS_DIR, "best_codes.json")

x = sympy.Symbol('x')


# ---------------------------------------------------------------------------
# Timeout-protected distance computation
# ---------------------------------------------------------------------------

class DistanceTimeout(Exception):
    pass


def _alarm_handler(signum, frame):
    raise DistanceTimeout()


def get_distance_safe(code, timeout_sec=30):
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
# Genome representation
# ---------------------------------------------------------------------------

class Individual:
    """A QC code candidate: (order, exponents_a, exponents_b)."""

    __slots__ = ['order', 'exp_a', 'exp_b', 'fitness', 'n', 'k', 'd', 'd_method']

    def __init__(self, order, exp_a, exp_b):
        self.order = order
        self.exp_a = sorted(set(exp_a))
        self.exp_b = sorted(set(exp_b))
        self.fitness = -1.0
        self.n = 0
        self.k = 0
        self.d = 0
        self.d_method = ""

    def poly_a(self):
        p = 1
        for e in self.exp_a:
            p += x ** e
        return p

    def poly_b(self):
        p = 1
        for e in self.exp_b:
            p += x ** e
        return p

    def desc(self):
        a_str = "+".join(["1"] + [f"x^{e}" for e in self.exp_a])
        b_str = "+".join(["1"] + [f"x^{e}" for e in self.exp_b])
        return f"QC({self.order}, {a_str}, {b_str})"

    def key(self):
        return (self.order, tuple(self.exp_a), tuple(self.exp_b))

    def copy(self):
        ind = Individual(self.order, list(self.exp_a), list(self.exp_b))
        ind.fitness = self.fitness
        ind.n = self.n
        ind.k = self.k
        ind.d = self.d
        ind.d_method = self.d_method
        return ind


# ---------------------------------------------------------------------------
# Fitness evaluation
# ---------------------------------------------------------------------------

def evaluate(ind, timeout_sec=30):
    """Evaluate an individual. Sets fitness and code params. Returns True if valid."""
    # Ensure exponents are valid
    ind.exp_a = sorted(set(e % ind.order for e in ind.exp_a if e % ind.order != 0))
    ind.exp_b = sorted(set(e % ind.order for e in ind.exp_b if e % ind.order != 0))

    if len(ind.exp_a) == 0 or len(ind.exp_b) == 0:
        ind.fitness = -1
        return False

    # Check for duplicate exponents after mod
    if len(ind.exp_a) != len(set(ind.exp_a)) or len(ind.exp_b) != len(set(ind.exp_b)):
        ind.fitness = -1
        return False

    try:
        code = QCCode([ind.order], ind.poly_a(), ind.poly_b())
        ind.n = code.num_qubits
        ind.k = code.dimension
    except Exception:
        ind.fitness = -1
        return False

    if ind.k < 1:
        ind.fitness = -1
        return False

    d, method = get_distance_safe(code, timeout_sec=timeout_sec)
    if d is None or d < 1:
        ind.fitness = -1
        ind.d_method = method
        return False

    ind.d = d
    ind.d_method = method

    # Fitness: primary = distance, secondary = prefer compact codes (smaller n for same d)
    # Also bonus for higher k
    ind.fitness = d * 100 - ind.n * 0.1 + ind.k * 0.5

    return True


# ---------------------------------------------------------------------------
# Genetic operators
# ---------------------------------------------------------------------------

def random_individual(order_range, weight, rng):
    """Create a random individual."""
    order = rng.randint(*order_range)
    num_exp = weight - 1  # weight includes the constant term 1
    exp_a = sorted(rng.sample(range(1, order), min(num_exp, order - 1)))
    exp_b = sorted(rng.sample(range(1, order), min(num_exp, order - 1)))
    return Individual(order, exp_a, exp_b)


def mutate(ind, order_range, weight, rng):
    """Mutate an individual. Returns a new Individual."""
    child = ind.copy()
    num_exp = weight - 1

    # Pick a mutation type
    mut_type = rng.random()

    if mut_type < 0.15:
        # Mutate order by +-1 or +-2
        delta = rng.choice([-2, -1, 1, 2])
        child.order = max(order_range[0], min(order_range[1] - 1, child.order + delta))
        # Clamp exponents to new order
        child.exp_a = [e % child.order for e in child.exp_a if e % child.order != 0]
        child.exp_b = [e % child.order for e in child.exp_b if e % child.order != 0]
        # Re-fill if needed
        while len(child.exp_a) < num_exp and child.order > num_exp:
            new_e = rng.randint(1, child.order - 1)
            if new_e not in child.exp_a:
                child.exp_a.append(new_e)
        while len(child.exp_b) < num_exp and child.order > num_exp:
            new_e = rng.randint(1, child.order - 1)
            if new_e not in child.exp_b:
                child.exp_b.append(new_e)

    elif mut_type < 0.55:
        # Mutate one exponent in poly_a
        if child.exp_a:
            idx = rng.randint(0, len(child.exp_a) - 1)
            op = rng.random()
            if op < 0.3:
                child.exp_a[idx] = max(1, min(child.order - 1, child.exp_a[idx] + rng.choice([-1, 1])))
            elif op < 0.6:
                child.exp_a[idx] = max(1, min(child.order - 1, child.exp_a[idx] + rng.choice([-2, 2])))
            else:
                child.exp_a[idx] = rng.randint(1, child.order - 1)

    elif mut_type < 0.95:
        # Mutate one exponent in poly_b
        if child.exp_b:
            idx = rng.randint(0, len(child.exp_b) - 1)
            op = rng.random()
            if op < 0.3:
                child.exp_b[idx] = max(1, min(child.order - 1, child.exp_b[idx] + rng.choice([-1, 1])))
            elif op < 0.6:
                child.exp_b[idx] = max(1, min(child.order - 1, child.exp_b[idx] + rng.choice([-2, 2])))
            else:
                child.exp_b[idx] = rng.randint(1, child.order - 1)

    else:
        # Swap poly_a and poly_b
        child.exp_a, child.exp_b = child.exp_b, child.exp_a

    # Deduplicate exponents
    child.exp_a = sorted(set(child.exp_a))
    child.exp_b = sorted(set(child.exp_b))

    return child


def crossover(parent1, parent2, rng):
    """Crossover two parents. Returns a new Individual."""
    cross_type = rng.random()

    if cross_type < 0.3:
        # Swap entire polynomial a
        child = Individual(parent1.order, list(parent2.exp_a), list(parent1.exp_b))
    elif cross_type < 0.6:
        # Swap entire polynomial b
        child = Individual(parent1.order, list(parent1.exp_a), list(parent2.exp_b))
    elif cross_type < 0.8:
        # Use parent2's order with parent1's exponents (clamped)
        order = parent2.order
        exp_a = [e % order for e in parent1.exp_a if e % order != 0]
        exp_b = [e % order for e in parent1.exp_b if e % order != 0]
        child = Individual(order, exp_a, exp_b)
    else:
        # Mix individual exponents
        exp_a = []
        for i in range(max(len(parent1.exp_a), len(parent2.exp_a))):
            if i < len(parent1.exp_a) and i < len(parent2.exp_a):
                exp_a.append(parent1.exp_a[i] if rng.random() < 0.5 else parent2.exp_a[i])
            elif i < len(parent1.exp_a):
                exp_a.append(parent1.exp_a[i])
            else:
                exp_a.append(parent2.exp_a[i])
        exp_b = []
        for i in range(max(len(parent1.exp_b), len(parent2.exp_b))):
            if i < len(parent1.exp_b) and i < len(parent2.exp_b):
                exp_b.append(parent1.exp_b[i] if rng.random() < 0.5 else parent2.exp_b[i])
            elif i < len(parent1.exp_b):
                exp_b.append(parent1.exp_b[i])
            else:
                exp_b.append(parent2.exp_b[i])
        child = Individual(parent1.order, exp_a, exp_b)

    child.exp_a = sorted(set(child.exp_a))
    child.exp_b = sorted(set(child.exp_b))
    return child


def tournament_select(population, k, rng):
    """Tournament selection: pick k random, return best."""
    contestants = rng.sample(population, min(k, len(population)))
    return max(contestants, key=lambda ind: ind.fitness)


# ---------------------------------------------------------------------------
# Seed population from known good codes
# ---------------------------------------------------------------------------

SEED_CODES = [
    # QC(24, 1+x^5+x^7, 1+x^7+x^8) -> [[48,4,8]]
    (24, [5, 7], [7, 8]),
    # QC(21, 1+x+x^3, 1+x^2+x^13) -> [[42,6,6]]
    (21, [1, 3], [2, 13]),
    # QC(15, 1+x+x^2, 1+x^2+x^7) -> [[30,4,6]]
    (15, [1, 2], [2, 7]),
]


def make_seeds(weight):
    """Create seed individuals from known codes."""
    seeds = []
    for order, exp_a, exp_b in SEED_CODES:
        if weight == 3 and len(exp_a) == 2:
            seeds.append(Individual(order, list(exp_a), list(exp_b)))
        elif weight == 4:
            # Pad weight-3 seeds with an extra exponent
            ea = list(exp_a)
            eb = list(exp_b)
            if len(ea) < 3:
                for candidate in range(1, order):
                    if candidate not in ea:
                        ea.append(candidate)
                        break
            if len(eb) < 3:
                for candidate in range(1, order):
                    if candidate not in eb:
                        eb.append(candidate)
                        break
            seeds.append(Individual(order, ea, eb))
    return seeds


# ---------------------------------------------------------------------------
# Best codes persistence
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


def update_best(existing, ind):
    """Update best codes list if this individual improves on existing. Returns (updated, list)."""
    by_nk = {}
    for e in existing:
        key = (e['n'], e['k'])
        if key not in by_nk or (e.get('d') or 0) > (by_nk[key].get('d') or 0):
            by_nk[key] = e

    key = (ind.n, ind.k)
    if key not in by_nk or ind.d > (by_nk[key].get('d') or 0):
        by_nk[key] = {
            'type': 'css',
            'n': ind.n,
            'k': ind.k,
            'd': ind.d,
            'strategy': f"evolve_qc: {ind.desc()}",
            'Hx': [],
            'Hz': [],
            'hash': f'evolve_qc_{ind.order}_{ind.n}_{ind.k}_{ind.d}',
            'max_weight': 0,
            'avg_weight': 0,
            'rate': ind.k / ind.n,
        }
        return True, list(sorted(by_nk.values(), key=lambda x: (x['n'], x['k'])))
    return False, existing


def check_bounds(n, k, d):
    """Check against known bounds. Returns (status, d_lower, d_upper)."""
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
# Main GA loop
# ---------------------------------------------------------------------------

def evolve(pop_size=50, num_gens=200, weight=3, order_range=(9, 31), seed=None):
    """Run the genetic algorithm."""
    rng = random.Random(seed)
    num_exp = weight - 1  # number of exponents (weight includes constant 1)

    print(f"Evolving QC code polynomials")
    print(f"  Population: {pop_size}, Generations: {num_gens}")
    print(f"  Polynomial weight: {weight} (1 + {num_exp} terms)")
    print(f"  Order range: {order_range[0]}-{order_range[1]-1} (n={order_range[0]*2}-{(order_range[1]-1)*2})")
    print(f"  Timeout: 30s (n<=40), 15s (n>40)")
    print("=" * 70)

    best_codes = load_best_codes()

    # Initialize population
    print("Initializing population...")
    population = []
    seen_keys = set()

    # Add seeds
    seeds = make_seeds(weight)
    for s in seeds:
        if s.key() not in seen_keys:
            timeout = 30 if s.order * 2 <= 40 else 15
            if evaluate(s, timeout_sec=timeout):
                population.append(s)
                seen_keys.add(s.key())
                status, dl, du = check_bounds(s.n, s.k, s.d)
                label = ""
                if status == "MATCHES":
                    label = f" (MATCHES known d={dl})"
                elif status == "IMPROVEMENT":
                    label = f" *** IMPROVEMENT over known d={dl} ***"
                print(f"  Seed: [[{s.n},{s.k},{s.d}]] {s.desc()}{label}")

    # Fill with random individuals
    init_attempts = 0
    while len(population) < pop_size and init_attempts < pop_size * 10:
        init_attempts += 1
        ind = random_individual(order_range, weight, rng)
        if ind.key() in seen_keys:
            continue
        timeout = 30 if ind.order * 2 <= 40 else 15
        if evaluate(ind, timeout_sec=timeout):
            population.append(ind)
            seen_keys.add(ind.key())

    print(f"Initialized {len(population)} valid individuals")

    if len(population) < 2:
        print("ERROR: Could not initialize enough valid individuals")
        return

    # Track overall best per (n,k)
    best_by_nk = {}
    for ind in population:
        key = (ind.n, ind.k)
        if key not in best_by_nk or ind.d > best_by_nk[key].d:
            best_by_nk[key] = ind.copy()

    # Evolution loop
    t_start = time.time()
    stagnation = 0
    prev_best_fitness = max(ind.fitness for ind in population)

    for gen in range(num_gens):
        gen_start = time.time()

        # Sort by fitness (descending)
        population.sort(key=lambda ind: ind.fitness, reverse=True)

        best = population[0]
        avg_fitness = sum(ind.fitness for ind in population) / len(population)
        avg_d = sum(ind.d for ind in population if ind.d > 0) / max(1, sum(1 for ind in population if ind.d > 0))

        # Check for stagnation
        if best.fitness > prev_best_fitness:
            stagnation = 0
            prev_best_fitness = best.fitness
        else:
            stagnation += 1

        # Adaptive mutation intensity
        mutation_boost = 1.0 + min(stagnation * 0.1, 2.0)

        # Print progress every 5 generations
        if gen % 5 == 0 or gen == num_gens - 1:
            elapsed = time.time() - t_start
            print(f"Gen {gen:4d} | best [[{best.n},{best.k},{best.d}]] "
                  f"fit={best.fitness:.1f} | avg_d={avg_d:.1f} "
                  f"| pop={len(population)} | stag={stagnation} "
                  f"| {elapsed:.0f}s")

        # Elitism: keep top 10%
        elite_count = max(2, pop_size // 10)
        new_population = [ind.copy() for ind in population[:elite_count]]
        new_keys = set(ind.key() for ind in new_population)

        # Fill the rest
        attempts = 0
        max_attempts = pop_size * 8
        while len(new_population) < pop_size and attempts < max_attempts:
            attempts += 1

            if rng.random() < 0.1 * mutation_boost:
                # Inject fresh random individual (exploration)
                child = random_individual(order_range, weight, rng)
            elif rng.random() < 0.7:
                # Crossover + mutation
                p1 = tournament_select(population, 3, rng)
                p2 = tournament_select(population, 3, rng)
                child = crossover(p1, p2, rng)
                if rng.random() < 0.8:
                    child = mutate(child, order_range, weight, rng)
            else:
                # Mutation only
                parent = tournament_select(population, 3, rng)
                child = mutate(parent, order_range, weight, rng)
                # Double mutate sometimes for bigger jumps
                if rng.random() < 0.3 * mutation_boost:
                    child = mutate(child, order_range, weight, rng)

            if child.key() in new_keys:
                continue

            timeout = 30 if child.order * 2 <= 40 else 15
            if evaluate(child, timeout_sec=timeout):
                new_population.append(child)
                new_keys.add(child.key())

                # Track best per (n,k) and check bounds
                key = (child.n, child.k)
                if key not in best_by_nk or child.d > best_by_nk[key].d:
                    best_by_nk[key] = child.copy()
                    status, dl, du = check_bounds(child.n, child.k, child.d)

                    if status == "MATCHES":
                        print(f"  >>> MATCHES known bound: [[{child.n},{child.k},{child.d}]] "
                              f"{child.desc()} (known d={dl})")
                    elif status == "IMPROVEMENT":
                        print(f"  ******* IMPROVEMENT: [[{child.n},{child.k},{child.d}]] "
                              f"{child.desc()} (known d={dl}, upper={du}) *******")

                    # Update persistent best codes
                    updated, best_codes = update_best(best_codes, child)
                    if updated:
                        save_best_codes(best_codes)

        population = new_population

    # Final summary
    elapsed = time.time() - t_start
    print("\n" + "=" * 70)
    print(f"Evolution complete: {num_gens} generations in {elapsed:.1f}s")
    print(f"\nBest codes found (by [n,k]):")

    match_count = 0
    improve_count = 0
    for key in sorted(best_by_nk.keys()):
        ind = best_by_nk[key]
        status, dl, du = check_bounds(ind.n, ind.k, ind.d)
        label = ""
        if status == "MATCHES":
            label = "  <-- MATCHES KNOWN"
            match_count += 1
        elif status == "IMPROVEMENT":
            label = f"  <-- *** IMPROVEMENT (was d={dl}) ***"
            improve_count += 1
        elif status == "BELOW":
            label = f"  (known d={dl})"
        print(f"  [[{ind.n},{ind.k},{ind.d}]] {ind.desc()}{label}")

    print(f"\nSummary: {match_count} matches, {improve_count} improvements")
    print(f"Best codes saved to {BEST_CODES_PATH}")

    # Final save
    for key in best_by_nk:
        ind = best_by_nk[key]
        _, best_codes = update_best(best_codes, ind)
    save_best_codes(best_codes)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Evolve QC code polynomials via genetic algorithm")
    parser.add_argument("--pop", type=int, default=50, help="Population size (default: 50)")
    parser.add_argument("--gens", type=int, default=200, help="Number of generations (default: 200)")
    parser.add_argument("--weight", type=int, default=3, choices=[3, 4],
                        help="Polynomial weight: 3 means 1+x^a+x^b, 4 means 1+x^a+x^b+x^c (default: 3)")
    parser.add_argument("--seed", type=int, default=None, help="Random seed")
    parser.add_argument("--min-order", type=int, default=9, help="Minimum cyclic group order (default: 9)")
    parser.add_argument("--max-order", type=int, default=31, help="Maximum cyclic group order (default: 31, exclusive)")
    args = parser.parse_args()

    evolve(
        pop_size=args.pop,
        num_gens=args.gens,
        weight=args.weight,
        order_range=(args.min_order, args.max_order),
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
