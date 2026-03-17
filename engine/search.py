"""
Search strategies for finding stabilizer codes with good parameters.

Strategies:
1. Random: generate random CSS codes, compute distance
2. Genetic: evolve populations of codes toward better distance
3. Algebraic: construct codes from hypergraph products of classical codes
"""

import numpy as np
import time
import json
import os
from dataclasses import dataclass, field
from typing import List, Optional, Callable
from codes import (CSSCode, random_css_code, random_css_from_classical, hypergraph_product,
                    random_self_orthogonal, random_bicycle_code, punctured_code, gf2_rank)
from distance import exact_distance_css, estimate_distance_random
from known_codes import KNOWN_BOUNDS


@dataclass
class SearchResult:
    """A single code found during search."""
    code: CSSCode
    n: int
    k: int
    d: int
    generation: int = 0
    strategy: str = "random"
    timestamp: float = field(default_factory=time.time)
    beats_known: bool = False  # True if d > known lower bound

    def to_dict(self):
        d = self.code.to_dict()
        d.update({
            "generation": self.generation,
            "strategy": self.strategy,
            "timestamp": self.timestamp,
            "beats_known": self.beats_known,
        })
        return d


@dataclass
class SearchStats:
    """Statistics from a search campaign."""
    total_codes_tested: int = 0
    valid_codes: int = 0
    codes_with_k_ge_1: int = 0
    best_distance: dict = field(default_factory=dict)  # (n,k) -> best d found
    best_codes: list = field(default_factory=list)  # List of SearchResult
    generations_completed: int = 0
    elapsed_sec: float = 0
    strategy: str = "random"


def random_search(n, target_k=1, num_samples=10000, distance_timeout=30,
                  callback=None, rng=None):
    """Random search for CSS codes with good distance.

    Args:
        n: number of physical qubits
        target_k: target number of logical qubits
        num_samples: number of random codes to try
        distance_timeout: seconds per distance computation
        callback: called with (i, code, distance) for each valid code
        rng: numpy random generator

    Returns:
        SearchStats with results
    """
    if rng is None:
        rng = np.random.default_rng()

    stats = SearchStats(strategy="random")
    start = time.time()
    seen_hashes = set()

    for i in range(num_samples):
        stats.total_codes_tested += 1

        # Rotate through construction methods
        method = i % 4
        if method == 0:
            code = random_css_code(n, target_k=target_k, density=rng.uniform(0.2, 0.5), rng=rng)
        elif method == 1:
            code = random_css_from_classical(n, rng=rng)
        elif method == 2:
            code = random_self_orthogonal(n, target_k=target_k, rng=rng)
        else:
            code = random_bicycle_code(n, rng=rng) if n % 2 == 0 else random_self_orthogonal(n, target_k=target_k, rng=rng)

        if code is None:
            continue

        valid, _ = code.is_valid()
        if not valid:
            continue
        stats.valid_codes += 1

        if code.k < 1:
            continue
        stats.codes_with_k_ge_1 += 1

        # Dedup
        h = code.code_hash()
        if h in seen_hashes:
            continue
        seen_hashes.add(h)

        # Compute distance
        d = exact_distance_css(code, timeout_sec=distance_timeout)
        if d is None or d < 1:
            continue

        code.distance = d
        key = (code.n, code.k)

        # Check against known bounds
        beats = False
        if key in KNOWN_BOUNDS:
            known_lower, known_upper = KNOWN_BOUNDS[key]
            if d > known_lower:
                beats = True

        result = SearchResult(
            code=code, n=code.n, k=code.k, d=d,
            generation=i, strategy="random", beats_known=beats,
        )

        # Track best per (n, k)
        if key not in stats.best_distance or d > stats.best_distance[key]:
            stats.best_distance[key] = d
            stats.best_codes.append(result)

        if callback:
            callback(i, code, d)

    stats.elapsed_sec = time.time() - start
    return stats


def genetic_search(n, target_k=1, pop_size=100, num_generations=200,
                   mutation_rate=0.05, tournament_size=3, elitism=5,
                   distance_timeout=30, callback=None, rng=None):
    """Genetic algorithm search for CSS codes.

    Genome: (Hx, Hz) pair of binary parity check matrices.
    Uses row-level mutations and commutativity-preserving crossover.

    Fitness: code distance (higher is better), with secondary objectives
    for stabilizer weight (lower is better).

    Args:
        n: number of physical qubits
        target_k: target logical qubits
        pop_size: population size
        num_generations: number of generations to evolve
        mutation_rate: base probability of row-level mutation
        tournament_size: tournament selection size
        elitism: number of best individuals preserved each generation
        distance_timeout: seconds per distance computation
        callback: called with (gen, best_code, best_d, pop_stats) each generation
        rng: numpy random generator

    Returns:
        SearchStats with results
    """
    if rng is None:
        rng = np.random.default_rng()

    stats = SearchStats(strategy="genetic")
    start = time.time()

    # Target number of Hx rows
    total_checks = n - target_k
    rx_target = total_checks // 2

    if rx_target < 1:
        print(f"Cannot build code with n={n}, k={target_k}: need at least 1 check row")
        return stats

    # Initialize population with exact distance computation
    population = []
    fitnesses = []
    seen_hashes = set()

    print(f"Initializing population of {pop_size} codes (n={n}, target k={target_k})...")

    for init_i in range(pop_size * 5):  # Over-generate to get enough valid ones
        # Rotate through construction methods for diversity
        method = init_i % 4
        if method == 0:
            code = random_css_code(n, target_k=target_k, density=rng.uniform(0.2, 0.5), rng=rng)
        elif method == 1:
            code = random_self_orthogonal(n, target_k=target_k, rng=rng)
        elif method == 2:
            code = random_bicycle_code(n, rng=rng) if n % 2 == 0 else random_css_code(n, target_k=target_k, rng=rng)
        else:
            code = random_css_from_classical(n, rng=rng)
        if code is None:
            continue
        if code.k < 1:
            continue

        stats.total_codes_tested += 1
        stats.valid_codes += 1
        stats.codes_with_k_ge_1 += 1

        # Dedup
        h = code.code_hash()
        if h in seen_hashes:
            continue
        seen_hashes.add(h)

        # Use exact distance for initialization (codes are small, it's fast)
        d = exact_distance_css(code, timeout_sec=distance_timeout)
        if d is None:
            d = 1
        code.distance = d

        population.append(code)
        fitnesses.append(d)

        if len(population) >= pop_size:
            break

    if len(population) < 2:
        print("Failed to initialize population")
        return stats

    print(f"Initialized {len(population)} codes (best init d={max(fitnesses)}). Starting evolution...")

    best_ever_d = 0
    best_ever_code = None
    stagnation_counter = 0
    current_mutation_rate = mutation_rate

    for gen in range(num_generations):
        stats.generations_completed = gen + 1

        # Evaluate: compute exact distance for top candidates
        # (top 20% get exact distance recomputed)
        top_k = max(2, len(population) // 5)
        top_indices = np.argsort(fitnesses)[-top_k:]

        for idx in top_indices:
            code = population[idx]
            # Recompute exact distance for promising codes
            d_exact = exact_distance_css(code, timeout_sec=distance_timeout)
            if d_exact is not None:
                code.distance = d_exact
                fitnesses[idx] = d_exact
                stats.total_codes_tested += 1
                stats.valid_codes += 1
                stats.codes_with_k_ge_1 += 1

        # Track best
        best_idx = np.argmax(fitnesses)
        gen_best_d = fitnesses[best_idx]
        gen_best_code = population[best_idx]

        # Adaptive mutation: increase when stagnating
        if gen_best_d > best_ever_d:
            best_ever_d = gen_best_d
            best_ever_code = gen_best_code
            stagnation_counter = 0
            current_mutation_rate = mutation_rate  # Reset to base

            key = (gen_best_code.n, gen_best_code.k)
            beats = False
            if key in KNOWN_BOUNDS:
                known_lower, _ = KNOWN_BOUNDS[key]
                if gen_best_d > known_lower:
                    beats = True

            result = SearchResult(
                code=gen_best_code, n=gen_best_code.n, k=gen_best_code.k,
                d=int(gen_best_d), generation=gen, strategy="genetic",
                beats_known=beats,
            )
            stats.best_codes.append(result)
            key_tuple = (gen_best_code.n, gen_best_code.k)
            stats.best_distance[key_tuple] = int(gen_best_d)
        else:
            stagnation_counter += 1
            # Ramp up mutation when stagnating
            if stagnation_counter > 5:
                current_mutation_rate = min(mutation_rate * (1 + stagnation_counter * 0.2), 0.5)

        # Population stats
        unique_hashes = len(set(c.code_hash() for c in population))
        pop_stats = {
            "generation": gen,
            "best_d": int(gen_best_d),
            "avg_d": float(np.mean(fitnesses)),
            "best_k": int(gen_best_code.k),
            "pop_size": len(population),
            "unique": unique_hashes,
            "mutation_rate": round(current_mutation_rate, 3),
            "best_weight": int(gen_best_code.max_stabilizer_weight),
        }

        if callback:
            callback(gen, gen_best_code, int(gen_best_d), pop_stats)

        if gen % 10 == 0:
            print(f"  Gen {gen:4d}: best d={int(gen_best_d)}, "
                  f"avg d={np.mean(fitnesses):.1f}, "
                  f"k={gen_best_code.k}, "
                  f"unique={unique_hashes}/{len(population)}, "
                  f"mut={current_mutation_rate:.3f}")

        # Selection + reproduction
        new_population = []
        new_fitnesses = []
        new_hashes = set()

        # Elitism: keep top codes
        elite_indices = np.argsort(fitnesses)[-elitism:]
        for idx in elite_indices:
            new_population.append(population[idx])
            new_fitnesses.append(fitnesses[idx])
            new_hashes.add(population[idx].code_hash())

        # Fill rest with offspring
        attempts = 0
        max_attempts = pop_size * 5
        while len(new_population) < pop_size and attempts < max_attempts:
            attempts += 1

            # Rank-based selection
            parent1 = _rank_select(population, fitnesses, rng)
            parent2 = _rank_select(population, fitnesses, rng)

            # Crossover: preserve compatible Hz rows
            child = _crossover_css(parent1, parent2, n, rng)
            if child is None:
                # Crossover failed, mutate a parent instead
                child = _mutate_css(parent1, n, current_mutation_rate, rng)

            if child is None:
                continue

            # Row-level mutation
            if rng.random() < 0.8:  # 80% chance of mutation
                mutant = _mutate_css(child, n, current_mutation_rate, rng)
                if mutant is not None:
                    child = mutant

            if child.k < 1:
                continue

            stats.total_codes_tested += 1
            stats.valid_codes += 1
            stats.codes_with_k_ge_1 += 1

            # Diversity: reject exact duplicates
            child_hash = child.code_hash()
            if child_hash in new_hashes:
                continue
            new_hashes.add(child_hash)

            # Use exact distance for small codes
            d_est = exact_distance_css(child, timeout_sec=max(1, distance_timeout // 2))
            if d_est is None:
                d_est = estimate_distance_random(child, num_samples=500, rng=rng)
            if d_est is None:
                d_est = 1
            child.distance = d_est

            new_population.append(child)
            new_fitnesses.append(d_est)

        # If population is too small, inject fresh random codes
        inject_attempts = 0
        while len(new_population) < pop_size and inject_attempts < pop_size * 3:
            inject_attempts += 1
            # Rotate methods for injection diversity
            if inject_attempts % 3 == 0:
                code = random_self_orthogonal(n, target_k=target_k, rng=rng)
            elif inject_attempts % 3 == 1:
                code = random_bicycle_code(n, rng=rng) if n % 2 == 0 else random_css_code(n, target_k=target_k, rng=rng)
            else:
                code = random_css_code(n, target_k=target_k, density=rng.uniform(0.2, 0.5), rng=rng)
            if code is None or code.k < 1:
                continue
            h = code.code_hash()
            if h in new_hashes:
                continue
            new_hashes.add(h)
            stats.total_codes_tested += 1
            stats.valid_codes += 1
            stats.codes_with_k_ge_1 += 1
            d = exact_distance_css(code, timeout_sec=max(1, distance_timeout // 2))
            if d is None:
                d = 1
            code.distance = d
            new_population.append(code)
            new_fitnesses.append(d)

        population = new_population
        fitnesses = new_fitnesses

    stats.elapsed_sec = time.time() - start

    # Final exact distance for best code
    if best_ever_code is not None:
        d_final = exact_distance_css(best_ever_code, timeout_sec=120)
        if d_final is not None:
            best_ever_code.distance = d_final
            key = (best_ever_code.n, best_ever_code.k)
            stats.best_distance[key] = d_final
            print(f"\nFinal best: [[{best_ever_code.n}, {best_ever_code.k}, {d_final}]]")

    return stats


def _rank_select(population, fitnesses, rng):
    """Rank-based selection: probability proportional to rank, not raw fitness.

    This maintains selection pressure while preserving diversity better than
    tournament selection alone.
    """
    n = len(population)
    if n == 0:
        return population[0]
    ranks = np.argsort(np.argsort(fitnesses)).astype(float) + 1  # 1-based ranks
    # Add small noise to break ties
    probs = ranks + rng.uniform(0, 0.1, size=n)
    probs = probs / probs.sum()
    idx = rng.choice(n, p=probs)
    return population[idx]


def _tournament_select(population, fitnesses, k, rng):
    """Tournament selection: pick k random individuals, return best."""
    indices = rng.choice(len(population), size=min(k, len(population)), replace=False)
    best_idx = indices[np.argmax([fitnesses[i] for i in indices])]
    return population[best_idx]


def _crossover_css(parent1, parent2, n, rng):
    """Crossover two CSS codes, preserving compatible Hz rows.

    Strategy:
    1. Mix Hx rows from both parents
    2. Keep Hz rows from both parents that commute with the new Hx
    3. Only generate new random Hz rows if needed to fill gaps
    """
    from codes import gf2_nullspace

    Hx1, Hx2 = parent1.Hx, parent2.Hx

    # Make same number of rows
    max_rows = max(Hx1.shape[0], Hx2.shape[0])
    if Hx1.shape[0] < max_rows:
        pad = np.zeros((max_rows - Hx1.shape[0], n), dtype=np.uint8)
        Hx1 = np.vstack([Hx1, pad])
    if Hx2.shape[0] < max_rows:
        pad = np.zeros((max_rows - Hx2.shape[0], n), dtype=np.uint8)
        Hx2 = np.vstack([Hx2, pad])

    # Uniform crossover on Hx rows
    child_Hx = np.zeros_like(Hx1)
    for i in range(max_rows):
        if rng.random() < 0.5:
            child_Hx[i] = Hx1[i]
        else:
            child_Hx[i] = Hx2[i]

    # Remove zero rows
    child_Hx = child_Hx[np.any(child_Hx, axis=1)]
    if child_Hx.shape[0] == 0:
        return None

    # Try to preserve Hz rows from both parents that commute with new Hx
    all_hz_candidates = np.vstack([parent1.Hz, parent2.Hz])
    # Remove duplicate rows
    unique_rows = set()
    kept_hz = []
    for row in all_hz_candidates:
        row_key = tuple(row)
        if row_key in unique_rows or not np.any(row):
            continue
        unique_rows.add(row_key)
        # Check commutativity: row @ child_Hx.T = 0 (mod 2)
        product = (child_Hx @ row.reshape(-1, 1)) % 2
        if not np.any(product):
            kept_hz.append(row)

    # Target number of Hz rows
    target_rz = (parent1.Hz.shape[0] + parent2.Hz.shape[0]) // 2
    target_rz = max(target_rz, 1)

    if len(kept_hz) >= target_rz:
        # We have enough compatible rows — use them
        child_Hz = np.array(kept_hz[:target_rz], dtype=np.uint8)
    else:
        # Need more rows — fill from nullspace
        null_hx = gf2_nullspace(child_Hx)
        if null_hx.shape[0] == 0 and len(kept_hz) == 0:
            return None

        if len(kept_hz) > 0:
            child_Hz_list = list(kept_hz)
        else:
            child_Hz_list = []

        if null_hx.shape[0] > 0:
            needed = target_rz - len(child_Hz_list)
            needed = min(needed, null_hx.shape[0])
            for _ in range(max(needed, 1)):
                coeffs = rng.integers(0, 2, size=null_hx.shape[0])
                if not np.any(coeffs):
                    coeffs[rng.integers(0, null_hx.shape[0])] = 1
                new_row = (coeffs @ null_hx) % 2
                if np.any(new_row):
                    child_Hz_list.append(new_row)

        if len(child_Hz_list) == 0:
            return None
        child_Hz = np.array(child_Hz_list, dtype=np.uint8)

    # Remove zero rows
    child_Hz = child_Hz[np.any(child_Hz, axis=1)]
    if child_Hz.shape[0] == 0:
        return None

    code = CSSCode(Hx=child_Hx, Hz=child_Hz)
    valid, _ = code.is_valid()
    return code if valid else None


def _mutate_css(code, n, rate, rng):
    """Row-level mutation of a CSS code.

    Instead of flipping individual bits (too destructive), perform
    structural row operations:
    1. Swap two rows within Hx or Hz
    2. Add one row to another (mod 2) — preserves commutativity within Hx/Hz
    3. Replace a row with a random vector from the nullspace

    These operations are much less likely to break the code.
    """
    from codes import gf2_nullspace

    Hx = code.Hx.copy()
    Hz = code.Hz.copy()

    # Pick a random row-level operation
    op = rng.choice(5)

    if op == 0 and Hx.shape[0] >= 2:
        # Swap two Hx rows (always valid, just reordering)
        i, j = rng.choice(Hx.shape[0], size=2, replace=False)
        Hx[[i, j]] = Hx[[j, i]]
        # After swapping Hx rows, Hz is still valid (commutativity unchanged)

    elif op == 1 and Hx.shape[0] >= 2:
        # Add one Hx row to another (mod 2) — preserves commutativity
        i, j = rng.choice(Hx.shape[0], size=2, replace=False)
        Hx[i] = (Hx[i] + Hx[j]) % 2

    elif op == 2 and Hz.shape[0] >= 2:
        # Add one Hz row to another (mod 2) — preserves commutativity
        i, j = rng.choice(Hz.shape[0], size=2, replace=False)
        Hz[i] = (Hz[i] + Hz[j]) % 2

    elif op == 3:
        # Replace one Hx row with a row from nullspace of Hz
        null_hz = gf2_nullspace(Hz)
        if null_hz.shape[0] > 0:
            i = rng.integers(0, Hx.shape[0])
            coeffs = rng.integers(0, 2, size=null_hz.shape[0])
            if not np.any(coeffs):
                coeffs[rng.integers(0, null_hz.shape[0])] = 1
            new_row = (coeffs @ null_hz) % 2
            if np.any(new_row):
                Hx[i] = new_row
                # This new Hx row is in nullspace of Hz, so commutativity holds

    elif op == 4:
        # Replace one Hz row with a row from nullspace of Hx
        null_hx = gf2_nullspace(Hx)
        if null_hx.shape[0] > 0:
            i = rng.integers(0, Hz.shape[0])
            coeffs = rng.integers(0, 2, size=null_hx.shape[0])
            if not np.any(coeffs):
                coeffs[rng.integers(0, null_hx.shape[0])] = 1
            new_row = (coeffs @ null_hx) % 2
            if np.any(new_row):
                Hz[i] = new_row

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


def algebraic_search(max_n=30, callback=None):
    """Search for codes via hypergraph products of small classical codes.

    Systematically constructs HGP codes from all pairs of small
    classical parity check matrices.
    """
    stats = SearchStats(strategy="algebraic")
    start = time.time()

    # Small classical parity check matrices
    classical_codes = []

    # Repetition codes
    for length in range(3, 8):
        H = np.zeros((length - 1, length), dtype=np.uint8)
        for i in range(length - 1):
            H[i, i] = 1
            H[i, i + 1] = 1
        classical_codes.append((f"rep({length})", H))

    # Hamming codes
    for r in range(2, 5):
        n_ham = 2**r - 1
        H = np.zeros((r, n_ham), dtype=np.uint8)
        for j in range(n_ham):
            col = j + 1
            for i in range(r):
                H[i, j] = (col >> i) & 1
        classical_codes.append((f"ham({r})", H))

    # Random small codes
    rng = np.random.default_rng(42)
    for _ in range(10):
        n_cls = rng.integers(3, 7)
        r_cls = rng.integers(1, n_cls)
        H = rng.integers(0, 2, size=(r_cls, n_cls)).astype(np.uint8)
        H = H[np.any(H, axis=1)]
        if H.shape[0] > 0:
            classical_codes.append((f"rand({r_cls}x{n_cls})", H))

    print(f"Testing {len(classical_codes)} classical codes in HGP...")

    for i, (name1, H1) in enumerate(classical_codes):
        for j, (name2, H2) in enumerate(classical_codes):
            if j < i:
                continue  # Symmetric

            code = hypergraph_product(H1, H2)
            stats.total_codes_tested += 1

            if code.n > max_n:
                continue

            valid, _ = code.is_valid()
            if not valid:
                continue
            stats.valid_codes += 1

            if code.k < 1:
                continue
            stats.codes_with_k_ge_1 += 1

            # Compute distance
            d = exact_distance_css(code, timeout_sec=60)
            if d is None or d < 1:
                continue

            code.distance = d
            key = (code.n, code.k)

            beats = False
            if key in KNOWN_BOUNDS:
                known_lower, _ = KNOWN_BOUNDS[key]
                if d > known_lower:
                    beats = True

            result = SearchResult(
                code=code, n=code.n, k=code.k, d=d,
                strategy=f"hgp({name1}x{name2})", beats_known=beats,
            )

            if key not in stats.best_distance or d > stats.best_distance[key]:
                stats.best_distance[key] = d
                stats.best_codes.append(result)
                print(f"  HGP({name1} x {name2}): [[{code.n}, {code.k}, {d}]]"
                      f"{' *** BEATS KNOWN' if beats else ''}")

            if callback:
                callback(stats.total_codes_tested, code, d)

    stats.elapsed_sec = time.time() - start
    return stats
