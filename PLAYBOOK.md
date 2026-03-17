# Quantum Error Correction Code Search

Systematic search for stabilizer codes with parameters [[n, k, d]] that match or improve on
best-known codes. New codes with better distance are publishable results.

**Owner:** Matt Loftus | **Dir:** `~/workspace/qec-search`

---

## Mission

Quantum error correction is THE bottleneck for practical quantum computing. The space of
possible stabilizer codes is enormous and under-explored. This project searches that space
computationally — using genetic algorithms, random search, and algebraic constructions — to
find codes with better parameters than currently known.

**Success criteria:** Find a stabilizer code [[n, k, d]] that improves on known lower bounds
in Grassl's codetables.de, or demonstrate a new construction technique.

---

## Stack

- **Python 3.9** (`/usr/bin/python3`) — search engine
- **stim** — fast stabilizer simulation, approximate distance (graphlike)
- **ldpc** — GF(2) linear algebra, exact distance, BP+OSD decoding
- **numpy / scipy** — matrix operations
- **React 19 + TypeScript + Vite + Zustand + Plotly + Tailwind** — web dashboard

---

## Project Structure

```
qec-search/
├── PLAYBOOK.md              # This file
├── engine/
│   ├── codes.py             # Stabilizer code representation (CSS, general)
│   ├── distance.py          # Distance computation (exact + approximate)
│   ├── validate.py          # Code property verification
│   ├── search.py            # Search strategies (random, genetic)
│   ├── benchmark.py         # Comparison against known code tables
│   ├── simulate.py          # Logical error rate Monte Carlo
│   ├── known_codes.py       # Reference implementations of known codes
│   ├── serve.py             # REST API for dashboard
│   └── evolver.py           # Main CLI entry point
├── web/                     # React dashboard
│   ├── src/
│   │   ├── App.tsx
│   │   ├── store.ts
│   │   ├── api.ts
│   │   ├── types.ts
│   │   └── components/
│   │       ├── Guide.tsx
│   │       ├── StatsPanel.tsx
│   │       ├── SearchProgress.tsx
│   │       ├── CodeViewer.tsx
│   │       ├── DistancePlot.tsx
│   │       └── BenchmarkTable.tsx
│   └── package.json
├── results/                 # Search campaign outputs
│   ├── campaigns/           # Per-campaign results
│   ├── best_codes.json      # All-time best codes found
│   └── benchmark.json       # Comparison vs codetables.de
└── scripts/
    └── campaign.sh          # Long-running search automation
```

---

## How to Run

```bash
cd ~/workspace/qec-search

# Validate against known codes (first thing to run)
/usr/bin/python3 engine/evolver.py --validate

# Random search for CSS codes (quick exploration)
/usr/bin/python3 engine/evolver.py --mode random --n 15 --k 1 --samples 10000

# Genetic algorithm search
/usr/bin/python3 engine/evolver.py --mode genetic --n 15 --k 1 --pop 200 --gens 500

# Search a range of parameters
/usr/bin/python3 engine/evolver.py --mode genetic --n-range 10,20 --k 1

# Compare results against known tables
/usr/bin/python3 engine/evolver.py --benchmark

# Simulate logical error rate for a found code
/usr/bin/python3 engine/evolver.py --simulate --code results/best_codes.json --index 0

# Start web dashboard API
/usr/bin/python3 engine/serve.py --port 8789

# Start web dashboard
cd web && npm run dev
```

---

## Stabilizer Codes: Quick Reference

An **[[n, k, d]]** stabilizer code encodes **k** logical qubits into **n** physical qubits
with minimum distance **d**. The code can correct floor((d-1)/2) arbitrary errors.

**Stabilizer group:** n-k independent commuting Pauli operators (generators).
The codespace is the simultaneous +1 eigenspace of all generators.

**CSS codes** are a subclass where generators split into pure-X and pure-Z types:
- X-checks: `Hx` (each row is a Z-type stabilizer measuring X errors)
- Z-checks: `Hz` (each row is an X-type stabilizer measuring Z errors)
- Commutativity constraint: `Hx @ Hz.T = 0 (mod 2)`
- Parameters: n = number of columns, k = n - rank(Hx) - rank(Hz)

**Distance:** The minimum weight of a nontrivial logical operator (Pauli that commutes
with all stabilizers but is not itself a stabilizer). Computing this is NP-hard in general,
but feasible for small codes by exhaustive search.

---

## Search Strategies

### 1. Random CSS Search
Generate random binary matrices Hx, Hz satisfying the commutativity constraint.
Simple but establishes baselines.

### 2. Genetic Algorithm
- **Genome:** Binary parity check matrix (flattened)
- **Fitness:** Code distance (primary), weight of stabilizers (secondary — lower is better for hardware)
- **Crossover:** Row exchange between parents
- **Mutation:** Bit flips in parity check matrix (must maintain commutativity)
- **Selection:** Tournament selection with elitism

### 3. Algebraic Constructions
- Hypergraph products of classical codes
- Bicycle and bivariate bicycle codes
- Lifted product codes
- Seed with algebraic construction, refine with local search

---

## Benchmarking

Reference: Markus Grassl's **codetables.de** — definitive database of best-known [[n,k,d]].

For each (n, k), the table gives lower and upper bounds on achievable distance d.
- If our code's d **exceeds the lower bound**: new best-known code (publishable)
- If our code's d **matches the lower bound**: confirms known result
- If our code's d **exceeds the upper bound**: we have a bug

The known bounds for small qubit codes (q=2) are stored in `engine/known_codes.py`.

---

## Validated Reference Codes

| Code | [[n, k, d]] | Type | Notes |
|------|-------------|------|-------|
| Perfect | [[5, 1, 3]] | General | Smallest single-error-correcting |
| Steane | [[7, 1, 3]] | CSS | Smallest CSS, transversal T |
| Shor | [[9, 1, 3]] | CSS | First QEC code, repetition-based |
| Surface-3 | [[9, 1, 3]] | CSS | Rotated surface code, d=3 |
| Surface-5 | [[25, 1, 5]] | CSS | Rotated surface code, d=5 |
| Reed-Muller | [[15, 1, 3]] | CSS | Punctured Reed-Muller |
| Hamming-CSS | [[15, 7, 3]] | CSS | From classical Hamming |
| Best known | [[11, 1, 5]] | General | Smallest 2-error-correcting |

---

## Milestones

- [x] Engine: code representation + validation
- [x] Engine: distance computation (exact, optimized)
- [x] Engine: validate known reference codes (4/4 pass)
- [x] Engine: random search (4 constructors: basic CSS, classical, self-orthogonal, bicycle)
- [x] Engine: genetic algorithm search (row-level mutation, adaptive, diversity-preserving)
- [x] Engine: algebraic search (hypergraph products)
- [x] Engine: benchmark against codetables.de (58 bounds, 33 open gaps)
- [x] Engine: Monte Carlo error correction simulation
- [x] Web: dashboard with stats, distance plot, best codes table
- [x] Web: Campaigns tab
- [x] Web: Bounds tab with gap highlighting
- [x] Web: Guide tab
- [x] First search campaigns (6 campaigns, 59K+ codes tested)
- [x] First code matching known bounds ([[9,1,3]] — Shor-equivalent via self-orthogonal construction)
- [x] Integrate qLDPC library (Python 3.12 venv, BBCode + HGPCode + QCCode + distance computation)
- [x] Test published Bravyi et al. BB codes (params match paper, distance too expensive for n>=72)
- [x] Systematic polynomial sweep: BB degree-3, QCCode weight-2/3 (consistently 40-60% of known bounds)
- [x] QCCode (generalized bicycle) search — best family: [[48,4,8]], [[42,6,6]], [[30,4,6]]
- [x] Tested published Bravyi/Wang-Mueller polynomials — confirmed [[54,8,6]] Wang-Mueller code
- [x] Nightly automation (launchd at 2 AM, ntfy notifications)
- [x] Nightly automation (launchd at 2 AM, ntfy, nice -n 15)
- [x] Expand known bounds table (154 entries, n=5-60, corrected from codetables.de)
- [x] Stim-based approximate distance (Monte Carlo, works for n=72-144+, verified on Bravyi codes)
- [x] Genetic search over polynomial exponents (evolve_polys.py — weight-3 maxes at d=8)

**Novel approaches:**
- [ ] ML-guided search: train model on (polynomial features → distance), score 100K candidates, compute exact distance on top predictions
- [ ] Puncture/shorten published large codes: systematically reduce Bravyi [[288,12,18]] and [[784,24,24]] to find good subcodes
- [ ] Group-space search: fix Bravyi polynomials, vary the group (dihedral, quaternion, semidirect products via LPCode)
- [ ] Reverse-engineer Bravyi polynomial algebra: factor over GF(2) extensions, find other polynomials with same algebraic fingerprint
- [ ] Cover code construction: find 3-covers of good codes, compute covers of non-Bravyi codes
- [ ] RL with stim distance as fast oracle: bandit/policy gradient over polynomial exponents
- [ ] Code concatenation/product search: tensor products of best small codes
- [ ] Weight-4 evolutionary search (evolve_polys.py --weight 4)

- [ ] First code exceeding known bounds

---

## Roadmap

### Near-term
- **ML-guided search** — train model on existing (polynomial, distance) data to predict promising candidates. Score 100K polynomials, compute exact distance on top 100. Highest expected value approach.
- **Puncture/shorten published codes** — Bravyi [[288,12,18]] and [[784,24,24]] contain good subcodes. Systematically puncture (remove qubits) and shorten (fix qubits). Use stim distance to evaluate.
- **Weight-4 evolutionary search** — current weight-3 maxes at d=8. Weight-4 polynomials are unexplored and could break through.

### Medium-term
- **Group-space search** — fix Bravyi polynomials, vary the group structure (dihedral, quaternion, semidirect products) via qLDPC's LPCode
- **Reverse-engineer Bravyi polynomial algebra** — factor over GF(2) extensions, find the algebraic fingerprint that makes (x³+y+y², y³+x+x²) special, search for other polynomials with same properties
- **Cover code construction** — the [[144,12,12]] is a 2-cover of [[72,12,6]]. Find 3-covers, or covers of other good codes
- **RL with stim distance oracle** — bandit/policy gradient over polynomial exponents, stim gives feedback in 10-40s
- **Code concatenation/product search** — tensor products of best small codes

### Long-term
- Non-CSS stabilizer codes (general symplectic representation)
- If novel code found: document for publication, contact quantum error correction community

---

## Lessons Learned

- Random CSS codes have 70% d=1 at n=10 — most randomly generated codes are useless
- Bicycle (circulant) codes are the best random constructor: 50% hit d=3 for n=10, always produce k=2
- Self-orthogonal codes can't produce k=1 at small n (constraints too tight)
- Self-orthogonal codes found [[16,4,4]] matching known best — good for k>1
- Known bounds for small codes (n<20) are very tight — random methods can't close ANY of the 33 open gaps even with 50K+ samples per gap
- Genetic search needs bicycle/self-orthogonal in initialization to work (pure random init → d=2 ceiling)
- Row-level mutation (swap rows, add rows mod 2, replace from nullspace) vastly outperforms bit-level
- Distance computation is the bottleneck (5-20x slower than code generation) — optimized with batch enumeration + RREF precomputation
- HGP of rep(4)×rep(4) recovers the surface code [[25,1,4]]
- Extending/puncturing known good codes produces nothing — the algebraic structure is too fragile
- qLDPC `get_distance()` is fast for n < ~40 but exponential — needs timeout for larger codes
- qLDPC tries to call GAP (algebra system) for automorphisms — noisy but harmless
- BB([3,m], y+1, xy+1) is a productive polynomial family: d=5 at n=30, d=6 at n=36-54, d=7 at n=60
- BB codes need Python 3.12 venv: `venv/bin/python3 engine/qldpc_search.py`
- HGP codes from qLDPC match our own implementation but with faster distance computation
- QCCode (univariate generalized bicycle) is the most productive code family in qLDPC — fast distance computation, diverse parameter coverage
- Published Bravyi et al. codes use (x^3+y+y^2, y^3+x+x^2) — only works when l,m are multiples of 3
- Published codes at n>=72 are beyond exact distance computation — need stim-based methods or get_distance(bound=N)
- Systematic sweep of weight-3 degree-3 polynomials consistently finds d at 40-60% of known lower bounds — higher degree/more terms needed
- The known bounds for n>=20 come from algebraic constructions over GF(4) computed in MAGMA — these are computationally intensive codes that random polynomial search can't reproduce
- Key realization: closing gaps likely requires either (a) much larger polynomial search space, (b) non-cyclic group structures, or (c) fundamentally new construction techniques
- Stim Monte Carlo distance: exact for Bravyi [[72,12,6]] (d=6 in 17s) and [[144,12,12]] (d=12 in 41s)
- Graphlike distance returns None/inf for BB/QC codes (high-weight stabilizers) — must use MC fallback
- Stim MC works by sampling random errors, checking if logical operators are hit, extracting min-weight from hits
- Verified all published Bravyi and Wang-Mueller distances: [[108,8,10]], [[90,8,10]], [[98,6,12]] all correct
- QC(30) with our best polynomials gives [[60,4,8]] — d=8 at n=60, known bound is d=14
- evolve_polys.py: genetic search over polynomial exponents, seeded from best known codes
- Weight-4 evolve found [[20,2,6]] MATCHING known bounds — then exhaustive search proved [[20,2,7]] doesn't exist as a CSS code (360K+ codes tested across QC, BB, and random CSS)
- Algebraic analysis (algebraic_analysis.py): Bravyi polynomials are UNITS in the group ring (0% root vanishing). This requires l,m divisible by 3 where x^l-1 factors into only (x+1) and (x²+x+1)
- Being a unit is necessary but not sufficient — (4,4), (5,5) have all-unit polynomials but produce k=0 codes. The relationship between A and B matters.
- Symmetric polynomial pairs B(x,y)=A(y,x) ensure dx=dz (equal X/Z distance)
- Wang-Mueller codes tolerate 12-22% vanishing but compensate with more CRT components

---

*Last updated: 2026-03-15*
