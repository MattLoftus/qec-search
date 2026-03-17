export function Guide() {
  return (
    <div className="max-w-4xl mx-auto space-y-8 p-6">
      <h1 className="text-2xl font-bold text-[var(--color-text-primary)]">
        QEC Code Search — Guide
      </h1>
      <p className="text-[var(--color-text-secondary)]">
        A systematic search for quantum error correcting codes with parameters
        that match or improve on the best known codes. New codes with better
        distance are publishable results.
      </p>

      {/* Overview */}
      <section>
        <h2 className="text-lg font-semibold text-[var(--color-accent-blue)] mb-3">
          Overview
        </h2>
        <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
          <div className="bg-white/[0.02] border border-[var(--color-border-dim)] rounded-lg p-4">
            <h3 className="font-medium text-[var(--color-text-primary)] mb-2">
              Code Construction
            </h3>
            <p className="text-sm text-[var(--color-text-secondary)]">
              Generate candidate stabilizer codes using bivariate bicycle codes
              (BBCode), generalized bicycle codes (QCCode), evolutionary polynomial
              search, hypergraph products, and random CSS construction. Powered by
              qLDPC (Google/IBM's quantum code library).
            </p>
          </div>
          <div className="bg-white/[0.02] border border-[var(--color-border-dim)] rounded-lg p-4">
            <h3 className="font-medium text-[var(--color-text-primary)] mb-2">
              Distance Computation
            </h3>
            <p className="text-sm text-[var(--color-text-secondary)]">
              Compute the code distance d — the minimum weight of a nontrivial
              logical operator. Three tiers: exact (n &lt; 40 via qLDPC), stim Monte
              Carlo (n &lt; 200 via BP+OSD decoding), and graphlike upper bound.
              Verified exact on Bravyi [[144,12,12]] in 41 seconds.
            </p>
          </div>
          <div className="bg-white/[0.02] border border-[var(--color-border-dim)] rounded-lg p-4">
            <h3 className="font-medium text-[var(--color-text-primary)] mb-2">
              Benchmarking
            </h3>
            <p className="text-sm text-[var(--color-text-secondary)]">
              Compare found codes against Markus Grassl's codetables.de — the
              definitive database of best-known quantum codes. A code that
              exceeds the known lower bound is a genuine new result.
            </p>
          </div>
        </div>
      </section>

      {/* Stabilizer Codes */}
      <section>
        <h2 className="text-lg font-semibold text-[var(--color-accent-blue)] mb-3">
          Stabilizer Codes
        </h2>
        <div className="bg-white/[0.02] border border-[var(--color-border-dim)] rounded-lg p-4 space-y-3">
          <p className="text-sm text-[var(--color-text-secondary)]">
            An <strong className="text-[var(--color-text-primary)]">[[n, k, d]]</strong> stabilizer
            code encodes <strong>k</strong> logical qubits into <strong>n</strong> physical
            qubits with minimum distance <strong>d</strong>. The code can correct
            ⌊(d-1)/2⌋ arbitrary single-qubit errors.
          </p>
          <p className="text-sm text-[var(--color-text-secondary)]">
            The code is defined by its <strong className="text-[var(--color-text-primary)]">stabilizer group</strong>:
            n-k independent commuting Pauli operators (generators). The codespace is the
            simultaneous +1 eigenspace of all generators.
          </p>
          <p className="text-sm text-[var(--color-text-secondary)]">
            <strong className="text-[var(--color-text-primary)]">CSS codes</strong> are a subclass
            where generators split into pure-X and pure-Z types, represented by binary parity
            check matrices Hx (detects X errors) and Hz (detects Z errors). The commutativity
            constraint Hx·Hz<sup>T</sup> = 0 (mod 2) ensures the generators commute.
          </p>
        </div>
      </section>

      {/* Key Parameters */}
      <section>
        <h2 className="text-lg font-semibold text-[var(--color-accent-blue)] mb-3">
          Key Parameters
        </h2>
        <div className="overflow-x-auto">
          <table className="w-full text-sm">
            <thead>
              <tr className="border-b border-[var(--color-border-dim)]">
                <th className="text-left p-2 text-[var(--color-text-primary)]">Parameter</th>
                <th className="text-left p-2 text-[var(--color-text-primary)]">Meaning</th>
                <th className="text-left p-2 text-[var(--color-text-primary)]">Good</th>
                <th className="text-left p-2 text-[var(--color-text-primary)]">Excellent</th>
              </tr>
            </thead>
            <tbody className="text-[var(--color-text-secondary)]">
              <tr className="border-b border-[var(--color-border-dim)]">
                <td className="p-2 font-mono">n</td>
                <td className="p-2">Physical qubits</td>
                <td className="p-2">10-20</td>
                <td className="p-2">5-10 (for given k,d)</td>
              </tr>
              <tr className="border-b border-[var(--color-border-dim)]">
                <td className="p-2 font-mono">k</td>
                <td className="p-2">Logical qubits encoded</td>
                <td className="p-2">&ge; 1</td>
                <td className="p-2">High k/n ratio</td>
              </tr>
              <tr className="border-b border-[var(--color-border-dim)]">
                <td className="p-2 font-mono">d</td>
                <td className="p-2">Code distance (errors correctable = ⌊(d-1)/2⌋)</td>
                <td className="p-2">3-5</td>
                <td className="p-2">&ge; 7</td>
              </tr>
              <tr className="border-b border-[var(--color-border-dim)]">
                <td className="p-2 font-mono">k/n</td>
                <td className="p-2">Code rate (encoding efficiency)</td>
                <td className="p-2">&gt; 0.1</td>
                <td className="p-2">&gt; 0.3</td>
              </tr>
              <tr>
                <td className="p-2 font-mono">weight</td>
                <td className="p-2">Max stabilizer weight (hardware feasibility)</td>
                <td className="p-2">&le; 6</td>
                <td className="p-2">&le; 4</td>
              </tr>
            </tbody>
          </table>
        </div>
      </section>

      {/* Search Strategies */}
      <section>
        <h2 className="text-lg font-semibold text-[var(--color-accent-blue)] mb-3">
          Search Strategies
        </h2>
        <div className="space-y-3">
          <div className="bg-white/[0.02] border border-[var(--color-border-dim)] rounded-lg p-4">
            <h3 className="font-medium text-[var(--color-accent-green)] mb-1">Evolutionary Polynomial Search</h3>
            <p className="text-sm text-[var(--color-text-secondary)]">
              Evolve QCCode polynomial exponents via genetic algorithm. The genome is
              (order, [a1,a2,...], [b1,b2,...]) defining polynomials 1+x<sup>a1</sup>+x<sup>a2</sup>+...
              in the cyclic group ring. Supports weight-3 and weight-4 polynomials. Found
              [[20,2,6]] matching known bounds. Runs nightly at 2 AM.
            </p>
          </div>
          <div className="bg-white/[0.02] border border-[var(--color-border-dim)] rounded-lg p-4">
            <h3 className="font-medium text-[var(--color-accent-green)] mb-1">Bivariate Bicycle Codes (BBCode)</h3>
            <p className="text-sm text-[var(--color-text-secondary)]">
              Construct codes from two polynomials in the bivariate group ring
              F<sub>2</sub>[Z<sub>l</sub> x Z<sub>m</sub>]. The Bravyi et al. family uses
              A=x<sup>3</sup>+y+y<sup>2</sup>, B=y<sup>3</sup>+x+x<sup>2</sup> to produce codes like
              [[72,12,6]] and [[144,12,12]]. These polynomials are units (invertible) in
              the group ring — a key algebraic property.
            </p>
          </div>
          <div className="bg-white/[0.02] border border-[var(--color-border-dim)] rounded-lg p-4">
            <h3 className="font-medium text-[var(--color-accent-green)] mb-1">Generalized Bicycle (QCCode)</h3>
            <p className="text-sm text-[var(--color-text-secondary)]">
              Univariate variant of BB codes over a single cyclic group Z<sub>n</sub>. Faster distance
              computation and the most productive code family overall. Best result: [[48,4,8]].
              Weight-3 polynomials max at d=8; weight-4 reaches d=9 and found the [[20,2,6]] match.
            </p>
          </div>
          <div className="bg-white/[0.02] border border-[var(--color-border-dim)] rounded-lg p-4">
            <h3 className="font-medium text-[var(--color-accent-green)] mb-1">Hypergraph Product (HGP)</h3>
            <p className="text-sm text-[var(--color-text-secondary)]">
              Construct quantum codes from pairs of classical codes. HGP(rep(5),rep(5)) produces
              the [[41,1,5]] code. Limited distance ceiling (d=3-5 for n&le;60) but fast
              construction and guaranteed CSS structure.
            </p>
          </div>
          <div className="bg-white/[0.02] border border-[var(--color-border-dim)] rounded-lg p-4">
            <h3 className="font-medium text-[var(--color-accent-green)] mb-1">Algebraic Analysis</h3>
            <p className="text-sm text-[var(--color-text-secondary)]">
              Reverse-engineering approach: factor x<sup>l</sup>-1 and y<sup>m</sup>-1 over GF(2),
              decompose the group ring via CRT, evaluate polynomial vanishing at roots. Good
              codes correspond to polynomial pairs that are near-units (vanish at few/no roots).
              Used to explain why Bravyi polynomials work and guide the search.
            </p>
          </div>
        </div>
      </section>

      {/* Distance Computation */}
      <section>
        <h2 className="text-lg font-semibold text-[var(--color-accent-blue)] mb-3">
          Distance Computation
        </h2>
        <div className="bg-white/[0.02] border border-[var(--color-border-dim)] rounded-lg p-4 space-y-3">
          <p className="text-sm text-[var(--color-text-secondary)]">
            The code distance is the minimum weight of a nontrivial logical operator — a Pauli
            that commutes with all stabilizers but is not itself a stabilizer. For CSS codes, the
            X-distance and Z-distance can be computed independently:
          </p>
          <ul className="text-sm text-[var(--color-text-secondary)] list-disc list-inside space-y-1">
            <li>
              <strong>d<sub>x</sub></strong> = min weight of vector in ker(Hz) \ rowspace(Hx)
            </li>
            <li>
              <strong>d<sub>z</sub></strong> = min weight of vector in ker(Hx) \ rowspace(Hz)
            </li>
            <li>
              <strong>d</strong> = min(d<sub>x</sub>, d<sub>z</sub>)
            </li>
          </ul>
          <p className="text-sm text-[var(--color-text-secondary)]">
            Three computation tiers:
          </p>
          <ul className="text-sm text-[var(--color-text-secondary)] list-disc list-inside space-y-1">
            <li><strong>Exact (n &lt; 40)</strong>: qLDPC&apos;s get_distance() — nullspace enumeration</li>
            <li><strong>Stim Monte Carlo (n &lt; 200)</strong>: Sample errors, decode with BP+OSD, find minimum-weight logical residuals. Verified exact on Bravyi [[72,12,6]] (d=6, 17s) and [[144,12,12]] (d=12, 41s)</li>
            <li><strong>Graphlike bound</strong>: Stim&apos;s shortest_graphlike_error — fast but only works for low-weight stabilizers (not BB/QC codes)</li>
          </ul>
        </div>
      </section>

      {/* Benchmarking */}
      <section>
        <h2 className="text-lg font-semibold text-[var(--color-accent-blue)] mb-3">
          Benchmarking Against Known Codes
        </h2>
        <div className="bg-white/[0.02] border border-[var(--color-border-dim)] rounded-lg p-4 space-y-3">
          <p className="text-sm text-[var(--color-text-secondary)]">
            Results are compared against <strong className="text-[var(--color-text-primary)]">codetables.de</strong> (Markus
            Grassl), the definitive reference for best-known quantum codes. For each (n, k), the
            table provides lower and upper bounds on achievable distance d.
          </p>
          <ul className="text-sm text-[var(--color-text-secondary)] list-disc list-inside space-y-1">
            <li>
              <span className="text-[var(--color-accent-green)]">IMPROVEMENT</span>: our d exceeds the known lower bound — a new best-known code
            </li>
            <li>
              <span className="text-[var(--color-accent-blue)]">MATCHES BEST</span>: our d equals the known lower bound — confirms known result
            </li>
            <li>
              <span className="text-[var(--color-accent-amber)]">BELOW BEST</span>: our d is below the known lower bound — room for improvement
            </li>
            <li>
              <span className="text-[var(--color-text-muted)]">NO REFERENCE</span>: no entry in codetables.de for these parameters
            </li>
          </ul>
        </div>
      </section>

      {/* Famous Codes */}
      <section>
        <h2 className="text-lg font-semibold text-[var(--color-accent-blue)] mb-3">
          Reference Codes
        </h2>
        <div className="overflow-x-auto">
          <table className="w-full text-sm">
            <thead>
              <tr className="border-b border-[var(--color-border-dim)]">
                <th className="text-left p-2 text-[var(--color-text-primary)]">Name</th>
                <th className="text-left p-2 text-[var(--color-text-primary)]">Parameters</th>
                <th className="text-left p-2 text-[var(--color-text-primary)]">Type</th>
                <th className="text-left p-2 text-[var(--color-text-primary)]">Significance</th>
              </tr>
            </thead>
            <tbody className="text-[var(--color-text-secondary)]">
              <tr className="border-b border-[var(--color-border-dim)]">
                <td className="p-2">Perfect Code</td>
                <td className="p-2 font-mono">[[5, 1, 3]]</td>
                <td className="p-2">General</td>
                <td className="p-2">Smallest code correcting 1 error</td>
              </tr>
              <tr className="border-b border-[var(--color-border-dim)]">
                <td className="p-2">Steane Code</td>
                <td className="p-2 font-mono">[[7, 1, 3]]</td>
                <td className="p-2">CSS</td>
                <td className="p-2">Smallest CSS code, transversal gates</td>
              </tr>
              <tr className="border-b border-[var(--color-border-dim)]">
                <td className="p-2">Shor Code</td>
                <td className="p-2 font-mono">[[9, 1, 3]]</td>
                <td className="p-2">CSS</td>
                <td className="p-2">First QEC code ever discovered</td>
              </tr>
              <tr className="border-b border-[var(--color-border-dim)]">
                <td className="p-2">Surface Code (d=3)</td>
                <td className="p-2 font-mono">[[9, 1, 3]]</td>
                <td className="p-2">CSS</td>
                <td className="p-2">Leading candidate for fault-tolerant QC</td>
              </tr>
              <tr className="border-b border-[var(--color-border-dim)]">
                <td className="p-2">[[11,1,5]] Code</td>
                <td className="p-2 font-mono">[[11, 1, 5]]</td>
                <td className="p-2">General</td>
                <td className="p-2">Smallest code correcting 2 errors</td>
              </tr>
              <tr className="border-b border-[var(--color-border-dim)]">
                <td className="p-2">Hamming CSS</td>
                <td className="p-2 font-mono">[[15, 7, 3]]</td>
                <td className="p-2">CSS</td>
                <td className="p-2">High-rate CSS from classical Hamming</td>
              </tr>
              <tr className="border-b border-[var(--color-border-dim)]">
                <td className="p-2">Bravyi BB</td>
                <td className="p-2 font-mono">[[72, 12, 6]]</td>
                <td className="p-2">BB/CSS</td>
                <td className="p-2">Smallest code in the Bravyi family</td>
              </tr>
              <tr className="border-b border-[var(--color-border-dim)]">
                <td className="p-2">Gross Code</td>
                <td className="p-2 font-mono">[[144, 12, 12]]</td>
                <td className="p-2">BB/CSS</td>
                <td className="p-2">Bravyi et al. 2024 — high-threshold QEC</td>
              </tr>
              <tr>
                <td className="p-2">Wang-Mueller</td>
                <td className="p-2 font-mono">[[98, 6, 12]]</td>
                <td className="p-2">BB/CSS</td>
                <td className="p-2">Coprime BB construction, d=12 at n=98</td>
              </tr>
            </tbody>
          </table>
        </div>
      </section>

      {/* How to Use */}
      <section>
        <h2 className="text-lg font-semibold text-[var(--color-accent-blue)] mb-3">
          Reading the Dashboard
        </h2>
        <div className="space-y-3">
          <div className="bg-white/[0.02] border border-[var(--color-border-dim)] rounded-lg p-4">
            <h3 className="font-medium text-[var(--color-text-primary)] mb-1">Stats Panel</h3>
            <p className="text-sm text-[var(--color-text-secondary)]">
              Top-level KPIs: total codes tested, best codes found, benchmark status.
              Color coding: green = matches/exceeds known bounds, amber = below, blue = informational.
            </p>
          </div>
          <div className="bg-white/[0.02] border border-[var(--color-border-dim)] rounded-lg p-4">
            <h3 className="font-medium text-[var(--color-text-primary)] mb-1">Best Codes Table</h3>
            <p className="text-sm text-[var(--color-text-secondary)]">
              All-time best codes found per (n, k) parameter pair. Shows distance, rate,
              stabilizer weight, and comparison against known bounds. Click a code to see
              its full parity check matrices.
            </p>
          </div>
          <div className="bg-white/[0.02] border border-[var(--color-border-dim)] rounded-lg p-4">
            <h3 className="font-medium text-[var(--color-text-primary)] mb-1">Campaigns Tab</h3>
            <p className="text-sm text-[var(--color-text-secondary)]">
              History of search campaigns — random, genetic, algebraic. Shows codes tested,
              best results, and elapsed time for each run. Use this to plan which parameter
              regimes need more search effort.
            </p>
          </div>
          <div className="bg-white/[0.02] border border-[var(--color-border-dim)] rounded-lg p-4">
            <h3 className="font-medium text-[var(--color-text-primary)] mb-1">Bounds Tab</h3>
            <p className="text-sm text-[var(--color-text-secondary)]">
              Full table of known bounds from codetables.de. Rows highlighted in amber have
              a gap between lower and upper bounds — these are open problems where finding a
              better code would be a publishable result.
            </p>
          </div>
        </div>
      </section>

      {/* Technical Stack */}
      <section>
        <h2 className="text-lg font-semibold text-[var(--color-accent-blue)] mb-3">
          Technical Stack
        </h2>
        <div className="bg-white/[0.02] border border-[var(--color-border-dim)] rounded-lg p-4">
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4 text-sm text-[var(--color-text-secondary)]">
            <div>
              <h4 className="font-medium text-[var(--color-text-primary)] mb-1">Engine (Python 3.9 + 3.12)</h4>
              <ul className="list-disc list-inside space-y-1">
                <li>qLDPC: BBCode, QCCode, HGPCode, LPCode</li>
                <li>stim + pymatching: Monte Carlo distance estimation</li>
                <li>Evolutionary polynomial search (evolve_polys.py)</li>
                <li>Algebraic analysis (CRT decomposition, root vanishing)</li>
                <li>GF(2) linear algebra, exact distance, REST API</li>
                <li>Nightly automation (launchd, 2 AM, ntfy notifications)</li>
              </ul>
            </div>
            <div>
              <h4 className="font-medium text-[var(--color-text-primary)] mb-1">Dashboard (React)</h4>
              <ul className="list-disc list-inside space-y-1">
                <li>React 19 + TypeScript + Vite</li>
                <li>Zustand for state management</li>
                <li>Plotly.js for scientific visualization</li>
                <li>Tailwind CSS for styling</li>
              </ul>
            </div>
          </div>
        </div>
      </section>
    </div>
  );
}
