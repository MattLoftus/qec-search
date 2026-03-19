"""
Concatenation search for quantum error-correcting codes.

Given inner code [[n1,k1,d1]] and outer code [[n2,k2,d2]]:
  - If k_inner = 1: concatenated code is [[n1*n2, k2, d1*d2]]
    (each logical qubit of outer code is encoded into inner code)
  - General case with k_inner = k_outer: [[n1*n2, k1*k2, >= d1*d2]]
    (but requires k_inner | k_outer for clean nesting)

We check d1*d2 against KNOWN_BOUNDS[(n1*n2, k_concat)] to find
improvements over known lower bounds.
"""

from known_codes import KNOWN_BOUNDS


# -----------------------------------------------------------------------
# Library of codes with known parameters: (name, n, k, d)
# -----------------------------------------------------------------------

CODE_LIBRARY = [
    # Well-known small codes
    ("[[4,2,2]] smallest",          4,  2, 2),
    ("[[5,1,3]] perfect",           5,  1, 3),
    ("[[7,1,3]] Steane",            7,  1, 3),
    ("[[9,1,3]] Shor",              9,  1, 3),
    ("[[15,7,3]] Hamming CSS",     15,  7, 3),

    # From our search results
    ("[[13,1,5]] QT GF(4)",       13,  1, 5),
    ("[[15,3,5]] QT GF(4)",       15,  3, 5),
    ("[[10,1,4]] non-CSS",        10,  1, 4),
    ("[[12,2,4]] non-CSS",        12,  2, 4),

    # Additional well-known codes for completeness
    ("[[6,4,2]] hexacode",         6,  4, 2),
    ("[[8,3,3]]",                  8,  3, 3),
    ("[[11,1,5]]",                11,  1, 5),
    ("[[17,1,7]]",                17,  1, 7),
    ("[[15,1,5]]",                15,  1, 5),
]


def concatenation_params(inner, outer):
    """Compute concatenated code parameters.

    inner = (n1, k1, d1), outer = (n2, k2, d2)

    Case 1: k_inner = 1
      Each logical qubit of outer is replaced by inner block.
      Result: [[n1 * n2, k2, d1 * d2]]

    Case 2: k_inner > 1, k_inner divides k_outer (or k_inner = k_outer)
      More nuanced — we use the standard result for k_inner = k_outer:
      Result: [[n1 * (n2/k_inner * k_inner ... ]]
      Actually for general concatenation with k_inner = k_outer:
        outer has n2 qubits, each encoded into inner → [[n1*n2, k2, >= d1*d2]]
        But this only works cleanly when k_inner = 1, or with matching encoding.

    For simplicity and correctness, we handle:
      - k_inner = 1: standard concatenation → [[n1*n2, k2, d1*d2]]
      - k_inner > 1, k_inner = k_outer: "level-matched" → [[n1*n2, k1*k2, d1*d2]]
        (each qubit of outer replaced by inner; outer's n2 physical qubits become
         n2 inner blocks, but inner encodes k1 logical qubits per block,
         so total logical = ... actually this gets complicated)

    Conservative approach: only claim d1*d2 for k_inner=1 concatenation.
    For k_inner > 1, we still try but flag them separately.
    """
    n1, k1, d1 = inner
    n2, k2, d2 = outer

    if k1 == 1:
        # Standard: each physical qubit of outer → one inner block
        return (n1 * n2, k2, d1 * d2)
    elif k1 == k2:
        # Level-matched: outer's k2 logical qubits encoded into inner
        # This is a recursive encoding. Actually the standard concatenated
        # code with k_inner = k_outer just gives [[n1*n2, k1, d1*d2]]
        # because we're encoding k1 logical qubits at the inner level
        # and those same k1 are the output of the outer encoding.
        # Wait — that's wrong. Let me think more carefully.
        #
        # Standard 2-level concatenation:
        #   Outer [[n2, k2, d2]] maps k2 logical → n2 physical
        #   Inner [[n1, k1, d1]] maps k1 logical → n1 physical
        #   If k1 = k2: we can feed outer's n2 physical qubits through
        #   ceil(n2/k1) inner blocks? No...
        #
        # The clean version: outer has n2 physical qubits.
        # If k_inner = 1: replace each physical qubit with an inner code block.
        # If k_inner > 1: replace each GROUP of k_inner physical qubits with
        #   one inner code block of n1 physical qubits.
        #   This requires k_inner | n2.
        #   Result: [[n1 * (n2 / k1), k2, d1 * d2]]
        if n2 % k1 == 0:
            return (n1 * (n2 // k1), k2, d1 * d2)
        else:
            return None
    else:
        # General case: k_inner divides n_outer's physical qubits
        if k1 > 0 and n2 % k1 == 0:
            return (n1 * (n2 // k1), k2, d1 * d2)
        else:
            return None


def run_search():
    """Run concatenation search over all code pairs."""

    print("=" * 72)
    print("QUANTUM CODE CONCATENATION SEARCH")
    print("=" * 72)
    print()
    print(f"Code library: {len(CODE_LIBRARY)} codes")
    print(f"Known bounds: {len(KNOWN_BOUNDS)} entries (up to n=60)")
    print()

    # Separate by k value for easy lookup
    k1_codes = [(name, n, k, d) for name, n, k, d in CODE_LIBRARY if k == 1]
    all_codes = CODE_LIBRARY

    results = []

    print("-" * 72)
    print("SECTION 1: k_inner = 1 concatenation (standard)")
    print("  inner [[n1,1,d1]] (x) outer [[n2,k2,d2]] -> [[n1*n2, k2, d1*d2]]")
    print("-" * 72)
    print()

    for iname, in_, ik, id_ in k1_codes:
        for oname, on, ok, od in all_codes:
            # Skip self-concatenation with trivial codes
            n_cat = in_ * on
            k_cat = ok
            d_cat = id_ * od

            # Skip if n too large for our bounds table
            if n_cat > 200:
                continue

            # Look up known bounds
            key = (n_cat, k_cat)
            if key in KNOWN_BOUNDS:
                d_lower, d_upper = KNOWN_BOUNDS[key]

                status = None
                if d_cat > d_upper:
                    # Exceeds theoretical upper bound — impossible, something wrong
                    status = "IMPOSSIBLE (exceeds upper bound)"
                elif d_cat > d_lower:
                    status = "IMPROVEMENT"
                elif d_cat == d_lower:
                    status = "MATCHES"
                elif d_cat >= d_lower - 1:
                    status = "NEAR-MISS"
                else:
                    status = None  # not interesting

                if status:
                    results.append((status, n_cat, k_cat, d_cat, d_lower, d_upper,
                                    iname, oname, "k1-standard"))
                    marker = "***" if status == "IMPROVEMENT" else "   "
                    print(f"  {marker} {iname} (x) {oname}")
                    print(f"      -> [[{n_cat}, {k_cat}, {d_cat}]]  "
                          f"bounds: d in [{d_lower}, {d_upper}]  => {status}")
            else:
                # No known bound — worth reporting if parameters are reasonable
                if n_cat <= 100 and k_cat >= 1:
                    results.append(("NO BOUND", n_cat, k_cat, d_cat, None, None,
                                    iname, oname, "k1-standard"))
                    print(f"      {iname} (x) {oname}")
                    print(f"      -> [[{n_cat}, {k_cat}, {d_cat}]]  "
                          f"(no known bound for (n={n_cat}, k={k_cat}))")

    print()
    print("-" * 72)
    print("SECTION 2: k_inner > 1 concatenation (group encoding)")
    print("  inner [[n1,k1,d1]] (x) outer [[n2,k2,d2]] where k1 | n2")
    print("  -> [[n1 * (n2/k1), k2, d1*d2]]")
    print("-" * 72)
    print()

    ki_gt1_codes = [(name, n, k, d) for name, n, k, d in CODE_LIBRARY if k > 1]

    for iname, in_, ik, id_ in ki_gt1_codes:
        for oname, on, ok, od in all_codes:
            if on % ik != 0:
                continue

            n_cat = in_ * (on // ik)
            k_cat = ok
            d_cat = id_ * od

            if n_cat > 200 or n_cat < 1:
                continue

            key = (n_cat, k_cat)
            if key in KNOWN_BOUNDS:
                d_lower, d_upper = KNOWN_BOUNDS[key]

                status = None
                if d_cat > d_upper:
                    status = "IMPOSSIBLE (exceeds upper bound)"
                elif d_cat > d_lower:
                    status = "IMPROVEMENT"
                elif d_cat == d_lower:
                    status = "MATCHES"
                elif d_cat >= d_lower - 1:
                    status = "NEAR-MISS"

                if status:
                    results.append((status, n_cat, k_cat, d_cat, d_lower, d_upper,
                                    iname, oname, "ki>1-group"))
                    marker = "***" if status == "IMPROVEMENT" else "   "
                    print(f"  {marker} {iname} (x) {oname}")
                    print(f"      -> [[{n_cat}, {k_cat}, {d_cat}]]  "
                          f"bounds: d in [{d_lower}, {d_upper}]  => {status}")
            else:
                if n_cat <= 100 and k_cat >= 1:
                    results.append(("NO BOUND", n_cat, k_cat, d_cat, None, None,
                                    iname, oname, "ki>1-group"))
                    print(f"      {iname} (x) {oname}")
                    print(f"      -> [[{n_cat}, {k_cat}, {d_cat}]]  "
                          f"(no known bound for (n={n_cat}, k={k_cat}))")

    # Summary
    print()
    print("=" * 72)
    print("SUMMARY")
    print("=" * 72)

    improvements = [r for r in results if r[0] == "IMPROVEMENT"]
    matches = [r for r in results if r[0] == "MATCHES"]
    near_misses = [r for r in results if r[0] == "NEAR-MISS"]
    impossible = [r for r in results if "IMPOSSIBLE" in r[0]]
    no_bound = [r for r in results if r[0] == "NO BOUND"]

    print(f"\n  IMPROVEMENTS over known lower bound: {len(improvements)}")
    for r in improvements:
        status, nc, kc, dc, dl, du, iname, oname, method = r
        print(f"    [[{nc},{kc},{dc}]] from {iname} (x) {oname}")
        print(f"      d_concat={dc} > d_lower={dl}  (upper={du}, method={method})")

    print(f"\n  MATCHES known lower bound: {len(matches)}")
    for r in matches:
        status, nc, kc, dc, dl, du, iname, oname, method = r
        print(f"    [[{nc},{kc},{dc}]] from {iname} (x) {oname}  [{method}]")

    print(f"\n  NEAR-MISSES (d_concat = d_lower - 1): {len(near_misses)}")
    for r in near_misses:
        status, nc, kc, dc, dl, du, iname, oname, method = r
        print(f"    [[{nc},{kc},{dc}]] from {iname} (x) {oname}  "
              f"(need {dl}, got {dc}) [{method}]")

    if impossible:
        print(f"\n  IMPOSSIBLE (exceeds upper bound — check logic): {len(impossible)}")
        for r in impossible:
            status, nc, kc, dc, dl, du, iname, oname, method = r
            print(f"    [[{nc},{kc},{dc}]] from {iname} (x) {oname}  "
                  f"(upper={du}) [{method}]")

    print(f"\n  NO KNOWN BOUND (n <= 100): {len(no_bound)}")
    for r in no_bound:
        status, nc, kc, dc, dl, du, iname, oname, method = r
        print(f"    [[{nc},{kc},{dc}]] from {iname} (x) {oname}  [{method}]")

    print()
    print("=" * 72)
    if improvements:
        print("RESULT: Found potential improvements via concatenation!")
        print("NOTE: These need independent verification — concatenation gives")
        print("a guaranteed LOWER bound on distance (d >= d1*d2), but the actual")
        print("code must be constructed to confirm.")
    else:
        print("RESULT: No improvements found over known bounds.")
        print("This is expected — concatenation is well-studied, and modern")
        print("algebraic codes generally outperform concatenated constructions")
        print("at small block lengths.")
    print("=" * 72)


if __name__ == "__main__":
    run_search()
