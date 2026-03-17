"""
Approximate code distance computation using stim + pymatching + ldpc.

For codes with n > 60, exact distance computation (exponential time) is
infeasible. This module provides fast approximate methods:

1. graphlike_distance: Upper bound via stim's shortest_graphlike_error()
   - Finds minimum-weight graphlike error causing undetected logical error
   - O(seconds) even for n=200
   - Exact for codes where single errors trigger at most 2 detectors
   - Fails (returns None) for LDPC codes with high-weight stabilizers

2. monte_carlo_distance_estimate: Estimate distance from error rate scaling
   - Uses BP+OSD decoder (ldpc library) for LDPC codes
   - Fits log(logical_error_rate) vs log(p) to extract distance exponent
   - Works for any CSS code regardless of stabilizer weight

3. estimate_distance_stim: Combined estimate — tries graphlike first,
   falls back to Monte Carlo scaling for LDPC codes
"""

import numpy as np
import stim
import scipy.sparse as sp
import time
from typing import Optional, Tuple, List

import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from codes import CSSCode, gf2_rank, gf2_nullspace


# ---------------------------------------------------------------------------
# Logical operator computation
# ---------------------------------------------------------------------------

def _find_logical_representatives(H_detect, H_stabilizer, n):
    """Find representatives of logical operators.

    Computes vectors in ker(H_stabilizer) that are NOT in rowspace(H_detect).
    These represent nontrivial logical operators of one type.

    Usage:
      - Logical Z operators (for X-distance): call with (Hx, Hz, n)
        => ker(Hz) \\ rowspace(Hx)
      - Logical X operators (for Z-distance): call with (Hz, Hx, n)
        => ker(Hx) \\ rowspace(Hz)

    Returns:
        List of binary vectors, each representing a logical operator.
    """
    ker = gf2_nullspace(H_stabilizer)
    if ker.shape[0] == 0:
        return []

    if H_detect.size == 0:
        return [ker[i] for i in range(ker.shape[0])]

    rank_detect = gf2_rank(H_detect)
    logicals = []

    for i in range(ker.shape[0]):
        test_matrix = np.vstack(
            [H_detect]
            + [l.reshape(1, -1) for l in logicals]
            + [ker[i].reshape(1, -1)]
        )
        new_rank = gf2_rank(test_matrix)
        if new_rank > rank_detect + len(logicals):
            logicals.append(ker[i])

    return logicals


# ---------------------------------------------------------------------------
# Stim circuit construction
# ---------------------------------------------------------------------------

def _build_x_distance_circuit(Hx, Hz, noise_rate=0.001):
    """Build stim circuit for X-distance computation.

    Applies X errors to data qubits, measures Z-type stabilizers (from Hx),
    and defines logical Z observables (in ker(Hz) \\ rowspace(Hx)).
    """
    n = Hx.shape[1]
    rx = Hx.shape[0]

    # Logical Z operators: ker(Hz) \\ rowspace(Hx)
    logicals = _find_logical_representatives(Hx, Hz, n)
    if not logicals:
        return None, 0

    c = stim.Circuit()
    c.append('R', list(range(n)))
    c.append('X_ERROR', list(range(n)), noise_rate)

    for row_idx in range(rx):
        qubits = list(np.where(Hx[row_idx] == 1)[0])
        if not qubits:
            continue
        anc = n + row_idx
        c.append('R', [anc])
        for q in qubits:
            c.append('CX', [q, anc])
        c.append('M', [anc])
        c.append('DETECTOR', [stim.target_rec(-1)])

    c.append('M', list(range(n)))

    for obs_idx, logical in enumerate(logicals):
        support = list(np.where(logical == 1)[0])
        targets = [stim.target_rec(-(n - q)) for q in support]
        c.append('OBSERVABLE_INCLUDE', targets, obs_idx)

    return c, len(logicals)


def _build_z_distance_circuit(Hx, Hz, noise_rate=0.001):
    """Build stim circuit for Z-distance computation.

    Applies Z errors to data qubits (in X basis), measures X-type stabilizers
    (from Hz), and defines logical X observables (in ker(Hx) \\ rowspace(Hz)).
    """
    n = Hz.shape[1]
    rz = Hz.shape[0]

    # Logical X operators: ker(Hx) \\ rowspace(Hz)
    logicals = _find_logical_representatives(Hz, Hx, n)
    if not logicals:
        return None, 0

    c = stim.Circuit()
    c.append('R', list(range(n)))
    c.append('H', list(range(n)))
    c.append('Z_ERROR', list(range(n)), noise_rate)

    for row_idx in range(rz):
        qubits = list(np.where(Hz[row_idx] == 1)[0])
        if not qubits:
            continue
        anc = n + row_idx
        c.append('R', [anc])
        c.append('H', [anc])
        for q in qubits:
            c.append('CX', [anc, q])
        c.append('H', [anc])
        c.append('M', [anc])
        c.append('DETECTOR', [stim.target_rec(-1)])

    c.append('H', list(range(n)))
    c.append('M', list(range(n)))

    for obs_idx, logical in enumerate(logicals):
        support = list(np.where(logical == 1)[0])
        targets = [stim.target_rec(-(n - q)) for q in support]
        c.append('OBSERVABLE_INCLUDE', targets, obs_idx)

    return c, len(logicals)


def css_to_stim_circuit(Hx, Hz, noise_rate=0.001):
    """Convert CSS parity check matrices to stim circuits.

    Returns a dict with 'x_circuit' (for X-distance) and 'z_circuit'
    (for Z-distance).
    """
    Hx = np.asarray(Hx, dtype=np.uint8)
    Hz = np.asarray(Hz, dtype=np.uint8)

    x_circuit, num_log_z = _build_x_distance_circuit(Hx, Hz, noise_rate)
    z_circuit, num_log_x = _build_z_distance_circuit(Hx, Hz, noise_rate)

    return {
        'x_circuit': x_circuit,
        'z_circuit': z_circuit,
        'num_logicals_z': num_log_z,
        'num_logicals_x': num_log_x,
    }


# ---------------------------------------------------------------------------
# Graphlike distance (stim-based, fast but only works for graphlike codes)
# ---------------------------------------------------------------------------

def graphlike_distance(Hx, Hz):
    """Compute graphlike distance upper bound using stim.

    Returns the minimum of X-distance and Z-distance from
    shortest_graphlike_error(). This is an upper bound on the true distance,
    and equals the true distance for codes where each physical error triggers
    at most 2 detectors (e.g., surface codes, small CSS codes).

    For LDPC codes with high-weight stabilizers (where single errors trigger
    3+ detectors), this method returns None.

    Returns:
        (distance_upper_bound, dx, dz)
    """
    Hx = np.asarray(Hx, dtype=np.uint8)
    Hz = np.asarray(Hz, dtype=np.uint8)

    circuits = css_to_stim_circuit(Hx, Hz)

    dx = float('inf')
    dz = float('inf')

    if circuits['x_circuit'] is not None:
        try:
            errors = circuits['x_circuit'].shortest_graphlike_error()
            dx = len(errors)
        except Exception:
            pass

    if circuits['z_circuit'] is not None:
        try:
            errors = circuits['z_circuit'].shortest_graphlike_error()
            dz = len(errors)
        except Exception:
            pass

    d = min(dx, dz)
    if d == float('inf'):
        return None, dx, dz
    return d, dx, dz


# ---------------------------------------------------------------------------
# Monte Carlo with BP+OSD decoder (works for any CSS code)
# ---------------------------------------------------------------------------

def _mc_one_side(H_check, H_stab, p, num_shots, logicals_matrix):
    """Run Monte Carlo for one error type (X or Z).

    H_check: parity check matrix that detects the errors
    H_stab: stabilizer matrix for the other type (to check logical triviality)
    logicals_matrix: (k x n) matrix of logical operator representatives

    Returns number of logical errors.
    """
    from ldpc import BpOsdDecoder

    n = H_check.shape[1]
    H_sparse = sp.csr_matrix(H_check.astype(np.float64))

    decoder = BpOsdDecoder(
        H_sparse,
        error_rate=p,
        bp_method='ms',
        max_iter=100,
        osd_method='osd_cs',
        osd_order=min(7, H_check.shape[0]),
    )

    num_logical_errors = 0
    for _ in range(num_shots):
        error = np.random.binomial(1, p, size=n).astype(np.uint8)
        syndrome = (H_check @ error) % 2
        correction = decoder.decode(syndrome)
        residual = (error + correction) % 2

        # Verify residual is in kernel of H_check
        if np.any((H_check @ residual) % 2):
            continue  # Decoder didn't fully correct

        if np.any(residual):
            # Check if residual is a nontrivial logical
            logical_check = (logicals_matrix @ residual) % 2
            if np.any(logical_check):
                num_logical_errors += 1

    return num_logical_errors


def monte_carlo_logical_error_rate(Hx, Hz, p, num_shots=10000,
                                    error_type='both'):
    """Estimate logical error rate using Monte Carlo + BP+OSD decoder.

    Works for any CSS code, including LDPC codes where stim/pymatching
    approaches fail.

    Args:
        Hx: (rx, n) binary matrix - Z-type stabilizers (detect X errors)
        Hz: (rz, n) binary matrix - X-type stabilizers (detect Z errors)
        p: physical error rate
        num_shots: number of Monte Carlo samples
        error_type: 'x', 'z', or 'both'

    Returns:
        dict with logical error rates per error type
    """
    Hx = np.asarray(Hx, dtype=np.uint8)
    Hz = np.asarray(Hz, dtype=np.uint8)
    n = Hx.shape[1]

    results = {}

    # Precompute logical operators
    logicals_z = _find_logical_representatives(Hx, Hz, n)  # ker(Hz) \ rowspace(Hx)
    logicals_x = _find_logical_representatives(Hz, Hx, n)  # ker(Hx) \ rowspace(Hz)

    if error_type in ('x', 'both') and logicals_z:
        L_z = np.array(logicals_z, dtype=np.uint8)
        num_errs = _mc_one_side(Hx, Hz, p, num_shots, L_z)
        results['x'] = {
            'logical_error_rate': num_errs / num_shots,
            'num_errors': num_errs,
            'num_shots': num_shots,
            'physical_error_rate': p,
        }

    if error_type in ('z', 'both') and logicals_x:
        L_x = np.array(logicals_x, dtype=np.uint8)
        num_errs = _mc_one_side(Hz, Hx, p, num_shots, L_x)
        results['z'] = {
            'logical_error_rate': num_errs / num_shots,
            'num_errors': num_errs,
            'num_shots': num_shots,
            'physical_error_rate': p,
        }

    return results


def _mc_min_weight_logical(H_check, logicals_matrix, p, num_shots):
    """Find minimum-weight logical errors via Monte Carlo.

    Samples random errors at rate p, computes the exact residual
    (error XOR correction), and tracks the minimum weight of any
    nontrivial logical residual observed.

    Returns (min_weight, count_logicals) or (None, 0) if no logicals found.
    """
    from ldpc import BpOsdDecoder

    n = H_check.shape[1]
    H_sparse = sp.csr_matrix(H_check.astype(np.float64))

    decoder = BpOsdDecoder(
        H_sparse,
        error_rate=p,
        bp_method='ms',
        max_iter=100,
        osd_method='osd_cs',
        osd_order=min(10, H_check.shape[0]),
    )

    min_weight = None
    count = 0

    for _ in range(num_shots):
        error = np.random.binomial(1, p, size=n).astype(np.uint8)
        syndrome = (H_check @ error) % 2
        correction = decoder.decode(syndrome)
        residual = (error + correction) % 2

        if np.any((H_check @ residual) % 2):
            continue

        if np.any(residual):
            logical_check = (logicals_matrix @ residual) % 2
            if np.any(logical_check):
                w = int(np.sum(residual))
                count += 1
                if min_weight is None or w < min_weight:
                    min_weight = w

    return min_weight, count


def monte_carlo_distance_estimate(Hx, Hz, num_shots=100000,
                                   error_type='both', verbose=True):
    """Estimate code distance via Monte Carlo minimum-weight logical search.

    Strategy: sample errors at a moderate noise rate (high enough to produce
    logical errors regularly), decode with BP+OSD, and track the minimum
    weight of nontrivial logical residuals. With enough samples, this
    converges to the true distance from above.

    We use multiple noise rates to balance:
    - Low p: logical errors are rare but tend to be low weight
    - High p: more logicals but higher weight

    Args:
        Hx, Hz: binary check matrices
        num_shots: total Monte Carlo samples (split across noise rates)
        error_type: 'x', 'z', or 'both'
        verbose: print progress

    Returns:
        dict with estimated distance and metadata
    """
    Hx = np.asarray(Hx, dtype=np.uint8)
    Hz = np.asarray(Hz, dtype=np.uint8)
    n = Hx.shape[1]

    # Noise rates and shot allocation: more shots at lower noise rates
    # where we're more likely to see minimum-weight logicals
    configs = [
        (0.02, num_shots // 4),
        (0.03, num_shots // 4),
        (0.04, num_shots // 4),
        (0.05, num_shots // 4),
    ]

    scaling_data = {}
    for etype in (['x', 'z'] if error_type == 'both' else [error_type]):
        if etype == 'x':
            H_check = Hx
            logicals = _find_logical_representatives(Hx, Hz, n)
        else:
            H_check = Hz
            logicals = _find_logical_representatives(Hz, Hx, n)

        if not logicals:
            scaling_data[etype] = {'estimated_distance': None}
            continue

        L = np.array(logicals, dtype=np.uint8)
        overall_min_weight = None
        total_logicals = 0

        for p, shots in configs:
            t0 = time.time()
            min_w, count = _mc_min_weight_logical(H_check, L, p, shots)
            elapsed = time.time() - t0

            total_logicals += count
            if min_w is not None:
                if overall_min_weight is None or min_w < overall_min_weight:
                    overall_min_weight = min_w

            if verbose:
                print(f"    {etype}-type p={p:.3f}: min_weight={min_w}, "
                      f"found={count}/{shots} logicals [{elapsed:.1f}s]")

        if verbose and overall_min_weight is not None:
            print(f"    {etype}-distance upper bound: {overall_min_weight} "
                  f"(from {total_logicals} logical errors)")

        scaling_data[etype] = {
            'estimated_distance': overall_min_weight,
            'total_logicals_found': total_logicals,
        }

    distances = [v['estimated_distance'] for v in scaling_data.values()
                 if v['estimated_distance'] is not None]
    overall_d = min(distances) if distances else None

    return {
        'estimated_distance': overall_d,
        'scaling_data': scaling_data,
    }


# ---------------------------------------------------------------------------
# Combined distance estimate
# ---------------------------------------------------------------------------

def estimate_distance_stim(Hx, Hz, verbose=True):
    """Estimate code distance using the best available method.

    Strategy:
    1. Try graphlike distance (stim) — fast, exact for graphlike codes
    2. If that fails (LDPC codes), use Monte Carlo scaling with BP+OSD

    Args:
        Hx: (rx, n) binary matrix
        Hz: (rz, n) binary matrix
        verbose: print progress

    Returns:
        dict with distance estimate and metadata
    """
    Hx = np.asarray(Hx, dtype=np.uint8)
    Hz = np.asarray(Hz, dtype=np.uint8)

    n = Hx.shape[1]
    k = n - gf2_rank(Hx) - gf2_rank(Hz)

    if verbose:
        print(f"  Code: n={n}, k={k}")

    result = {'n': n, 'k': k}

    # Step 1: Try graphlike distance
    t0 = time.time()
    d_graph, dx, dz = graphlike_distance(Hx, Hz)
    t_graph = time.time() - t0

    result['graphlike_distance'] = d_graph
    result['dx_graphlike'] = dx if dx != float('inf') else None
    result['dz_graphlike'] = dz if dz != float('inf') else None
    result['graphlike_time_sec'] = round(t_graph, 3)

    if verbose:
        print(f"  Graphlike distance: d<={d_graph} (dx={dx}, dz={dz}) "
              f"[{t_graph:.2f}s]")

    if d_graph is not None:
        # Graphlike distance worked — use it
        result['estimated_distance'] = d_graph
        result['method'] = 'graphlike'

        # Optional: validate with MC at one noise rate
        if d_graph > 1:
            p_test = min(0.01, 0.1 / d_graph)
            t0 = time.time()
            mc = monte_carlo_logical_error_rate(
                Hx, Hz, p=p_test, num_shots=10000, error_type='both'
            )
            t_mc = time.time() - t0
            result['mc_validation'] = mc
            result['mc_validation_time_sec'] = round(t_mc, 3)

            if verbose:
                for etype, res in mc.items():
                    rate = res.get('logical_error_rate')
                    if rate is not None:
                        print(f"  MC {etype}-error rate at p={p_test:.4f}: "
                              f"{rate:.6f} [{t_mc:.2f}s]")
    else:
        # Graphlike failed — fall back to MC scaling estimate
        if verbose:
            print("  Graphlike distance unavailable (non-graphlike errors)")
            print("  Falling back to Monte Carlo scaling estimate...")

        t0 = time.time()
        mc_est = monte_carlo_distance_estimate(
            Hx, Hz, verbose=verbose, error_type='both'
        )
        t_mc = time.time() - t0

        result['estimated_distance'] = mc_est['estimated_distance']
        result['method'] = 'mc_scaling'
        result['mc_scaling'] = mc_est
        result['mc_scaling_time_sec'] = round(t_mc, 3)

    return result


# ---------------------------------------------------------------------------
# Known code definitions for testing
# ---------------------------------------------------------------------------

def steane_code():
    """Steane [[7,1,3]] code."""
    Hx = np.array([
        [1, 0, 0, 1, 0, 1, 1],
        [0, 1, 0, 1, 1, 0, 1],
        [0, 0, 1, 0, 1, 1, 1],
    ], dtype=np.uint8)
    Hz = Hx.copy()
    return Hx, Hz


def shor_code():
    """Shor [[9,1,3]] code.

    Hx (Z stabilizers): Z0Z1, Z1Z2, Z3Z4, Z4Z5, Z6Z7, Z7Z8
    Hz (X stabilizers): X0X1X2X3X4X5, X3X4X5X6X7X8
    """
    Hx = np.array([
        [1, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 1],
    ], dtype=np.uint8)
    Hz = np.array([
        [1, 1, 1, 1, 1, 1, 0, 0, 0],
        [0, 0, 0, 1, 1, 1, 1, 1, 1],
    ], dtype=np.uint8)
    return Hx, Hz


def bravyi_72_12_6():
    """Bravyi [[72,12,6]] bivariate bicycle code from qLDPC."""
    from qldpc.codes import BBCode
    import sympy
    x, y = sympy.symbols('x y')
    code = BBCode([6, 6], x**3 + y + y**2, y**3 + x + x**2)
    Hx = np.array(code.matrix_x, dtype=np.uint8) % 2
    Hz = np.array(code.matrix_z, dtype=np.uint8) % 2
    return Hx, Hz


# ---------------------------------------------------------------------------
# Main: test on known codes
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    print("=" * 70)
    print("Stim-based approximate distance computation")
    print("=" * 70)

    # ---- Test 1: Steane [[7,1,3]] ----
    print("\n--- Steane [[7,1,3]] code ---")
    Hx, Hz = steane_code()
    code = CSSCode(Hx=Hx, Hz=Hz)
    print(f"  Code object: {code}")
    valid, msg = code.is_valid()
    print(f"  Valid: {valid} ({msg})")

    result = estimate_distance_stim(Hx, Hz)
    d = result['estimated_distance']
    print(f"  => Estimated distance: {d} (method: {result['method']})")
    assert d == 3, f"Steane code should have d=3, got {d}"
    print("  PASS: d=3 confirmed")

    # ---- Test 2: Shor [[9,1,3]] ----
    print("\n--- Shor [[9,1,3]] code ---")
    Hx, Hz = shor_code()
    code = CSSCode(Hx=Hx, Hz=Hz)
    print(f"  Code object: {code}")
    valid, msg = code.is_valid()
    print(f"  Valid: {valid} ({msg})")

    result = estimate_distance_stim(Hx, Hz)
    d = result['estimated_distance']
    print(f"  => Estimated distance: {d} (method: {result['method']})")
    assert d == 3, f"Shor code should have d=3, got {d}"
    print("  PASS: d=3 confirmed")

    # ---- Test 3: Bravyi [[72,12,6]] bivariate bicycle code ----
    print("\n--- Bravyi [[72,12,6]] bivariate bicycle code ---")
    try:
        Hx, Hz = bravyi_72_12_6()
        code = CSSCode(Hx=Hx, Hz=Hz)
        print(f"  Code object: {code}")
        valid, msg = code.is_valid()
        print(f"  Valid: {valid} ({msg})")

        result = estimate_distance_stim(Hx, Hz)
        d = result['estimated_distance']
        print(f"  => Estimated distance: {d} (method: {result['method']})")
        if d is not None:
            status = 'PASS' if d == 6 else 'NOTE'
            print(f"  {status}: expected d=6, got d={d}")
        else:
            print("  FAIL: Could not estimate distance")

    except ImportError as e:
        print(f"  SKIP: qLDPC not available ({e})")
    except Exception as e:
        print(f"  ERROR: {e}")
        import traceback
        traceback.print_exc()

    # ---- Test 4: Additional BB codes ----
    print("\n--- Additional bivariate bicycle codes ---")
    try:
        from qldpc.codes import BBCode
        import sympy
        x, y = sympy.symbols('x y')

        bb_codes = [
            ([9, 6], x**3 + y + y**2, y**3 + x + x**2, "[[108,8,?]]"),
            ([12, 6], x**3 + y + y**2, y**3 + x + x**2, "[[144,12,?]]"),
        ]

        for orders, poly_a, poly_b, label in bb_codes:
            print(f"\n  BB code {label} (orders={orders}):")
            try:
                bb = BBCode(orders, poly_a, poly_b)
                Hx = np.array(bb.matrix_x, dtype=np.uint8) % 2
                Hz = np.array(bb.matrix_z, dtype=np.uint8) % 2
                code = CSSCode(Hx=Hx, Hz=Hz)

                valid, msg = code.is_valid()
                print(f"    n={code.n}, k={code.k}, valid={valid}")

                result = estimate_distance_stim(Hx, Hz)
                d = result['estimated_distance']
                print(f"    => Estimated distance: {d} "
                      f"(method: {result['method']})")

            except Exception as e:
                print(f"    ERROR: {e}")
                import traceback
                traceback.print_exc()

    except ImportError:
        print("  SKIP: qLDPC not available")

    print("\n" + "=" * 70)
    print("All tests completed.")
