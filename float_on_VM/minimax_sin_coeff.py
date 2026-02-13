#!/usr/bin/env python3
"""Compute minimax polynomial coefficients for sin(x) using Remez exchange.

Model:
  sin(x) ~= c0*x + c1*x^3 + c2*x^5 + ... + cn*x^(2n+1)

Default interval:
  x in [-pi/2, pi/2]
"""

from __future__ import annotations

import argparse
import math
from typing import List, Sequence, Tuple


def solve_linear_system(a: List[List[float]], b: List[float]) -> List[float]:
    """Gaussian elimination with partial pivoting."""
    n = len(a)
    for i in range(n):
        # Pivot
        piv = i
        max_abs = abs(a[i][i])
        for r in range(i + 1, n):
            v = abs(a[r][i])
            if v > max_abs:
                max_abs = v
                piv = r
        if max_abs == 0.0:
            raise ValueError("Singular matrix while solving Remez system")
        if piv != i:
            a[i], a[piv] = a[piv], a[i]
            b[i], b[piv] = b[piv], b[i]

        # Normalize row i
        diag = a[i][i]
        for c in range(i, n):
            a[i][c] /= diag
        b[i] /= diag

        # Eliminate
        for r in range(n):
            if r == i:
                continue
            f = a[r][i]
            if f == 0.0:
                continue
            for c in range(i, n):
                a[r][c] -= f * a[i][c]
            b[r] -= f * b[i]

    return b


def poly_odd_eval(x: float, coeffs: Sequence[float]) -> float:
    """Evaluate c0*x + c1*x^3 + ... using Horner on x^2."""
    z = x * x
    p = coeffs[-1]
    for c in reversed(coeffs[:-1]):
        p = p * z + c
    return x * p


def make_initial_points(a: float, m: int) -> List[float]:
    """Chebyshev-like initial reference points on [0, a]."""
    pts = []
    for i in range(m):
        # m points including 0 and a
        t = i / (m - 1)
        x = 0.5 * a * (1.0 - math.cos(math.pi * t))
        pts.append(x)
    return pts


def local_extrema_candidates(xs: Sequence[float], es: Sequence[float]) -> List[Tuple[float, float]]:
    out: List[Tuple[float, float]] = []
    n = len(xs)
    # Include endpoints
    out.append((xs[0], es[0]))
    for i in range(1, n - 1):
        e0 = es[i - 1]
        e1 = es[i]
        e2 = es[i + 1]
        if (e1 >= e0 and e1 >= e2) or (e1 <= e0 and e1 <= e2):
            out.append((xs[i], e1))
    out.append((xs[-1], es[-1]))
    return out


def alternating_subsequence(cands: Sequence[Tuple[float, float]], m: int) -> List[Tuple[float, float]]:
    """Build an alternating-sign list, then choose length-m contiguous best block."""
    if not cands:
        return []

    # Greedy collapse to alternating signs, keep larger |e| for same-sign runs.
    alt: List[Tuple[float, float]] = []
    for x, e in cands:
        if e == 0.0:
            continue
        s = 1.0 if e > 0.0 else -1.0
        if not alt:
            alt.append((x, e))
            continue
        last_s = 1.0 if alt[-1][1] > 0.0 else -1.0
        if s == last_s:
            if abs(e) > abs(alt[-1][1]):
                alt[-1] = (x, e)
        else:
            alt.append((x, e))

    if len(alt) < m:
        # Fallback: take m largest |e| and sort by x (may lose strict alternation).
        top = sorted(cands, key=lambda t: abs(t[1]), reverse=True)[:m]
        return sorted(top, key=lambda t: t[0])

    if len(alt) == m:
        return alt

    # Choose the best contiguous window of length m by maximin criterion.
    best_i = 0
    best_score = -1.0
    for i in range(0, len(alt) - m + 1):
        window = alt[i : i + m]
        score = min(abs(e) for _, e in window)
        if score > best_score:
            best_score = score
            best_i = i
    return alt[best_i : best_i + m]


def solve_remez_step(points: Sequence[float], degree_n: int) -> Tuple[List[float], float]:
    """Solve coefficients and equioscillation error for current reference points."""
    # Unknowns: c0..cn and E => n+2 unknowns
    m = degree_n + 2
    a_mat: List[List[float]] = []
    b_vec: List[float] = []
    for i, x in enumerate(points):
        row = []
        x_pow = x
        for _ in range(degree_n + 1):
            row.append(x_pow)
            x_pow *= x * x
        row.append(1.0 if (i % 2 == 0) else -1.0)  # alternating error sign
        a_mat.append(row)
        b_vec.append(math.sin(x))

    sol = solve_linear_system(a_mat, b_vec)
    coeffs = sol[:-1]
    err = sol[-1]
    return coeffs, err


def max_error_on_grid(coeffs: Sequence[float], a: float, grid_n: int) -> Tuple[float, float]:
    max_abs = -1.0
    at_x = 0.0
    for i in range(grid_n + 1):
        x = a * i / grid_n
        e = poly_odd_eval(x, coeffs) - math.sin(x)
        ae = abs(e)
        if ae > max_abs:
            max_abs = ae
            at_x = x
    return max_abs, at_x


def remez_sin_odd(
    degree_n: int,
    interval_a: float,
    max_iter: int,
    grid_n: int,
    tol: float,
) -> Tuple[List[float], float]:
    """Return odd polynomial coeffs for sin on [-a, a]."""
    m = degree_n + 2
    pts = make_initial_points(interval_a, m)

    last_max = None
    coeffs: List[float] = []
    e_ref = 0.0

    for _ in range(max_iter):
        coeffs, e_ref = solve_remez_step(pts, degree_n)

        # Build dense-grid error and new reference points
        xs = [interval_a * i / grid_n for i in range(grid_n + 1)]
        es = [poly_odd_eval(x, coeffs) - math.sin(x) for x in xs]
        cands = local_extrema_candidates(xs, es)
        new_pts_with_e = alternating_subsequence(cands, m)
        new_pts = [x for x, _ in new_pts_with_e]
        new_pts.sort()
        pts = new_pts

        max_abs = max(abs(e) for e in es)
        if last_max is not None and abs(max_abs - last_max) < tol:
            break
        last_max = max_abs

    return coeffs, e_ref


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute minimax coefficients for sin(x).")
    parser.add_argument("--degree-n", type=int, default=6, help="n in odd polynomial up to x^(2n+1)")
    parser.add_argument(
        "--interval-a",
        type=float,
        default=math.pi / 2.0,
        help="Approximation interval is [-a, a]",
    )
    parser.add_argument("--max-iter", type=int, default=20)
    parser.add_argument("--grid-n", type=int, default=40000)
    parser.add_argument("--tol", type=float, default=1e-15)
    args = parser.parse_args()

    coeffs, e_ref = remez_sin_odd(
        degree_n=args.degree_n,
        interval_a=args.interval_a,
        max_iter=args.max_iter,
        grid_n=args.grid_n,
        tol=args.tol,
    )
    max_err, at_x = max_error_on_grid(coeffs, args.interval_a, args.grid_n)

    print("Model:")
    print("  sin(x) ~= c0*x + c1*x^3 + ... + cn*x^(2n+1)")
    print(f"degree_n = {args.degree_n} (highest power = x^{2*args.degree_n + 1})")
    print(f"interval = [-{args.interval_a}, {args.interval_a}]")
    print(f"reference_error_E ~= {e_ref:+.18e}")
    print(f"max_abs_error_on_grid ~= {max_err:.18e} at x={at_x:.12f}")
    print("")
    print("coefficients:")
    for i, c in enumerate(coeffs):
        p = 2 * i + 1
        print(f"  c{i} (x^{p}): {c:+.18e}")

    print("")
    print("comparison with Taylor coefficients:")
    print("  power | minimax | taylor | diff(%)")
    for i, c in enumerate(coeffs):
        p = 2 * i + 1
        taylor = ((-1.0) ** i) / float(math.factorial(p))
        diff_pct = (c - taylor) / taylor * 100.0
        print(f"  x^{p:<2d} | {c:+.18e} | {taylor:+.18e} | {diff_pct:+.9e}%")


if __name__ == "__main__":
    main()
