#!/usr/bin/env python3
"""CORDIC implementation for basic trigonometric functions.

Usage:
  python cordic.py 0.5
  python cordic.py 30 --deg
"""

from __future__ import annotations

import argparse
import sys
from typing import Dict, Tuple

PI = 3.1415926535897932384626433832795028841971
TWO_PI = 2.0 * PI
HALF_PI = PI / 2.0

# atan(2^-i), i in [0, 31]
CORDIC_ANGLES = [
    0.7853981633974483,
    0.4636476090008061,
    0.24497866312686414,
    0.12435499454676144,
    0.06241880999595735,
    0.031239833430268277,
    0.015623728620476831,
    0.007812341060101111,
    0.0039062301319669718,
    0.0019531225164788188,
    0.0009765621895593195,
    0.0004882812111948983,
    0.00024414062014936177,
    0.00012207031189367021,
    6.103515617420877e-05,
    3.0517578115526096e-05,
    1.5258789061315762e-05,
    7.62939453110197e-06,
    3.814697265606496e-06,
    1.907348632810187e-06,
    9.536743164059608e-07,
    4.7683715820308884e-07,
    2.3841857910155797e-07,
    1.1920928955078068e-07,
    5.960464477539055e-08,
    2.9802322387695303e-08,
    1.4901161193847655e-08,
    7.450580596923828e-09,
    3.725290298461914e-09,
    1.862645149230957e-09,
    9.313225746154785e-10,
    4.656612873077393e-10,
]

# Product_i(1 / sqrt(1 + 2^(-2i))) for i in [0, 31]
CORDIC_K_32 = 0.6072529350088814
LN2 = 0.6931471805599453
HALF_LN2 = LN2 / 2.0
SQRT_HALF = 0.7071067811865476
SQRT_TWO = 1.4142135623730951


def _base_ops() -> Dict[str, int]:
    return {"add_sub": 0, "multiply": 0, "divide": 0}


def _merge_stats(*parts: Dict[str, int]) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for part in parts:
        for key, value in part.items():
            out[key] = out.get(key, 0) + value
    return out


def _horner(x: float, coeffs: Tuple[float, ...], stats: Dict[str, int]) -> float:
    out = coeffs[0]
    for c in coeffs[1:]:
        out = out * x + c
        stats["multiply"] += 1
        stats["add_sub"] += 1
    return out


def _normalize_angle(theta: float) -> Tuple[float, int]:
    """Map theta to [-pi/2, pi/2] and return sign correction for sin/cos."""
    theta = (theta + PI) % TWO_PI - PI
    sign = 1
    if theta > HALF_PI:
        theta -= PI
        sign = -1
    elif theta < -HALF_PI:
        theta += PI
        sign = -1
    return theta, sign


def cordic_sin_cos(
    theta: float, iterations: int = 32, return_stats: bool = False
) -> Tuple[float, float] | Tuple[Tuple[float, float], Dict[str, int]]:
    """Return (sin(theta), cos(theta)) using CORDIC rotation mode."""
    if iterations < 1 or iterations > len(CORDIC_ANGLES):
        raise ValueError(f"iterations must be in [1, {len(CORDIC_ANGLES)}]")
    theta_n, sign = _normalize_angle(theta)

    x = CORDIC_K_32
    y = 0.0
    z = theta_n

    for i in range(iterations):
        angle = CORDIC_ANGLES[i]
        di = 1.0 if z >= 0.0 else -1.0
        shift = 2.0**-i
        x_new = x - di * y * shift
        y_new = y + di * x * shift
        z -= di * angle
        x, y = x_new, y_new

    out = (sign * y, sign * x)
    if return_stats:
        stats = _base_ops()
        stats["cordic_rotation_loops"] = iterations
        stats["multiply"] += 5 * iterations + 2
        stats["add_sub"] += 3 * iterations
        return out, stats
    return out


def cordic_tan(
    theta: float, iterations: int = 32, return_stats: bool = False
) -> float | Tuple[float, Dict[str, int]]:
    (s, c), stats = cordic_sin_cos(theta, iterations=iterations, return_stats=True)
    if abs(c) < 1e-15:
        raise ZeroDivisionError("tan(theta) undefined because cos(theta) is too close to 0")
    out = s / c
    if return_stats:
        stats = dict(stats)
        stats["divide"] += 1
        return out, stats
    return out


def _floor(v: float) -> int:
    i = int(v)
    return i if i <= v else i - 1


def _sqrt_newton(v: float, iterations: int = 30, return_stats: bool = False) -> float | Tuple[float, Dict[str, int]]:
    if v < 0.0:
        raise ValueError("sqrt input must be non-negative")
    if v == 0.0:
        if return_stats:
            stats = _base_ops()
            stats["sqrt_newton_loops"] = 0
            return 0.0, stats
        return 0.0
    x = v if v >= 1.0 else 1.0
    for _ in range(iterations):
        x = 0.5 * (x + v / x)
    if return_stats:
        stats = _base_ops()
        stats["sqrt_newton_loops"] = iterations
        stats["divide"] += iterations
        stats["add_sub"] += iterations
        stats["multiply"] += iterations
        return x, stats
    return x


def cordic_atan(
    value: float, iterations: int = 32, return_stats: bool = False
) -> float | Tuple[float, Dict[str, int]]:
    """atan(value) using CORDIC vectoring mode."""
    if iterations < 1 or iterations > len(CORDIC_ANGLES):
        raise ValueError(f"iterations must be in [1, {len(CORDIC_ANGLES)}]")

    x = 1.0
    y = value
    z = 0.0

    for i in range(iterations):
        angle = CORDIC_ANGLES[i]
        di = 1.0 if y >= 0.0 else -1.0
        shift = 2.0**-i
        x_new = x + di * y * shift
        y_new = y - di * x * shift
        z += di * angle
        x, y = x_new, y_new

    if return_stats:
        stats = _base_ops()
        stats["cordic_vectoring_loops"] = iterations
        stats["multiply"] += 5 * iterations
        stats["add_sub"] += 3 * iterations
        return z, stats
    return z


def cordic_atan2(
    y: float, x: float, iterations: int = 32, return_stats: bool = False
) -> float | Tuple[float, Dict[str, int]]:
    if x > 0.0:
        out, stats = cordic_atan(y / x, iterations=iterations, return_stats=True)
        if return_stats:
            stats = dict(stats)
            stats["divide"] += 1
            return out, stats
        return out
    if x < 0.0 and y >= 0.0:
        out, stats = cordic_atan(y / x, iterations=iterations, return_stats=True)
        out += PI
        if return_stats:
            stats = dict(stats)
            stats["divide"] += 1
            stats["add_sub"] += 1
            return out, stats
        return out
    if x < 0.0 and y < 0.0:
        out, stats = cordic_atan(y / x, iterations=iterations, return_stats=True)
        out -= PI
        if return_stats:
            stats = dict(stats)
            stats["divide"] += 1
            stats["add_sub"] += 1
            return out, stats
        return out
    if x == 0.0 and y > 0.0:
        if return_stats:
            stats = _base_ops()
            stats["cordic_vectoring_loops"] = 0
            return HALF_PI, stats
        return HALF_PI
    if x == 0.0 and y < 0.0:
        if return_stats:
            stats = _base_ops()
            stats["cordic_vectoring_loops"] = 0
            return -HALF_PI, stats
        return -HALF_PI
    if return_stats:
        stats = _base_ops()
        stats["cordic_vectoring_loops"] = 0
        return 0.0, stats
    return 0.0


def cordic_asin(
    value: float, iterations: int = 32, return_stats: bool = False
) -> float | Tuple[float, Dict[str, int]]:
    if value < -1.0 or value > 1.0:
        raise ValueError("asin input must be in [-1, 1]")
    base = 1.0 - value * value
    root, sqrt_stats = _sqrt_newton(base, iterations=30, return_stats=True)
    out, atan_stats = cordic_atan2(value, root, iterations=iterations, return_stats=True)
    if return_stats:
        base_stats = _base_ops()
        base_stats["multiply"] += 1
        base_stats["add_sub"] += 1
        return out, _merge_stats(base_stats, sqrt_stats, atan_stats)
    return out


def cordic_acos(
    value: float, iterations: int = 32, return_stats: bool = False
) -> float | Tuple[float, Dict[str, int]]:
    asin_out = cordic_asin(value, iterations=iterations, return_stats=return_stats)
    if return_stats:
        out, stats = asin_out
        stats = dict(stats)
        stats["add_sub"] += 1
        return HALF_PI - out, stats
    return HALF_PI - asin_out


def basic_exp(x: float, terms: int = 40, return_stats: bool = False) -> float | Tuple[float, Dict[str, int]]:
    """exp(x) via range reduction + Taylor on small interval."""
    n = _floor(x / LN2)
    r = x - n * LN2

    term = 1.0
    s = 1.0
    used = 0
    for k in range(1, terms + 1):
        used += 1
        term *= r / k
        s += term
        if abs(term) < 1e-16:
            break

    out = s * (2.0**n)
    if return_stats:
        stats = _base_ops()
        stats["exp_series_loops"] = used
        stats["divide"] += 1 + used
        stats["multiply"] += 1 + used + 1
        stats["add_sub"] += 1 + used
        return out, stats
    return out


def basic_log(x: float, terms: int = 60, return_stats: bool = False) -> float | Tuple[float, Dict[str, int]]:
    """ln(x) via normalization + atanh series."""
    if x <= 0.0:
        raise ValueError("log input must be > 0")

    k = 0
    m = x
    norm_down = 0
    norm_up = 0
    while m > 1.5:
        m *= 0.5
        k += 1
        norm_down += 1
    while m < 0.75:
        m *= 2.0
        k -= 1
        norm_up += 1

    t = (m - 1.0) / (m + 1.0)
    t2 = t * t
    cur = t
    s = 0.0
    n = 1
    used = 0
    for _ in range(terms):
        used += 1
        s += cur / n
        cur *= t2
        n += 2

    out = 2.0 * s + k * LN2
    if return_stats:
        stats = _base_ops()
        stats["log_norm_down_loops"] = norm_down
        stats["log_norm_up_loops"] = norm_up
        stats["log_series_loops"] = used
        stats["multiply"] += norm_down + norm_up + 1 + used + 2
        stats["add_sub"] += norm_down + norm_up + 2 + used + used + 1
        stats["divide"] += 1 + used
        return out, stats
    return out


def minimax_sin_cos(
    theta: float, return_stats: bool = False
) -> Tuple[float, float] | Tuple[Tuple[float, float], Dict[str, int]]:
    """sin/cos via pre-fit minimax-style kernel polynomials."""
    stats = _base_ops()
    x = (theta + PI) % TWO_PI - PI

    sign_s = 1.0
    sign_c = 1.0
    if x > HALF_PI:
        x = PI - x
        sign_c = -1.0
        stats["add_sub"] += 1
    elif x < -HALF_PI:
        x = -PI - x
        sign_s = -1.0
        sign_c = -1.0
        stats["add_sub"] += 1

    z = x * x
    stats["multiply"] += 1

    # fdlibm-style kernel_sin coefficients
    sin_coeffs = (
        1.58969099521155010221e-10,
        -2.50507602534068634195e-08,
        2.75573137070700676789e-06,
        -1.98412698298579493134e-04,
        8.33333333332248946124e-03,
        -1.66666666666666324348e-01,
    )
    sin_poly = _horner(z, sin_coeffs, stats)
    sin_core = x + x * z * sin_poly
    stats["multiply"] += 2
    stats["add_sub"] += 1
    s = sign_s * sin_core
    if sign_s < 0:
        stats["multiply"] += 1

    # fdlibm-style kernel_cos coefficients
    cos_coeffs = (
        -1.13596475577881948265e-11,
        2.08757232129817482790e-09,
        -2.75573143513906633035e-07,
        2.48015872894767294178e-05,
        -1.38888888888741095749e-03,
        4.16666666666666019037e-02,
    )
    cos_poly = _horner(z, cos_coeffs, stats)
    cos_core = 1.0 - 0.5 * z + z * z * cos_poly
    stats["multiply"] += 3
    stats["add_sub"] += 2
    c = sign_c * cos_core
    if sign_c < 0:
        stats["multiply"] += 1

    if return_stats:
        stats["minimax_trig_loops"] = len(sin_coeffs) - 1 + len(cos_coeffs) - 1
        return (s, c), stats
    return s, c


def minimax_tan(theta: float, return_stats: bool = False) -> float | Tuple[float, Dict[str, int]]:
    (s, c), stats = minimax_sin_cos(theta, return_stats=True)
    if abs(c) < 1e-15:
        raise ZeroDivisionError("tan(theta) undefined because cos(theta) is too close to 0")
    out = s / c
    stats["divide"] += 1
    if return_stats:
        return out, stats
    return out


def minimax_atan(value: float, return_stats: bool = False) -> float | Tuple[float, Dict[str, int]]:
    """atan via range transform + odd minimax-style polynomial on [0,1]."""
    stats = _base_ops()
    sign = 1.0
    x = value
    if x < 0.0:
        x = -x
        sign = -1.0
        stats["multiply"] += 1

    invert = x > 1.0
    if invert:
        x = 1.0 / x
        stats["divide"] += 1

    z = x * x
    stats["multiply"] += 1
    # Pre-fit minimax-like coefficients on [0, 1] for atan(x)/x
    coeffs = (
        5.265332e-02,
        -1.1643287e-01,
        1.9354346e-01,
        -3.3262347e-01,
        9.9997726e-01,
    )
    poly = _horner(z, coeffs, stats)
    out = x * poly
    stats["multiply"] += 1

    if invert:
        out = HALF_PI - out
        stats["add_sub"] += 1
    if sign < 0.0:
        out = -out
        stats["multiply"] += 1

    if return_stats:
        stats["minimax_atan_loops"] = len(coeffs) - 1
        return out, stats
    return out


def minimax_asin(value: float, return_stats: bool = False) -> float | Tuple[float, Dict[str, int]]:
    if value < -1.0 or value > 1.0:
        raise ValueError("asin input must be in [-1, 1]")
    base_stats = _base_ops()
    base = 1.0 - value * value
    base_stats["multiply"] += 1
    base_stats["add_sub"] += 1
    root, sqrt_stats = _sqrt_newton(base, iterations=30, return_stats=True)
    atan_stats = _base_ops()
    x = value / root if root != 0.0 else (1e30 if value >= 0.0 else -1e30)
    atan_stats["divide"] += 1 if root != 0.0 else 0
    out, mini_stats = minimax_atan(x, return_stats=True)
    if return_stats:
        return out, _merge_stats(base_stats, sqrt_stats, atan_stats, mini_stats)
    return out


def minimax_acos(value: float, return_stats: bool = False) -> float | Tuple[float, Dict[str, int]]:
    asin_out = minimax_asin(value, return_stats=return_stats)
    if return_stats:
        out, stats = asin_out
        stats = dict(stats)
        stats["add_sub"] += 1
        return HALF_PI - out, stats
    return HALF_PI - asin_out


def minimax_exp(x: float, return_stats: bool = False) -> float | Tuple[float, Dict[str, int]]:
    stats = _base_ops()
    n = _floor(x / LN2)
    stats["divide"] += 1
    r = x - n * LN2
    stats["multiply"] += 1
    stats["add_sub"] += 1

    # Pre-fit minimax-like polynomial on [-ln2/2, ln2/2]
    coeffs = (
        1.3888949086377719e-03,
        8.333367984342196e-03,
        4.166662829120047e-02,
        1.666666676600111e-01,
        4.999999999998344e-01,
        1.000000000000000e+00,
        1.000000000000000e+00,
    )
    out = _horner(r, coeffs, stats)
    out *= 2.0**n
    stats["multiply"] += 1
    if return_stats:
        stats["minimax_exp_loops"] = len(coeffs) - 1
        return out, stats
    return out


def minimax_log(x: float, return_stats: bool = False) -> float | Tuple[float, Dict[str, int]]:
    if x <= 0.0:
        raise ValueError("log input must be > 0")
    stats = _base_ops()
    k = 0
    m = x
    down = 0
    up = 0
    while m > SQRT_TWO:
        m *= 0.5
        k += 1
        down += 1
        stats["multiply"] += 1
        stats["add_sub"] += 1
    while m < SQRT_HALF:
        m *= 2.0
        k -= 1
        up += 1
        stats["multiply"] += 1
        stats["add_sub"] += 1

    u = m - 1.0
    stats["add_sub"] += 1
    # Pre-fit minimax-like coeffs for ln(1+u) on small interval around 0
    coeffs = (
        5.72287066e-02,
        -1.367338464e-01,
        1.991076053e-01,
        -2.499922755e-01,
        3.333314528e-01,
        -5.000000000e-01,
        1.000000000e+00,
        0.0,
    )
    log_m = _horner(u, coeffs, stats)
    out = log_m + k * LN2
    stats["multiply"] += 1
    stats["add_sub"] += 1
    if return_stats:
        stats["log_norm_down_loops"] = down
        stats["log_norm_up_loops"] = up
        stats["minimax_log_loops"] = len(coeffs) - 1
        return out, stats
    return out


def _print_stats(stats: Dict[str, int]) -> None:
    if not stats:
        return
    print("stats:")
    for key, value in stats.items():
        print(f"  {key} = {value}")


def _print_func_stats(name: str, stats: Dict[str, int]) -> None:
    print(f"stats[{name}]:")
    for key, value in stats.items():
        print(f"  {key} = {value}")


def main() -> None:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--check", action="store_true", help="Use math only for verification")
    args, unknown = parser.parse_known_args()
    if unknown:
        print("Usage: python3 cordic.py [--check]")
        print("No positional args are required.")
        sys.exit(1)

    iterations = 32
    trig_deg = 30.0
    atan_in = 0.5
    asin_in = 0.5
    acos_in = 0.5
    exp_in = 1.2
    log_in = 1.2

    theta = trig_deg * PI / 180.0
    (s, c), trig_stats = cordic_sin_cos(theta, iterations=iterations, return_stats=True)
    tan_val = s / c
    tan_stats = dict(trig_stats)
    tan_stats["divide"] += 1

    print(f"angle(deg) = {trig_deg:.6f}")
    print(f"sin  CORDIC = {s:.16f}")
    _print_func_stats("sin", trig_stats)
    print(f"cos  CORDIC = {c:.16f}")
    _print_func_stats("cos", trig_stats)
    print(f"tan  CORDIC = {tan_val:.16f}")
    _print_func_stats("tan", tan_stats)

    atan_out, atan_stats = cordic_atan(atan_in, iterations=iterations, return_stats=True)
    print(f"atan CORDIC = {atan_out:.16f} (x={atan_in})")
    _print_func_stats("atan", atan_stats)

    asin_out, asin_stats = cordic_asin(asin_in, iterations=iterations, return_stats=True)
    print(f"asin CORDIC = {asin_out:.16f} (x={asin_in})")
    _print_func_stats("asin", asin_stats)

    acos_out, acos_stats = cordic_acos(acos_in, iterations=iterations, return_stats=True)
    print(f"acos CORDIC = {acos_out:.16f} (x={acos_in})")
    _print_func_stats("acos", acos_stats)

    exp_out, exp_stats = basic_exp(exp_in, return_stats=True)
    print(f"exp  BASIC  = {exp_out:.16f} (x={exp_in})")
    _print_func_stats("exp", exp_stats)

    log_out, log_stats = basic_log(log_in, return_stats=True)
    print(f"log  BASIC  = {log_out:.16f} (x={log_in})")
    _print_func_stats("log", log_stats)

    print("--- minimax polynomial method ---")
    (ms, mc), mtrig_stats = minimax_sin_cos(theta, return_stats=True)
    mt, mtan_stats = minimax_tan(theta, return_stats=True)
    print(f"sin  MINIMAX = {ms:.16f}")
    _print_func_stats("minimax_sin", mtrig_stats)
    print(f"cos  MINIMAX = {mc:.16f}")
    _print_func_stats("minimax_cos", mtrig_stats)
    print(f"tan  MINIMAX = {mt:.16f}")
    _print_func_stats("minimax_tan", mtan_stats)

    matan, matan_stats = minimax_atan(atan_in, return_stats=True)
    print(f"atan MINIMAX = {matan:.16f} (x={atan_in})")
    _print_func_stats("minimax_atan", matan_stats)

    masin, masin_stats = minimax_asin(asin_in, return_stats=True)
    print(f"asin MINIMAX = {masin:.16f} (x={asin_in})")
    _print_func_stats("minimax_asin", masin_stats)

    macos, macos_stats = minimax_acos(acos_in, return_stats=True)
    print(f"acos MINIMAX = {macos:.16f} (x={acos_in})")
    _print_func_stats("minimax_acos", macos_stats)

    mexp, mexp_stats = minimax_exp(exp_in, return_stats=True)
    print(f"exp  MINIMAX = {mexp:.16f} (x={exp_in})")
    _print_func_stats("minimax_exp", mexp_stats)

    mlog, mlog_stats = minimax_log(log_in, return_stats=True)
    print(f"log  MINIMAX = {mlog:.16f} (x={log_in})")
    _print_func_stats("minimax_log", mlog_stats)

    if args.check:
        import math

        print(f"sin  math   = {math.sin(theta):.16f} | err = {abs(s - math.sin(theta)):.3e}")
        print(f"cos  math   = {math.cos(theta):.16f} | err = {abs(c - math.cos(theta)):.3e}")
        print(f"tan  math   = {math.tan(theta):.16f} | err = {abs(tan_val - math.tan(theta)):.3e}")
        print(f"atan math   = {math.atan(atan_in):.16f} | err = {abs(atan_out - math.atan(atan_in)):.3e}")
        print(f"asin math   = {math.asin(asin_in):.16f} | err = {abs(asin_out - math.asin(asin_in)):.3e}")
        print(f"acos math   = {math.acos(acos_in):.16f} | err = {abs(acos_out - math.acos(acos_in)):.3e}")
        print(f"exp  math   = {math.exp(exp_in):.16f} | err = {abs(exp_out - math.exp(exp_in)):.3e}")
        print(f"log  math   = {math.log(log_in):.16f} | err = {abs(log_out - math.log(log_in)):.3e}")
        print(f"sin  minimax err = {abs(ms - math.sin(theta)):.3e}")
        print(f"cos  minimax err = {abs(mc - math.cos(theta)):.3e}")
        print(f"tan  minimax err = {abs(mt - math.tan(theta)):.3e}")
        print(f"atan minimax err = {abs(matan - math.atan(atan_in)):.3e}")
        print(f"asin minimax err = {abs(masin - math.asin(asin_in)):.3e}")
        print(f"acos minimax err = {abs(macos - math.acos(acos_in)):.3e}")
        print(f"exp  minimax err = {abs(mexp - math.exp(exp_in)):.3e}")
        print(f"log  minimax err = {abs(mlog - math.log(log_in)):.3e}")


if __name__ == "__main__":
    main()
