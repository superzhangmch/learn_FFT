#!/usr/bin/env python3

from __future__ import annotations

import random
import struct

from vm16.demos import soft64_add_bits_vm, soft64_div_bits_vm, soft64_mul_bits_vm


def f2u(x: float) -> int:
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def u2f(u: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", u & 0xFFFFFFFFFFFFFFFF))[0]


def is_nan_bits(u: int) -> bool:
    exp = (u >> 52) & 0x7FF
    frac = u & ((1 << 52) - 1)
    return exp == 0x7FF and frac != 0


def expected_div_bits(a: float, b: float) -> int:
    ua = f2u(a)
    ub = f2u(b)
    sa = (ua >> 63) & 1
    sb = (ub >> 63) & 1
    ea = (ua >> 52) & 0x7FF
    eb = (ub >> 52) & 0x7FF
    fa = ua & ((1 << 52) - 1)
    fb = ub & ((1 << 52) - 1)
    is_na = ea == 0x7FF and fa != 0
    is_nb = eb == 0x7FF and fb != 0
    is_ia = ea == 0x7FF and fa == 0
    is_ib = eb == 0x7FF and fb == 0
    is_za = ea == 0 and fa == 0
    is_zb = eb == 0 and fb == 0
    sign = sa ^ sb

    if is_na or is_nb:
        return f2u(float("nan"))
    if is_ia and is_ib:
        return f2u(float("nan"))
    if is_za and is_zb:
        return f2u(float("nan"))
    if is_ia:
        return f2u(-float("inf") if sign else float("inf"))
    if is_ib:
        return f2u(-0.0 if sign else 0.0)
    if is_zb:
        return f2u(-float("inf") if sign else float("inf"))
    if is_za:
        return f2u(-0.0 if sign else 0.0)
    return f2u(a / b)


def check_bits(op: str, got: int, want: int) -> bool:
    if is_nan_bits(got) and is_nan_bits(want):
        return True
    return got == want


def run_edge_tests() -> None:
    edge_bits = [
        0x0000000000000000,  # +0
        0x8000000000000000,  # -0
        0x0000000000000001,  # min subnormal
        0x000FFFFFFFFFFFFF,  # max subnormal
        0x0010000000000000,  # min normal
        0x7FEFFFFFFFFFFFFF,  # max finite
        0x3FF0000000000000,  # 1.0
        0xBFF0000000000000,  # -1.0
        0x7FF0000000000000,  # +inf
        0xFFF0000000000000,  # -inf
        0x7FF8000000000001,  # qNaN
    ]
    vals = [u2f(u) for u in edge_bits]

    total = 0
    for a in vals:
        for b in vals:
            total += 3
            got_add, _, _ = soft64_add_bits_vm(f2u(a), f2u(b))
            got_mul, _, _, _ = soft64_mul_bits_vm(f2u(a), f2u(b))
            got_div, _, _, _ = soft64_div_bits_vm(f2u(a), f2u(b))

            want_add = f2u(a + b)
            want_mul = f2u(a * b)
            want_div = expected_div_bits(a, b)

            if not check_bits("add", got_add, want_add):
                raise AssertionError(f"edge add mismatch: a={a} b={b} got=0x{got_add:016X} want=0x{want_add:016X}")
            if not check_bits("mul", got_mul, want_mul):
                raise AssertionError(f"edge mul mismatch: a={a} b={b} got=0x{got_mul:016X} want=0x{want_mul:016X}")
            if not check_bits("div", got_div, want_div):
                raise AssertionError(f"edge div mismatch: a={a} b={b} got=0x{got_div:016X} want=0x{want_div:016X}")
    print(f"edge_cases: PASS ({total} checks)")


def run_random_tests(n: int = 800) -> None:
    random.seed(42)
    bad_add = 0
    bad_mul = 0
    bad_div = 0
    for _ in range(n):
        ua = random.getrandbits(64)
        ub = random.getrandbits(64)
        a = u2f(ua)
        b = u2f(ub)

        got_add, _, _ = soft64_add_bits_vm(ua, ub)
        got_mul, _, _, _ = soft64_mul_bits_vm(ua, ub)
        got_div, _, _, _ = soft64_div_bits_vm(ua, ub)

        want_add = f2u(a + b)
        want_mul = f2u(a * b)
        want_div = expected_div_bits(a, b)

        if not check_bits("add", got_add, want_add):
            bad_add += 1
        if not check_bits("mul", got_mul, want_mul):
            bad_mul += 1
        if not check_bits("div", got_div, want_div):
            bad_div += 1

    print(f"random_add_mismatch={bad_add}/{n}")
    print(f"random_mul_mismatch={bad_mul}/{n}")
    print(f"random_div_mismatch={bad_div}/{n}")
    if bad_add or bad_mul or bad_div:
        raise AssertionError("random fp64 tests failed")


def main() -> None:
    run_edge_tests()
    run_random_tests()
    print("vm16_fp64_test: PASS")


if __name__ == "__main__":
    main()
