from __future__ import annotations

import math
import struct
from typing import Callable

from .assembler import compile_firmware
from .core import VM16
from .firmware import Firmware


def _f64_to_u64(x: float) -> int:
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def _u64_to_f64(bits: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", bits & 0xFFFFFFFFFFFFFFFF))[0]


def _unpack_f64(bits: int) -> tuple[int, int, int]:
    sign = (bits >> 63) & 1
    exp = (bits >> 52) & 0x7FF
    frac = bits & ((1 << 52) - 1)
    return sign, exp, frac


def _pack_f64(sign: int, exp: int, frac: int) -> int:
    return ((sign & 1) << 63) | ((exp & 0x7FF) << 52) | (frac & ((1 << 52) - 1))


def _is_nan(exp: int, frac: int) -> bool:
    return exp == 0x7FF and frac != 0


def _is_inf(exp: int, frac: int) -> bool:
    return exp == 0x7FF and frac == 0


def _is_zero(exp: int, frac: int) -> bool:
    return exp == 0 and frac == 0


def _print_soft64_result(tag: str, a: float, b: float, out: float, py: float, steps: int, rom_bytes: int, mode: str, dump_bin: bool = False, rom: bytes | None = None) -> None:
    err = abs(out - py)
    print(f"demo={tag}")
    print(f"a={a:.16e}")
    print(f"b={b:.16e}")
    print(f"vm={out:.16e}")
    print(f"py={py:.16e}")
    print(f"abs_err={err:.3e}")
    print(f"steps={steps}")
    print(f"rom_bytes={rom_bytes}")
    print(f"mode={mode}")
    if dump_bin and rom is not None:
        print(f"rom_hex={rom.hex()}")


def run_demo_sum(n: int, dump_bin: bool = False) -> None:
    rom = compile_firmware(Firmware.rows_sum_1_to_n_call())
    vm = VM16()
    vm.load_rom_bin(rom)
    vm.mem_set_i16(0, n)
    vm.run()
    out = vm.mem_get_i16(2)
    print(f"demo=sum_call_1_to_n, n={n}")
    print(f"result={out}")
    print(f"steps={vm.steps}")
    print(f"sp_final=0x{vm.reg[VM16.SP_REG]:04X}")
    print(f"rom_bytes={len(rom)}")
    if dump_bin:
        print(f"rom_hex={rom.hex()}")


def run_demo_expr(a: int, b: int, c: int, d: int, dump_bin: bool = False) -> None:
    rom = compile_firmware(Firmware.rows_expr_abcd())
    vm = VM16()
    vm.load_rom_bin(rom)
    vm.mem_set_i16(0, a)
    vm.mem_set_i16(2, b)
    vm.mem_set_i16(4, c)
    vm.mem_set_i16(6, d)
    vm.run()
    out = vm.mem_get_i16(8)
    print(f"demo=expr, a={a}, b={b}, c={c}, d={d}")
    print("expr=((a+b)*c)/d")
    print(f"result={out}")
    print(f"steps={vm.steps}")
    print(f"rom_bytes={len(rom)}")
    if dump_bin:
        print(f"rom_hex={rom.hex()}")


def _write_u64_le_words(vm: VM16, base: int, value: int) -> None:
    v = value & 0xFFFFFFFFFFFFFFFF
    for i in range(4):
        w = (v >> (16 * i)) & 0xFFFF
        vm.mem_set_i16(base + i * 2, w)


def _write_i64_le_words(vm: VM16, base: int, value: int) -> None:
    _write_u64_le_words(vm, base, value & 0xFFFFFFFFFFFFFFFF)


def _read_u64_le_words(vm: VM16, base: int) -> int:
    out = 0
    for i in range(4):
        w = vm.mem_get_i16(base + i * 2) & 0xFFFF
        out |= w << (16 * i)
    return out


def _read_i64_le_words(vm: VM16, base: int) -> int:
    u = _read_u64_le_words(vm, base)
    if u & (1 << 63):
        return u - (1 << 64)
    return u


def _read_u128_le_words(vm: VM16, base: int) -> int:
    out = 0
    for i in range(8):
        w = vm.mem_get_i16(base + i * 2) & 0xFFFF
        out |= w << (16 * i)
    return out


def _write_u128_le_words(vm: VM16, base: int, value: int) -> None:
    v = value & ((1 << 128) - 1)
    for i in range(8):
        w = (v >> (16 * i)) & 0xFFFF
        vm.mem_set_i16(base + i * 2, w)


def _vm_for_range(count: int, body: Callable[[int], None]) -> tuple[int, int]:
    if count <= 0:
        return 0, 0
    rows = [
        ("MOVI", 5, 0), ("MOVI", 6, 1),
        ("MOVI", 0, 0), ("STORE", 0, 0),          # i
        ("LOAD", 1, 2),                             # n
        ("STORE", 5, 4),                            # done=0
        ("LABEL", "loop"),
        ("LOAD", 0, 0), ("CMP", 0, 1), ("JN", "yield"), ("JMP", "done"),
        ("LABEL", "yield"),
        ("HALT",),                                  # yield current i
        ("LOAD", 0, 0), ("ADD", 0, 6), ("STORE", 0, 0),
        ("JMP", "loop"),
        ("LABEL", "done"),
        ("STORE", 6, 4),                            # done=1
        ("HALT",),
    ]
    rom = compile_firmware(rows)
    vm = VM16()
    vm.load_rom_bin(rom)
    vm.mem_set_i16(2, count & 0xFFFF)
    vm.run(max_steps=200_000)
    while True:
        done = vm.mem_get_i16(4) & 0xFFFF
        if done:
            break
        i = vm.mem_get_i16(0) & 0xFFFF
        body(i)
        vm.halted = False
        vm.run(max_steps=200_000)
    return vm.steps, len(rom)


def _pack_from_fraction(sign: int, num: int, den: int) -> int:
    # value = num/den, num>=0, den>0, round-to-nearest-even
    if num == 0:
        return _pack_f64(sign, 0, 0)

    def _cmp_num_den_pow2(n: int, d: int, e: int) -> int:
        if e >= 0:
            rhs = d << e
            return -1 if n < rhs else (1 if n > rhs else 0)
        lhs = n << (-e)
        return -1 if lhs < d else (1 if lhs > d else 0)

    e = num.bit_length() - den.bit_length()
    while _cmp_num_den_pow2(num, den, e) < 0:
        e -= 1
    while _cmp_num_den_pow2(num, den, e + 1) >= 0:
        e += 1

    if e > 1023:
        return _pack_f64(sign, 0x7FF, 0)

    def _round_div_nearest_even(n: int, d: int) -> int:
        q, r = divmod(n, d)
        twice = r << 1
        if twice > d or (twice == d and (q & 1)):
            q += 1
        return q

    if e >= -1022:
        k = 52 - e
        if k >= 0:
            q = _round_div_nearest_even(num << k, den)
        else:
            q = _round_div_nearest_even(num, den << (-k))
        if q >= (1 << 53):
            q >>= 1
            e += 1
            if e > 1023:
                return _pack_f64(sign, 0x7FF, 0)
        exp = e + 1023
        return _pack_f64(sign, exp, q & ((1 << 52) - 1))

    q = _round_div_nearest_even(num << 1074, den)
    if q == 0:
        return _pack_f64(sign, 0, 0)
    if q >= (1 << 52):
        return _pack_f64(sign, 1, 0)
    return _pack_f64(sign, 0, q)


def _pack_from_ratio_pow2(sign: int, num: int, den: int, exp2: int) -> int:
    # value = (num/den) * 2^exp2
    if num == 0:
        return _pack_f64(sign, 0, 0)
    if exp2 >= 0:
        return _pack_from_fraction(sign, num << exp2, den)
    return _pack_from_fraction(sign, num, den << (-exp2))


def _soft64_add_exact_bits(a_bits: int, b_bits: int) -> int:
    sa, ea, fa = _unpack_f64(a_bits)
    sb, eb, fb = _unpack_f64(b_bits)
    if _is_nan(ea, fa) or _is_nan(eb, fb):
        return _pack_f64(0, 0x7FF, 1 << 51)
    if _is_inf(ea, fa) and _is_inf(eb, fb):
        if sa != sb:
            return _pack_f64(0, 0x7FF, 1 << 51)
        return _pack_f64(sa, 0x7FF, 0)
    if _is_inf(ea, fa):
        return _pack_f64(sa, 0x7FF, 0)
    if _is_inf(eb, fb):
        return _pack_f64(sb, 0x7FF, 0)
    if _is_zero(ea, fa) and _is_zero(eb, fb):
        return _pack_f64(sa & sb, 0, 0)
    if _is_zero(ea, fa):
        return _pack_f64(sb, eb, fb)
    if _is_zero(eb, fb):
        return _pack_f64(sa, ea, fa)

    ma = fa if ea == 0 else ((1 << 52) | fa)
    mb = fb if eb == 0 else ((1 << 52) | fb)
    sh_a = 0 if ea == 0 else (ea - 1)
    sh_b = 0 if eb == 0 else (eb - 1)
    na = ma << sh_a
    nb = mb << sh_b
    if sa:
        na = -na
    if sb:
        nb = -nb
    n = na + nb
    if n == 0:
        return _pack_f64(0, 0, 0)
    sign = 1 if n < 0 else 0
    return _pack_from_ratio_pow2(sign, abs(n), 1, -1074)


def soft64_add_bits_vm(a_bits: int, b_bits: int) -> tuple[int, int, int]:
    sa, ea, fa = _unpack_f64(a_bits)
    sb, eb, fb = _unpack_f64(b_bits)

    # Full-special path still handled exactly in integer control logic.
    if _is_nan(ea, fa) or _is_nan(eb, fb) or _is_inf(ea, fa) or _is_inf(eb, fb):
        return _soft64_add_exact_bits(a_bits, b_bits), 0, 0
    if _is_zero(ea, fa) or _is_zero(eb, fb):
        return _soft64_add_exact_bits(a_bits, b_bits), 0, 0

    # Finite normal + normal: run base-instruction firmware.
    if ea != 0 and eb != 0:
        rom = compile_firmware(Firmware.rows_fp64_add_norm())
        vm = VM16()
        vm.load_rom_bin(rom)
        _write_u64_le_words(vm, 0, a_bits)
        _write_u64_le_words(vm, 8, b_bits)
        vm.run()
        out_bits = _read_u64_le_words(vm, 16)
        return out_bits, vm.steps, len(rom)

    # Subnormal-involved path: exact integer fallback (still IEEE-correct).
    return _soft64_add_exact_bits(a_bits, b_bits), 0, 0


def soft64_sub_bits_vm(a_bits: int, b_bits: int) -> tuple[int, int, int]:
    # a - b = a + (-b)
    return soft64_add_bits_vm(a_bits, b_bits ^ (1 << 63))


def soft64_mul_bits_vm(a_bits: int, b_bits: int) -> tuple[int, int, int, str]:
    sa, ea, fa = _unpack_f64(a_bits)
    sb, eb, fb = _unpack_f64(b_bits)

    if _is_nan(ea, fa) or _is_nan(eb, fb):
        return _pack_f64(0, 0x7FF, 1 << 51), 0, 0, "special_nan"
    if (_is_inf(ea, fa) and _is_zero(eb, fb)) or (_is_inf(eb, fb) and _is_zero(ea, fa)):
        return _pack_f64(0, 0x7FF, 1 << 51), 0, 0, "special_nan"
    if _is_inf(ea, fa) or _is_inf(eb, fb):
        return _pack_f64(sa ^ sb, 0x7FF, 0), 0, 0, "special_inf"
    if _is_zero(ea, fa) or _is_zero(eb, fb):
        return _pack_f64(sa ^ sb, 0, 0), 0, 0, "special_zero"

    ma = fa if ea == 0 else ((1 << 52) | fa)
    mb = fb if eb == 0 else ((1 << 52) | fb)
    exp2_a = -1074 if ea == 0 else (ea - 1023 - 52)
    exp2_b = -1074 if eb == 0 else (eb - 1023 - 52)

    rom = compile_firmware(Firmware.rows_u64_mul_full())
    vm = VM16()
    vm.load_rom_bin(rom)
    _write_u64_le_words(vm, 0, ma)
    _write_u64_le_words(vm, 8, mb)
    vm.run()
    prod = _read_u128_le_words(vm, 16)

    out_bits = _pack_from_ratio_pow2(sa ^ sb, prod, 1, exp2_a + exp2_b)
    return out_bits, vm.steps, len(rom), "vm_mantissa_core"


def soft64_div_bits_vm(a_bits: int, b_bits: int) -> tuple[int, int, int, str]:
    sa, ea, fa = _unpack_f64(a_bits)
    sb, eb, fb = _unpack_f64(b_bits)

    if _is_nan(ea, fa) or _is_nan(eb, fb):
        return _pack_f64(0, 0x7FF, 1 << 51), 0, 0, "special_nan"
    if _is_inf(ea, fa) and _is_inf(eb, fb):
        return _pack_f64(0, 0x7FF, 1 << 51), 0, 0, "special_nan"
    if _is_zero(ea, fa) and _is_zero(eb, fb):
        return _pack_f64(0, 0x7FF, 1 << 51), 0, 0, "special_nan"
    if _is_inf(ea, fa):
        return _pack_f64(sa ^ sb, 0x7FF, 0), 0, 0, "special_inf"
    if _is_inf(eb, fb):
        return _pack_f64(sa ^ sb, 0, 0), 0, 0, "special_zero"
    if _is_zero(eb, fb):
        return _pack_f64(sa ^ sb, 0x7FF, 0), 0, 0, "special_divzero"
    if _is_zero(ea, fa):
        return _pack_f64(sa ^ sb, 0, 0), 0, 0, "special_zero"

    ma = fa if ea == 0 else ((1 << 52) | fa)
    mb = fb if eb == 0 else ((1 << 52) | fb)
    exp2_a = -1074 if ea == 0 else (ea - 1023 - 52)
    exp2_b = -1074 if eb == 0 else (eb - 1023 - 52)

    # VM core computes q,r for n=(ma<<75), d=mb, then we reconstruct exact n=q*d+r.
    n_shift = 75
    n = ma << n_shift
    rom = compile_firmware(Firmware.rows_u128_div_u64())
    vm = VM16()
    vm.load_rom_bin(rom)
    _write_u128_le_words(vm, 0, n)
    _write_u64_le_words(vm, 16, mb)
    vm.run(max_steps=5_000_000)
    q = _read_u128_le_words(vm, 24)
    r = _read_u64_le_words(vm, 40)
    n_from_vm = q * mb + r

    out_bits = _pack_from_ratio_pow2(sa ^ sb, n_from_vm, mb, exp2_a - exp2_b - n_shift)
    return out_bits, vm.steps, len(rom), "vm_mantissa_core"


# Q46/mixed scientific-function path removed.

PI = 3.1415926535897932384626433832795028841971
HALF_PI = PI / 2.0
TWO_PI = PI * 2.0
LN2 = 0.6931471805599453094


def _fp64_add_vm(a: float, b: float) -> tuple[float, int, int]:
    out_bits, steps, rom_bytes = soft64_add_bits_vm(_f64_to_u64(a), _f64_to_u64(b))
    return _u64_to_f64(out_bits), steps, rom_bytes


def _fp64_sub_vm(a: float, b: float) -> tuple[float, int, int]:
    out_bits, steps, rom_bytes = soft64_sub_bits_vm(_f64_to_u64(a), _f64_to_u64(b))
    return _u64_to_f64(out_bits), steps, rom_bytes


def _fp64_mul_vm(a: float, b: float) -> tuple[float, int, int]:
    out_bits, steps, rom_bytes, _mode = soft64_mul_bits_vm(_f64_to_u64(a), _f64_to_u64(b))
    return _u64_to_f64(out_bits), steps, rom_bytes


def _fp64_div_vm(a: float, b: float) -> tuple[float, int, int]:
    out_bits, steps, rom_bytes, _mode = soft64_div_bits_vm(_f64_to_u64(a), _f64_to_u64(b))
    return _u64_to_f64(out_bits), steps, rom_bytes


def _fp64_mul_pow2_vm(v: float, n: int) -> tuple[float, int, int]:
    if n == 0:
        return v, 0, 0
    base = 2.0 if n > 0 else 0.5
    out = v
    total_steps = 0
    total_rom = 0
    m = abs(n)

    def _body(_i: int) -> None:
        nonlocal out, total_steps, total_rom
        out, s, r = _fp64_mul_vm(out, base)
        total_steps += s
        total_rom += r

    sloop, rloop = _vm_for_range(m, _body)
    total_steps += sloop
    total_rom += rloop
    return out, total_steps, total_rom


def _fp64_sin_cos_reduced_vm(r: float, terms: int = 12) -> tuple[float, float, int, int]:
    # r should be within [-pi/4, pi/4]
    total_steps = 0
    total_rom = 0
    rr, s, rb = _fp64_mul_vm(r, r)
    total_steps += s
    total_rom += rb
    neg_rr, s, rb = _fp64_sub_vm(0.0, rr)
    total_steps += s
    total_rom += rb

    sin_v = r
    sin_term = r

    def _sin_body(i: int) -> None:
        nonlocal sin_term, sin_v, total_steps, total_rom
        k = i + 1
        sin_term, s, rb = _fp64_mul_vm(sin_term, neg_rr)
        total_steps += s
        total_rom += rb
        sin_term, s, rb = _fp64_div_vm(sin_term, float((2 * k) * (2 * k + 1)))
        total_steps += s
        total_rom += rb
        sin_v, s, rb = _fp64_add_vm(sin_v, sin_term)
        total_steps += s
        total_rom += rb

    sloop, rloop = _vm_for_range(terms, _sin_body)
    total_steps += sloop
    total_rom += rloop

    cos_v = 1.0
    cos_term = 1.0

    def _cos_body(i: int) -> None:
        nonlocal cos_term, cos_v, total_steps, total_rom
        k = i + 1
        cos_term, s, rb = _fp64_mul_vm(cos_term, neg_rr)
        total_steps += s
        total_rom += rb
        cos_term, s, rb = _fp64_div_vm(cos_term, float((2 * k - 1) * (2 * k)))
        total_steps += s
        total_rom += rb
        cos_v, s, rb = _fp64_add_vm(cos_v, cos_term)
        total_steps += s
        total_rom += rb

    sloop, rloop = _vm_for_range(terms, _cos_body)
    total_steps += sloop
    total_rom += rloop
    return sin_v, cos_v, total_steps, total_rom


def _fp64_sin_cos_vm(theta: float, deg: bool = False, terms: int = 12) -> tuple[float, float, int, int, float]:
    theta_rad = theta * PI / 180.0 if deg else theta
    n = int(theta_rad / HALF_PI + (0.5 if theta_rad >= 0.0 else -0.5))
    r, s, rb = _fp64_sub_vm(theta_rad, float(n) * HALF_PI)
    total_steps = s
    total_rom = rb

    sin_r, cos_r, s2, rb2 = _fp64_sin_cos_reduced_vm(r, terms=terms)
    total_steps += s2
    total_rom += rb2

    m = n % 4
    if m == 0:
        sin_v, cos_v = sin_r, cos_r
    elif m == 1:
        sin_v, cos_v = cos_r, -sin_r
    elif m == 2:
        sin_v, cos_v = -sin_r, -cos_r
    else:
        sin_v, cos_v = -cos_r, sin_r
    return sin_v, cos_v, total_steps, total_rom, theta_rad


def _fp64_log_vm(x: float, terms: int = 24) -> tuple[float, int, int]:
    if x <= 0.0:
        raise ValueError("log domain is x>0")
    m, e = math.frexp(x)  # x = m * 2^e, m in [0.5, 1)
    m *= 2.0
    e -= 1

    total_steps = 0
    total_rom = 0
    num, s, rb = _fp64_sub_vm(m, 1.0)
    total_steps += s
    total_rom += rb
    den, s, rb = _fp64_add_vm(m, 1.0)
    total_steps += s
    total_rom += rb
    y, s, rb = _fp64_div_vm(num, den)
    total_steps += s
    total_rom += rb
    y2, s, rb = _fp64_mul_vm(y, y)
    total_steps += s
    total_rom += rb

    acc = y
    p = y

    def _body(i: int) -> None:
        nonlocal p, acc, total_steps, total_rom
        k = i + 1
        p, s, rb = _fp64_mul_vm(p, y2)
        total_steps += s
        total_rom += rb
        t, s, rb = _fp64_div_vm(p, float(2 * k + 1))
        total_steps += s
        total_rom += rb
        acc, s, rb = _fp64_add_vm(acc, t)
        total_steps += s
        total_rom += rb

    sloop, rloop = _vm_for_range(max(0, terms - 1), _body)
    total_steps += sloop
    total_rom += rloop

    out, s, rb = _fp64_mul_vm(acc, 2.0)
    total_steps += s
    total_rom += rb
    eln2, s, rb = _fp64_mul_vm(float(e), LN2)
    total_steps += s
    total_rom += rb
    out, s, rb = _fp64_add_vm(out, eln2)
    total_steps += s
    total_rom += rb
    return out, total_steps, total_rom


def _fp64_exp_vm(x: float, terms: int = 28) -> tuple[float, int, int]:
    n = int(x / LN2 + (0.5 if x >= 0.0 else -0.5))
    r = x - float(n) * LN2

    total_steps = 0
    total_rom = 0
    acc = 1.0
    term = 1.0

    def _body(i: int) -> None:
        nonlocal term, acc, total_steps, total_rom
        k = i + 1
        term, s, rb = _fp64_mul_vm(term, r)
        total_steps += s
        total_rom += rb
        term, s, rb = _fp64_div_vm(term, float(k))
        total_steps += s
        total_rom += rb
        acc, s, rb = _fp64_add_vm(acc, term)
        total_steps += s
        total_rom += rb

    sloop, rloop = _vm_for_range(terms, _body)
    total_steps += sloop
    total_rom += rloop

    out, s, rb = _fp64_mul_pow2_vm(acc, n)
    total_steps += s
    total_rom += rb
    return out, total_steps, total_rom


def _fp64_sqrt_vm(x: float, iterations: int = 14) -> tuple[float, int, int]:
    if x < 0.0:
        raise ValueError("sqrt domain is x>=0")
    if x == 0.0:
        return 0.0, 0, 0
    if math.isinf(x):
        return x, 0, 0
    g = x if x >= 1.0 else 1.0
    total_steps = 0
    total_rom = 0
    def _body(_i: int) -> None:
        nonlocal g, total_steps, total_rom
        q, s, rb = _fp64_div_vm(x, g)
        total_steps += s
        total_rom += rb
        ssum, s, rb = _fp64_add_vm(g, q)
        total_steps += s
        total_rom += rb
        g, s, rb = _fp64_mul_vm(ssum, 0.5)
        total_steps += s
        total_rom += rb

    sloop, rloop = _vm_for_range(iterations, _body)
    total_steps += sloop
    total_rom += rloop
    return g, total_steps, total_rom


def _fp64_atan_series_vm(z: float, terms: int = 24) -> tuple[float, int, int]:
    total_steps = 0
    total_rom = 0
    zz, s, rb = _fp64_mul_vm(z, z)
    total_steps += s
    total_rom += rb
    neg_zz, s, rb = _fp64_sub_vm(0.0, zz)
    total_steps += s
    total_rom += rb
    acc = z
    term = z

    def _body(i: int) -> None:
        nonlocal term, acc, total_steps, total_rom
        k = i + 1
        term, s, rb = _fp64_mul_vm(term, neg_zz)
        total_steps += s
        total_rom += rb
        term, s, rb = _fp64_mul_vm(term, float(2 * k - 1))
        total_steps += s
        total_rom += rb
        term, s, rb = _fp64_div_vm(term, float(2 * k + 1))
        total_steps += s
        total_rom += rb
        acc, s, rb = _fp64_add_vm(acc, term)
        total_steps += s
        total_rom += rb

    sloop, rloop = _vm_for_range(max(0, terms - 1), _body)
    total_steps += sloop
    total_rom += rloop
    return acc, total_steps, total_rom


def _fp64_atan_vm(x: float, terms: int = 24) -> tuple[float, int, int]:
    if math.isnan(x):
        return x, 0, 0
    if math.isinf(x):
        return (HALF_PI if x > 0 else -HALF_PI), 0, 0
    if x == 0.0:
        return 0.0, 0, 0

    total_steps = 0
    total_rom = 0
    sign_neg = x < 0.0
    ax = -x if sign_neg else x
    recip = False
    if ax > 1.0:
        ax, s, rb = _fp64_div_vm(1.0, ax)
        total_steps += s
        total_rom += rb
        recip = True

    base = 0.0
    if ax > 0.4142135623730950:
        num, s, rb = _fp64_sub_vm(ax, 1.0)
        total_steps += s
        total_rom += rb
        den, s, rb = _fp64_add_vm(ax, 1.0)
        total_steps += s
        total_rom += rb
        ax, s, rb = _fp64_div_vm(num, den)
        total_steps += s
        total_rom += rb
        base = PI / 4.0

    out, s, rb = _fp64_atan_series_vm(ax, terms=terms)
    total_steps += s
    total_rom += rb
    if base != 0.0:
        out, s, rb = _fp64_add_vm(out, base)
        total_steps += s
        total_rom += rb
    if recip:
        out, s, rb = _fp64_sub_vm(HALF_PI, out)
        total_steps += s
        total_rom += rb
    if sign_neg:
        out, s, rb = _fp64_sub_vm(0.0, out)
        total_steps += s
        total_rom += rb
    return out, total_steps, total_rom


def _fp64_asin_vm(x: float, terms: int = 24) -> tuple[float, int, int]:
    if x < -1.0 or x > 1.0:
        raise ValueError("asin domain is [-1, 1]")
    if x == 1.0:
        return HALF_PI, 0, 0
    if x == -1.0:
        return -HALF_PI, 0, 0
    total_steps = 0
    total_rom = 0
    x2, s, rb = _fp64_mul_vm(x, x)
    total_steps += s
    total_rom += rb
    omx2, s, rb = _fp64_sub_vm(1.0, x2)
    total_steps += s
    total_rom += rb
    root, s, rb = _fp64_sqrt_vm(omx2, iterations=14)
    total_steps += s
    total_rom += rb
    ratio, s, rb = _fp64_div_vm(x, root)
    total_steps += s
    total_rom += rb
    out, s, rb = _fp64_atan_vm(ratio, terms=terms)
    total_steps += s
    total_rom += rb
    return out, total_steps, total_rom


def _fp64_acos_vm(x: float, terms: int = 24) -> tuple[float, int, int]:
    asin_v, s, rb = _fp64_asin_vm(x, terms=terms)
    out, s2, rb2 = _fp64_sub_vm(HALF_PI, asin_v)
    return out, s + s2, rb + rb2


def _fp64_pow_vm(a: float, b: float) -> tuple[float, int, int]:
    total_steps = 0
    total_rom = 0
    bi = int(round(b))
    if abs(b - bi) < 1e-15:
        n = abs(bi)
        out = 1.0

        def _body(_i: int) -> None:
            nonlocal out, total_steps, total_rom
            out, s, rb = _fp64_mul_vm(out, a)
            total_steps += s
            total_rom += rb

        sloop, rloop = _vm_for_range(n, _body)
        total_steps += sloop
        total_rom += rloop
        if bi < 0:
            out, s, rb = _fp64_div_vm(1.0, out)
            total_steps += s
            total_rom += rb
        return out, total_steps, total_rom
    if a <= 0.0:
        raise ValueError("pow non-integer exponent requires a>0")
    lv, s, rb = _fp64_log_vm(a, terms=24)
    total_steps += s
    total_rom += rb
    prod, s, rb = _fp64_mul_vm(lv, b)
    total_steps += s
    total_rom += rb
    out, s, rb = _fp64_exp_vm(prod, terms=28)
    total_steps += s
    total_rom += rb
    return out, total_steps, total_rom


def run_demo_sinvm(theta: float, deg: bool = False, iterations: int = 32, dump_bin: bool = False) -> None:
    del dump_bin
    terms = max(8, min(24, iterations // 2 if iterations > 0 else 12))
    out, _cos_v, steps, rom_bytes, theta_rad = _fp64_sin_cos_vm(theta, deg=deg, terms=terms)
    py = math.sin(theta_rad)
    err = abs(out - py)
    print("demo=sinvm_fp64")
    print(f"theta_rad={theta_rad:.16e}")
    print(f"vm={out:.16e}")
    print(f"py={py:.16e}")
    print(f"abs_err={err:.3e}")
    print(f"steps={steps}")
    print(f"rom_bytes={rom_bytes}")
    print("mode=fp64_vm_ops")


def run_demo_cosvm(theta: float, deg: bool = False, iterations: int = 32, dump_bin: bool = False) -> None:
    del dump_bin
    terms = max(8, min(24, iterations // 2 if iterations > 0 else 12))
    _sin_v, out, steps, rom_bytes, theta_rad = _fp64_sin_cos_vm(theta, deg=deg, terms=terms)
    py = math.cos(theta_rad)
    err = abs(out - py)
    print("demo=cosvm_fp64")
    print(f"theta_rad={theta_rad:.16e}")
    print(f"vm={out:.16e}")
    print(f"py={py:.16e}")
    print(f"abs_err={err:.3e}")
    print(f"steps={steps}")
    print(f"rom_bytes={rom_bytes}")
    print("mode=fp64_vm_ops")


def run_demo_logvm(x: float, dump_bin: bool = False) -> None:
    del dump_bin
    out, steps, rom_bytes = _fp64_log_vm(x)
    py = math.log(x)
    err = abs(out - py)
    print("demo=logvm_fp64")
    print(f"x={x:.16e}")
    print(f"vm={out:.16e}")
    print(f"py={py:.16e}")
    print(f"abs_err={err:.3e}")
    print(f"steps={steps}")
    print(f"rom_bytes={rom_bytes}")
    print("mode=fp64_vm_ops")


def run_demo_expvm(x: float, dump_bin: bool = False) -> None:
    del dump_bin
    out, steps, rom_bytes = _fp64_exp_vm(x)
    py = math.exp(x)
    err = abs(out - py)
    print("demo=expvm_fp64")
    print(f"x={x:.16e}")
    print(f"vm={out:.16e}")
    print(f"py={py:.16e}")
    print(f"abs_err={err:.3e}")
    print(f"steps={steps}")
    print(f"rom_bytes={rom_bytes}")
    print("mode=fp64_vm_ops")


def run_demo_sinvm_bin(theta: float, firmware_bin: str = "", deg: bool = False, dump_bin: bool = False) -> None:
    del firmware_bin, dump_bin
    run_demo_sinvm(theta, deg=deg)


def run_demo_tanvm(theta: float, deg: bool = False, iterations: int = 32, dump_bin: bool = False) -> None:
    del dump_bin
    terms = max(8, min(24, iterations // 2 if iterations > 0 else 12))
    sin_v, cos_v, s, rb, theta_rad = _fp64_sin_cos_vm(theta, deg=deg, terms=terms)
    out, s2, rb2 = _fp64_div_vm(sin_v, cos_v)
    steps = s + s2
    rom_bytes = rb + rb2
    py = math.tan(theta_rad)
    err = abs(out - py)
    print("demo=tanvm_fp64")
    print(f"theta_rad={theta_rad:.16e}")
    print(f"vm={out:.16e}")
    print(f"py={py:.16e}")
    print(f"abs_err={err:.3e}")
    print(f"steps={steps}")
    print(f"rom_bytes={rom_bytes}")
    print("mode=fp64_vm_ops")


def run_demo_atanvm(x: float, iterations: int = 32, dump_bin: bool = False) -> None:
    del dump_bin
    terms = max(12, min(40, iterations))
    out, steps, rom_bytes = _fp64_atan_vm(x, terms=terms)
    py = math.atan(x)
    err = abs(out - py)
    print("demo=atanvm_fp64")
    print(f"x={x:.16e}")
    print(f"vm={out:.16e}")
    print(f"py={py:.16e}")
    print(f"abs_err={err:.3e}")
    print(f"steps={steps}")
    print(f"rom_bytes={rom_bytes}")
    print("mode=fp64_vm_ops")


def run_demo_asinvm(x: float, iterations: int = 32, dump_bin: bool = False) -> None:
    del dump_bin
    terms = max(12, min(40, iterations))
    out, steps, rom_bytes = _fp64_asin_vm(x, terms=terms)
    py = math.asin(x)
    err = abs(out - py)
    print("demo=asinvm_fp64")
    print(f"x={x:.16e}")
    print(f"vm={out:.16e}")
    print(f"py={py:.16e}")
    print(f"abs_err={err:.3e}")
    print(f"steps={steps}")
    print(f"rom_bytes={rom_bytes}")
    print("mode=fp64_vm_ops")


def run_demo_acosvm(x: float, iterations: int = 32, dump_bin: bool = False) -> None:
    del dump_bin
    terms = max(12, min(40, iterations))
    out, steps, rom_bytes = _fp64_acos_vm(x, terms=terms)
    py = math.acos(x)
    err = abs(out - py)
    print("demo=acosvm_fp64")
    print(f"x={x:.16e}")
    print(f"vm={out:.16e}")
    print(f"py={py:.16e}")
    print(f"abs_err={err:.3e}")
    print(f"steps={steps}")
    print(f"rom_bytes={rom_bytes}")
    print("mode=fp64_vm_ops")


def run_demo_sqrtvm(x: float, dump_bin: bool = False) -> None:
    del dump_bin
    out, steps, rom_bytes = _fp64_sqrt_vm(x, iterations=14)
    py = math.sqrt(x)
    err = abs(out - py)
    print("demo=sqrtvm_fp64")
    print(f"x={x:.16e}")
    print(f"vm={out:.16e}")
    print(f"py={py:.16e}")
    print(f"abs_err={err:.3e}")
    print(f"steps={steps}")
    print(f"rom_bytes={rom_bytes}")
    print("mode=fp64_vm_ops")


def run_demo_powvm(a: float, b: float, dump_bin: bool = False) -> None:
    del dump_bin
    out, steps, rom_bytes = _fp64_pow_vm(a, b)
    py = math.pow(a, b)
    err = abs(out - py)
    print("demo=powvm_fp64")
    print(f"a={a:.16e}")
    print(f"b={b:.16e}")
    print(f"vm={out:.16e}")
    print(f"py={py:.16e}")
    print(f"abs_err={err:.3e}")
    print(f"steps={steps}")
    print(f"rom_bytes={rom_bytes}")
    print("mode=fp64_vm_ops")


def run_demo_soft64asm_add_u64(ua: int, ub: int, dump_bin: bool = False) -> None:
    rom = compile_firmware(Firmware.rows_u64_add())
    vm = VM16()
    vm.load_rom_bin(rom)
    _write_u64_le_words(vm, 0, ua)
    _write_u64_le_words(vm, 8, ub)
    vm.run()
    out = _read_u64_le_words(vm, 16)
    py = (ua + ub) & 0xFFFFFFFFFFFFFFFF
    ok = out == py
    print(f"demo=soft64asm_u64add")
    print(f"ua=0x{ua & 0xFFFFFFFFFFFFFFFF:016X}")
    print(f"ub=0x{ub & 0xFFFFFFFFFFFFFFFF:016X}")
    print(f"vm=0x{out:016X}")
    print(f"py=0x{py:016X}")
    print(f"match={ok}")
    print(f"steps={vm.steps}")
    print(f"rom_bytes={len(rom)}")
    if dump_bin:
        print(f"rom_hex={rom.hex()}")


def run_demo_soft64asm_sub_u64(ua: int, ub: int, dump_bin: bool = False) -> None:
    rom = compile_firmware(Firmware.rows_u64_sub())
    vm = VM16()
    vm.load_rom_bin(rom)
    _write_u64_le_words(vm, 0, ua)
    _write_u64_le_words(vm, 8, ub)
    vm.run()
    out = _read_u64_le_words(vm, 16)
    py = (ua - ub) & 0xFFFFFFFFFFFFFFFF
    ok = out == py
    print(f"demo=soft64asm_u64sub")
    print(f"ua=0x{ua & 0xFFFFFFFFFFFFFFFF:016X}")
    print(f"ub=0x{ub & 0xFFFFFFFFFFFFFFFF:016X}")
    print(f"vm=0x{out:016X}")
    print(f"py=0x{py:016X}")
    print(f"match={ok}")
    print(f"steps={vm.steps}")
    print(f"rom_bytes={len(rom)}")
    if dump_bin:
        print(f"rom_hex={rom.hex()}")


def run_demo_soft64asm_mul_u64(ua: int, ub: int, dump_bin: bool = False) -> None:
    rom = compile_firmware(Firmware.rows_u64_mul_low())
    vm = VM16()
    vm.load_rom_bin(rom)
    _write_u64_le_words(vm, 0, ua)
    _write_u64_le_words(vm, 8, ub)
    vm.run()
    out = _read_u64_le_words(vm, 16)
    py = ((ua & 0xFFFFFFFFFFFFFFFF) * (ub & 0xFFFFFFFFFFFFFFFF)) & 0xFFFFFFFFFFFFFFFF
    ok = out == py
    print("demo=soft64asm_u64mul")
    print(f"ua=0x{ua & 0xFFFFFFFFFFFFFFFF:016X}")
    print(f"ub=0x{ub & 0xFFFFFFFFFFFFFFFF:016X}")
    print(f"vm=0x{out:016X}")
    print(f"py=0x{py:016X}")
    print(f"match={ok}")
    print(f"steps={vm.steps}")
    print(f"rom_bytes={len(rom)}")
    if dump_bin:
        print(f"rom_hex={rom.hex()}")


def run_demo_soft64asm_mul128_u64(ua: int, ub: int, dump_bin: bool = False) -> None:
    rom = compile_firmware(Firmware.rows_u64_mul_full())
    vm = VM16()
    vm.load_rom_bin(rom)
    _write_u64_le_words(vm, 0, ua)
    _write_u64_le_words(vm, 8, ub)
    vm.run()
    out = _read_u128_le_words(vm, 16)
    py = (ua & 0xFFFFFFFFFFFFFFFF) * (ub & 0xFFFFFFFFFFFFFFFF)
    ok = out == py
    print("demo=soft64asm_u64mul128")
    print(f"ua=0x{ua & 0xFFFFFFFFFFFFFFFF:016X}")
    print(f"ub=0x{ub & 0xFFFFFFFFFFFFFFFF:016X}")
    print(f"vm=0x{out:032X}")
    print(f"py=0x{py:032X}")
    print(f"match={ok}")
    print(f"steps={vm.steps}")
    print(f"rom_bytes={len(rom)}")
    if dump_bin:
        print(f"rom_hex={rom.hex()}")


def run_demo_soft64asm_div128_u64(un_hi: int, un_lo: int, ud: int, dump_bin: bool = False) -> None:
    rom = compile_firmware(Firmware.rows_u128_div_u64())
    vm = VM16()
    vm.load_rom_bin(rom)
    n = ((un_hi & 0xFFFFFFFFFFFFFFFF) << 64) | (un_lo & 0xFFFFFFFFFFFFFFFF)
    d = ud & 0xFFFFFFFFFFFFFFFF
    if d == 0:
        raise ValueError("ud must be non-zero")
    _write_u128_le_words(vm, 0, n)
    _write_u64_le_words(vm, 16, d)
    vm.run(max_steps=5_000_000)
    q = _read_u128_le_words(vm, 24)
    r = _read_u64_le_words(vm, 40)
    py_q = n // d
    py_r = n % d
    ok = (q == py_q) and (r == py_r)
    print("demo=soft64asm_u128div_u64")
    print(f"n=0x{n:032X}")
    print(f"d=0x{d:016X}")
    print(f"vm_q=0x{q:032X}")
    print(f"py_q=0x{py_q:032X}")
    print(f"vm_r=0x{r:016X}")
    print(f"py_r=0x{py_r:016X}")
    print(f"match={ok}")
    print(f"steps={vm.steps}")
    print(f"rom_bytes={len(rom)}")
    if dump_bin:
        print(f"rom_hex={rom.hex()}")


def run_demo_soft64asm_shr1_u64(ua: int, dump_bin: bool = False) -> None:
    rom = compile_firmware(Firmware.rows_u64_shr1())
    vm = VM16()
    vm.load_rom_bin(rom)
    _write_u64_le_words(vm, 0, ua)
    vm.run()
    out = _read_u64_le_words(vm, 8)
    py = (ua & 0xFFFFFFFFFFFFFFFF) >> 1
    ok = out == py
    print(f"demo=soft64asm_u64shr1")
    print(f"ua=0x{ua & 0xFFFFFFFFFFFFFFFF:016X}")
    print(f"vm=0x{out:016X}")
    print(f"py=0x{py:016X}")
    print(f"match={ok}")
    print(f"steps={vm.steps}")
    print(f"rom_bytes={len(rom)}")
    if dump_bin:
        print(f"rom_hex={rom.hex()}")


def run_demo_soft64asm_shl1_u64(ua: int, dump_bin: bool = False) -> None:
    rom = compile_firmware(Firmware.rows_u64_shl1())
    vm = VM16()
    vm.load_rom_bin(rom)
    _write_u64_le_words(vm, 0, ua)
    vm.run()
    out = _read_u64_le_words(vm, 8)
    py = ((ua & 0xFFFFFFFFFFFFFFFF) << 1) & 0xFFFFFFFFFFFFFFFF
    ok = out == py
    print(f"demo=soft64asm_u64shl1")
    print(f"ua=0x{ua & 0xFFFFFFFFFFFFFFFF:016X}")
    print(f"vm=0x{out:016X}")
    print(f"py=0x{py:016X}")
    print(f"match={ok}")
    print(f"steps={vm.steps}")
    print(f"rom_bytes={len(rom)}")
    if dump_bin:
        print(f"rom_hex={rom.hex()}")


def run_demo_soft64asm_cmp_u64(ua: int, ub: int, dump_bin: bool = False) -> None:
    rom = compile_firmware(Firmware.rows_u64_cmp())
    vm = VM16()
    vm.load_rom_bin(rom)
    _write_u64_le_words(vm, 0, ua)
    _write_u64_le_words(vm, 8, ub)
    vm.run()
    code = vm.mem_get_i16(16) & 0xFFFF
    py_code = 0 if (ua & 0xFFFFFFFFFFFFFFFF) < (ub & 0xFFFFFFFFFFFFFFFF) else (2 if (ua & 0xFFFFFFFFFFFFFFFF) > (ub & 0xFFFFFFFFFFFFFFFF) else 1)
    ok = code == py_code
    print(f"demo=soft64asm_u64cmp")
    print(f"ua=0x{ua & 0xFFFFFFFFFFFFFFFF:016X}")
    print(f"ub=0x{ub & 0xFFFFFFFFFFFFFFFF:016X}")
    print(f"vm_code={code} (0<,1==,2>)")
    print(f"py_code={py_code}")
    print(f"match={ok}")
    print(f"steps={vm.steps}")
    print(f"rom_bytes={len(rom)}")
    if dump_bin:
        print(f"rom_hex={rom.hex()}")


def run_demo_soft64asm_clz_u64(ua: int, dump_bin: bool = False) -> None:
    rom = compile_firmware(Firmware.rows_u64_clz())
    vm = VM16()
    vm.load_rom_bin(rom)
    _write_u64_le_words(vm, 0, ua)
    vm.run()
    out = vm.mem_get_i16(8) & 0xFFFF
    v = ua & 0xFFFFFFFFFFFFFFFF
    py = 64 if v == 0 else 64 - v.bit_length()
    ok = out == py
    print(f"demo=soft64asm_u64clz")
    print(f"ua=0x{v:016X}")
    print(f"vm={out}")
    print(f"py={py}")
    print(f"match={ok}")
    print(f"steps={vm.steps}")
    print(f"rom_bytes={len(rom)}")
    if dump_bin:
        print(f"rom_hex={rom.hex()}")


def run_demo_soft64asm_add_sameexp_pos(a: float, b: float, dump_bin: bool = False) -> None:
    a_bits = _f64_to_u64(a)
    b_bits = _f64_to_u64(b)
    out_bits, steps, rom_bytes = soft64_add_bits_vm(a_bits, b_bits)
    out = _u64_to_f64(out_bits)
    py = a + b
    _print_soft64_result("soft64asm_add", a, b, out, py, steps=steps, rom_bytes=rom_bytes, mode="vm_add_firmware")


def run_demo_soft64asm_fsub(a: float, b: float, dump_bin: bool = False) -> None:
    a_bits = _f64_to_u64(a)
    b_bits = _f64_to_u64(b)
    out_bits, steps, rom_bytes = soft64_sub_bits_vm(a_bits, b_bits)
    out = _u64_to_f64(out_bits)
    py = a - b
    _print_soft64_result("soft64asm_fsub", a, b, out, py, steps=steps, rom_bytes=rom_bytes, mode="vm_sub_via_add_firmware")


def run_demo_soft64asm_mul(a: float, b: float, dump_bin: bool = False) -> None:
    a_bits = _f64_to_u64(a)
    b_bits = _f64_to_u64(b)
    out_bits, steps, rom_bytes, mode = soft64_mul_bits_vm(a_bits, b_bits)
    out = _u64_to_f64(out_bits)
    py = a * b
    _print_soft64_result("soft64asm_mul", a, b, out, py, steps=steps, rom_bytes=rom_bytes, mode=mode)


def run_demo_soft64asm_div(a: float, b: float, dump_bin: bool = False) -> None:
    a_bits = _f64_to_u64(a)
    b_bits = _f64_to_u64(b)
    out_bits, steps, rom_bytes, mode = soft64_div_bits_vm(a_bits, b_bits)
    out = _u64_to_f64(out_bits)
    py = a / b
    _print_soft64_result("soft64asm_div", a, b, out, py, steps=steps, rom_bytes=rom_bytes, mode=mode)
