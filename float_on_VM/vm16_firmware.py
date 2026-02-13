#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
from pathlib import Path

from vm16.assembler import Assembler
from vm16.core import VM16
from vm16.firmware import Firmware


def _i48_to_words(v: int) -> tuple[int, int, int]:
    u = v & ((1 << 48) - 1)
    return u & 0xFFFF, (u >> 16) & 0xFFFF, (u >> 32) & 0xFFFF


def _words_to_i48(w0: int, w1: int, w2: int) -> int:
    u = (w0 & 0xFFFF) | ((w1 & 0xFFFF) << 16) | ((w2 & 0xFFFF) << 32)
    if u & (1 << 47):
        u -= 1 << 48
    return u


def _write_i48(vm: VM16, base: int, v: int) -> None:
    w0, w1, w2 = _i48_to_words(v)
    vm.mem_set_i16(base, w0)
    vm.mem_set_i16(base + 2, w1)
    vm.mem_set_i16(base + 4, w2)


def _read_i48(vm: VM16, base: int) -> int:
    w0 = vm.mem_get_i16(base) & 0xFFFF
    w1 = vm.mem_get_i16(base + 2) & 0xFFFF
    w2 = vm.mem_get_i16(base + 4) & 0xFFFF
    return _words_to_i48(w0, w1, w2)


def _normalize_angle(theta: float) -> tuple[float, int]:
    pi = math.pi
    two_pi = 2.0 * pi
    half_pi = pi / 2.0
    t = (theta + pi) % two_pi - pi
    sign = 1
    if t > half_pi:
        t -= pi
        sign = -1
    elif t < -half_pi:
        t += pi
        sign = -1
    return t, sign


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

CORDIC_K_32 = 0.6072529350088814


def build_sin_firmware(asm_path: Path, bin_path: Path, iterations: int = 32) -> None:
    rows = Firmware.rows_sin_cordic48(iterations=iterations, angle_base=32)
    words = Assembler.assemble_rows(rows)
    asm_text = "\n".join(Assembler.disassemble_words(words)) + "\n"
    rom = Assembler.words_to_bin(words)
    asm_path.parent.mkdir(parents=True, exist_ok=True)
    bin_path.parent.mkdir(parents=True, exist_ok=True)
    asm_path.write_text(asm_text, encoding="utf-8")
    bin_path.write_bytes(rom)
    print(f"built asm={asm_path} words={len(words)}")
    print(f"built bin={bin_path} bytes={len(rom)}")


def build_atan_firmware(asm_path: Path, bin_path: Path, iterations: int = 32) -> None:
    rows = Firmware.rows_atan_cordic48(iterations=iterations, angle_base=32)
    words = Assembler.assemble_rows(rows)
    asm_text = "\n".join(Assembler.disassemble_words(words)) + "\n"
    rom = Assembler.words_to_bin(words)
    asm_path.parent.mkdir(parents=True, exist_ok=True)
    bin_path.parent.mkdir(parents=True, exist_ok=True)
    asm_path.write_text(asm_text, encoding="utf-8")
    bin_path.write_bytes(rom)
    print(f"built asm={asm_path} words={len(words)}")
    print(f"built bin={bin_path} bytes={len(rom)}")


def build_sqrt_firmware(asm_path: Path, bin_path: Path, iterations: int = 16) -> None:
    rows = Firmware.rows_sqrt_q46_newton(iterations=iterations)
    words = Assembler.assemble_rows(rows)
    asm_text = "\n".join(Assembler.disassemble_words(words)) + "\n"
    rom = Assembler.words_to_bin(words)
    asm_path.parent.mkdir(parents=True, exist_ok=True)
    bin_path.parent.mkdir(parents=True, exist_ok=True)
    asm_path.write_text(asm_text, encoding="utf-8")
    bin_path.write_bytes(rom)
    print(f"built asm={asm_path} words={len(words)}")
    print(f"built bin={bin_path} bytes={len(rom)}")


def run_sin_from_bin(bin_path: Path, x: float, deg: bool = False) -> None:
    rom = bin_path.read_bytes()
    vm = VM16()
    vm.load_rom_bin(rom)

    theta = math.radians(x) if deg else x
    t_norm, sign = _normalize_angle(theta)
    scale = 1 << 46
    _write_i48(vm, 0, int(round(CORDIC_K_32 * scale)))
    _write_i48(vm, 6, 0)
    _write_i48(vm, 12, int(round(t_norm * scale)))
    for i in range(32):
        _write_i48(vm, 32 + i * 6, int(round(CORDIC_ANGLES[i] * scale)))

    vm.run(max_steps=3_000_000)
    y_q46 = _read_i48(vm, 6)
    s = sign * (y_q46 / float(scale))
    py = math.sin(theta)
    print("demo=sinvm_bin")
    print(f"fw_bin={bin_path}")
    print(f"x={'deg' if deg else 'rad'}={x:.16e}")
    print(f"vm={s:.16e}")
    print(f"py={py:.16e}")
    print(f"abs_err={abs(s - py):.3e}")
    print(f"steps={vm.steps}")
    print(f"rom_bytes={len(rom)}")


def run_atan_from_bin(bin_path: Path, x: float) -> None:
    rom = bin_path.read_bytes()
    vm = VM16()
    vm.load_rom_bin(rom)

    scale = 1 << 46
    vm_x = x
    adjust = 0.0
    if abs(x) > 1.0:
        vm_x = 1.0 / x
        adjust = math.pi / 2.0 if x > 0 else -math.pi / 2.0

    _write_i48(vm, 0, scale)
    _write_i48(vm, 6, int(round(vm_x * scale)))
    _write_i48(vm, 12, 0)
    for i in range(32):
        _write_i48(vm, 32 + i * 6, int(round(CORDIC_ANGLES[i] * scale)))

    vm.run(max_steps=3_000_000)
    z_q46 = _read_i48(vm, 12)
    vm_out = z_q46 / float(scale)
    out = adjust - vm_out if abs(x) > 1.0 else vm_out
    py = math.atan(x)
    print("demo=atanvm_bin")
    print(f"fw_bin={bin_path}")
    print(f"x={x:.16e}")
    print(f"vm={out:.16e}")
    print(f"py={py:.16e}")
    print(f"abs_err={abs(out - py):.3e}")
    print(f"steps={vm.steps}")
    print(f"rom_bytes={len(rom)}")


def run_sqrt_from_bin(bin_path: Path, x: float) -> None:
    if x < 0:
        raise ValueError("sqrt domain is x>=0")
    rom = bin_path.read_bytes()
    vm = VM16()
    vm.load_rom_bin(rom)
    scale = 1 << 46
    x_q = int(round(x * scale))
    for i in range(4):
        vm.mem_set_i16(100 + i * 2, (x_q >> (16 * i)) & 0xFFFF)
    vm.run(max_steps=3_000_000)
    out_q = 0
    for i in range(4):
        out_q |= (vm.mem_get_i16(8 + i * 2) & 0xFFFF) << (16 * i)
    out = out_q / float(scale)
    py = math.sqrt(x)
    print("demo=sqrtvm_bin")
    print(f"fw_bin={bin_path}")
    print(f"x={x:.16e}")
    print(f"vm={out:.16e}")
    print(f"py={py:.16e}")
    print(f"abs_err={abs(out - py):.3e}")
    print(f"steps={vm.steps}")
    print(f"rom_bytes={len(rom)}")


def main() -> None:
    parser = argparse.ArgumentParser(description="vm16 firmware build/run tool")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_build = sub.add_parser("build-sin", help="generate sin_cordic48 asm+bin firmware")
    p_build.add_argument("--asm", type=Path, default=Path("firmware/sin_cordic48.asm"))
    p_build.add_argument("--bin", type=Path, default=Path("firmware/sin_cordic48.bin"))
    p_build.add_argument("--iter", type=int, default=32)

    p_run = sub.add_parser("run-sin", help="run sin firmware .bin in VM")
    p_run.add_argument("x", type=float)
    p_run.add_argument("--deg", action="store_true")
    p_run.add_argument("--bin", type=Path, default=Path("firmware/sin_cordic48.bin"))

    p_build_atan = sub.add_parser("build-atan", help="generate atan_cordic48 asm+bin firmware")
    p_build_atan.add_argument("--asm", type=Path, default=Path("firmware/atan_cordic48.asm"))
    p_build_atan.add_argument("--bin", type=Path, default=Path("firmware/atan_cordic48.bin"))
    p_build_atan.add_argument("--iter", type=int, default=32)

    p_run_atan = sub.add_parser("run-atan", help="run atan firmware .bin in VM")
    p_run_atan.add_argument("x", type=float)
    p_run_atan.add_argument("--bin", type=Path, default=Path("firmware/atan_cordic48.bin"))

    p_build_sqrt = sub.add_parser("build-sqrt", help="generate sqrt_q46_newton asm+bin firmware")
    p_build_sqrt.add_argument("--asm", type=Path, default=Path("firmware/sqrt_q46_newton.asm"))
    p_build_sqrt.add_argument("--bin", type=Path, default=Path("firmware/sqrt_q46_newton.bin"))
    p_build_sqrt.add_argument("--iter", type=int, default=16)

    p_run_sqrt = sub.add_parser("run-sqrt", help="run sqrt firmware .bin in VM")
    p_run_sqrt.add_argument("x", type=float)
    p_run_sqrt.add_argument("--bin", type=Path, default=Path("firmware/sqrt_q46_newton.bin"))

    args = parser.parse_args()
    if args.cmd == "build-sin":
        build_sin_firmware(args.asm, args.bin, iterations=args.iter)
    elif args.cmd == "run-sin":
        run_sin_from_bin(args.bin, args.x, deg=args.deg)
    elif args.cmd == "build-atan":
        build_atan_firmware(args.asm, args.bin, iterations=args.iter)
    elif args.cmd == "run-atan":
        run_atan_from_bin(args.bin, args.x)
    elif args.cmd == "build-sqrt":
        build_sqrt_firmware(args.asm, args.bin, iterations=args.iter)
    else:
        run_sqrt_from_bin(args.bin, args.x)


if __name__ == "__main__":
    main()
