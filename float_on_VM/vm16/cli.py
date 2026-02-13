from __future__ import annotations

import argparse

from .demos import (
    run_demo_acosvm,
    run_demo_asinvm,
    run_demo_atanvm,
    run_demo_cosvm,
    run_demo_expvm,
    run_demo_expr,
    run_demo_logvm,
    run_demo_powvm,
    run_demo_sinvm,
    run_demo_sqrtvm,
    run_demo_tanvm,
    run_demo_soft64asm_add_u64,
    run_demo_soft64asm_clz_u64,
    run_demo_soft64asm_cmp_u64,
    run_demo_soft64asm_add_sameexp_pos,
    run_demo_soft64asm_fsub,
    run_demo_soft64asm_shr1_u64,
    run_demo_soft64asm_shl1_u64,
    run_demo_soft64asm_sub_u64,
    run_demo_soft64asm_mul_u64,
    run_demo_soft64asm_mul128_u64,
    run_demo_soft64asm_div128_u64,
    run_demo_soft64asm_mul,
    run_demo_soft64asm_div,
    run_demo_sum,
)

def _auto_int(v: str) -> int:
    return int(v, 0)


def main() -> None:
    parser = argparse.ArgumentParser(description="16-bit virtual hardware calculator (fp64 primitives)")
    parser.add_argument("--demo", choices=["sum", "expr", "sinvm", "cosvm", "tanvm", "atanvm", "asinvm", "acosvm", "sqrtvm", "logvm", "expvm", "powvm", "soft64asm", "soft64asm_sub", "soft64asm_mul64", "soft64asm_mul128", "soft64asm_div128", "soft64asm_shr1", "soft64asm_shl1", "soft64asm_cmp", "soft64asm_clz", "soft64asm_add", "soft64asm_fsub", "soft64asm_mul", "soft64asm_div"], default="sum")
    parser.add_argument("--n", type=int, default=10, help="for sum demo: compute 1+...+n")
    parser.add_argument("--a", type=int, default=12, help="for expr demo")
    parser.add_argument("--b", type=int, default=8, help="for expr demo")
    parser.add_argument("--c", type=int, default=5, help="for expr demo")
    parser.add_argument("--d", type=int, default=2, help="for expr demo")
    parser.add_argument("--x", type=float, default=0.5, help="for sin/cos/tan/atan/asin/acos/sqrt/log/exp demos")
    parser.add_argument("--deg", action="store_true", help="sin/cos/tan: input x as degrees")
    parser.add_argument("--iter", type=int, default=32, help="series/iteration budget for trig and inverse-trig demos")
    parser.add_argument("--fa", type=float, default=1.23456789012345, help="for soft64asm add/mul/div demos")
    parser.add_argument("--fb", type=float, default=9.87654321098765, help="for soft64asm add/mul/div demos")
    parser.add_argument("--ua", type=_auto_int, default=0x0123456789ABCDEF, help="for soft64asm demo (u64, supports 0x...)")
    parser.add_argument("--ub", type=_auto_int, default=0x1111222233334444, help="for soft64asm demo (u64, supports 0x...)")
    parser.add_argument("--un-hi", type=_auto_int, default=0x0123456789ABCDEF, help="for soft64asm_div128: numerator high 64")
    parser.add_argument("--un-lo", type=_auto_int, default=0x1111222233334444, help="for soft64asm_div128: numerator low 64")
    parser.add_argument("--ud", type=_auto_int, default=0x0000000000001234, help="for soft64asm_div128: denominator u64")
    parser.add_argument("--dump-bin", action="store_true", help="print ROM hex")
    args = parser.parse_args()

    if args.demo == "sum":
        run_demo_sum(args.n, dump_bin=args.dump_bin)
    elif args.demo == "expr":
        run_demo_expr(args.a, args.b, args.c, args.d, dump_bin=args.dump_bin)
    elif args.demo == "sinvm":
        run_demo_sinvm(args.x, deg=args.deg, dump_bin=args.dump_bin)
    elif args.demo == "cosvm":
        run_demo_cosvm(args.x, deg=args.deg, dump_bin=args.dump_bin)
    elif args.demo == "tanvm":
        run_demo_tanvm(args.x, deg=args.deg, iterations=args.iter, dump_bin=args.dump_bin)
    elif args.demo == "atanvm":
        run_demo_atanvm(args.x, iterations=args.iter, dump_bin=args.dump_bin)
    elif args.demo == "asinvm":
        run_demo_asinvm(args.x, iterations=args.iter, dump_bin=args.dump_bin)
    elif args.demo == "acosvm":
        run_demo_acosvm(args.x, iterations=args.iter, dump_bin=args.dump_bin)
    elif args.demo == "sqrtvm":
        run_demo_sqrtvm(args.x, dump_bin=args.dump_bin)
    elif args.demo == "logvm":
        run_demo_logvm(args.x, dump_bin=args.dump_bin)
    elif args.demo == "expvm":
        run_demo_expvm(args.x, dump_bin=args.dump_bin)
    elif args.demo == "powvm":
        run_demo_powvm(args.fa, args.fb, dump_bin=args.dump_bin)
    elif args.demo == "soft64asm":
        run_demo_soft64asm_add_u64(args.ua, args.ub, dump_bin=args.dump_bin)
    elif args.demo == "soft64asm_sub":
        run_demo_soft64asm_sub_u64(args.ua, args.ub, dump_bin=args.dump_bin)
    elif args.demo == "soft64asm_mul64":
        run_demo_soft64asm_mul_u64(args.ua, args.ub, dump_bin=args.dump_bin)
    elif args.demo == "soft64asm_mul128":
        run_demo_soft64asm_mul128_u64(args.ua, args.ub, dump_bin=args.dump_bin)
    elif args.demo == "soft64asm_div128":
        run_demo_soft64asm_div128_u64(args.un_hi, args.un_lo, args.ud, dump_bin=args.dump_bin)
    elif args.demo == "soft64asm_shr1":
        run_demo_soft64asm_shr1_u64(args.ua, dump_bin=args.dump_bin)
    elif args.demo == "soft64asm_shl1":
        run_demo_soft64asm_shl1_u64(args.ua, dump_bin=args.dump_bin)
    elif args.demo == "soft64asm_cmp":
        run_demo_soft64asm_cmp_u64(args.ua, args.ub, dump_bin=args.dump_bin)
    elif args.demo == "soft64asm_clz":
        run_demo_soft64asm_clz_u64(args.ua, dump_bin=args.dump_bin)
    elif args.demo == "soft64asm_add":
        run_demo_soft64asm_add_sameexp_pos(args.fa, args.fb, dump_bin=args.dump_bin)
    elif args.demo == "soft64asm_fsub":
        run_demo_soft64asm_fsub(args.fa, args.fb, dump_bin=args.dump_bin)
    elif args.demo == "soft64asm_mul":
        run_demo_soft64asm_mul(args.fa, args.fb, dump_bin=args.dump_bin)
    elif args.demo == "soft64asm_div":
        run_demo_soft64asm_div(args.fa, args.fb, dump_bin=args.dump_bin)


if __name__ == "__main__":
    main()
