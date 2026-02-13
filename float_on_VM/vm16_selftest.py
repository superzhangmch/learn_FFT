#!/usr/bin/env python3

from __future__ import annotations

from vm16.demos import (
    run_demo_soft64asm_add_u64,
    run_demo_soft64asm_cmp_u64,
    run_demo_soft64asm_div128_u64,
    run_demo_soft64asm_mul128_u64,
    run_demo_soft64asm_mul_u64,
    run_demo_soft64asm_sub_u64,
)


def test_vm_u64_primitives() -> None:
    run_demo_soft64asm_add_u64(0x0123456789ABCDEF, 0x1111222233334444)
    run_demo_soft64asm_sub_u64(0x0, 0x1)
    run_demo_soft64asm_mul_u64(0x123456789ABCDEF0, 0x0FEDCBA987654321)
    run_demo_soft64asm_mul128_u64(0x123456789ABCDEF0, 0x0FEDCBA987654321)
    run_demo_soft64asm_div128_u64(0x0123456789ABCDEF, 0x1111222233334444, 0x1234)
    run_demo_soft64asm_cmp_u64(0xFFFFFFFFFFFFFFFF, 0x1)


def main() -> None:
    test_vm_u64_primitives()
    print("selftest: PASS")


if __name__ == "__main__":
    main()
