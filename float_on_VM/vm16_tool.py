#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from vm16.assembler import Assembler


def cmd_asm(input_path: Path, output_path: Path) -> None:
    text = input_path.read_text(encoding="utf-8")
    rows = Assembler.parse_asm_text(text)
    words = Assembler.assemble_rows(rows)
    data = Assembler.words_to_bin(words)
    output_path.write_bytes(data)
    print(f"assembled: {input_path} -> {output_path} ({len(data)} bytes)")


def cmd_disasm(input_path: Path, output_path: Path | None) -> None:
    data = input_path.read_bytes()
    words = Assembler.bin_to_words(data)
    lines = Assembler.disassemble_words(words)
    text = "\n".join(lines) + "\n"
    if output_path:
        output_path.write_text(text, encoding="utf-8")
        print(f"disassembled: {input_path} -> {output_path} ({len(words)} words)")
    else:
        print(text, end="")


def main() -> None:
    parser = argparse.ArgumentParser(description="vm16 asm/bin tool")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_asm = sub.add_parser("asm", help="assemble .asm to .bin")
    p_asm.add_argument("input", type=Path)
    p_asm.add_argument("output", type=Path)

    p_dis = sub.add_parser("disasm", help="disassemble .bin")
    p_dis.add_argument("input", type=Path)
    p_dis.add_argument("--out", type=Path, default=None)

    args = parser.parse_args()
    if args.cmd == "asm":
        cmd_asm(args.input, args.output)
    else:
        cmd_disasm(args.input, args.out)


if __name__ == "__main__":
    main()

