from __future__ import annotations

from typing import Dict, Tuple


OPCODE: Dict[str, int] = {
    "HALT": 0,
    "MOVI": 1,
    "MOV": 2,
    "LOAD": 3,
    "STORE": 4,
    "ADD": 5,
    "SUB": 6,
    "MUL": 7,
    "DIV": 8,
    "CMP": 9,
    "JMP": 10,
    "JZ": 11,
    "JNZ": 12,
    "JN": 13,
    "JP": 14,
    "PUSH": 15,
    "POP": 16,
    "CALL": 17,
    "RET": 18,
    "ADDC": 19,
    "SUBB": 20,
    "AND": 21,
    "OR": 22,
    "XOR": 23,
    "SHL": 24,
    "SHR": 25,
    "JC": 26,
    "JNC": 27,
    "SAR": 28,
    "UMUL": 29,
}
OPNAME = {v: k for k, v in OPCODE.items()}


class ISA:
    @staticmethod
    def _op(op: str) -> int:
        if op not in OPCODE:
            raise ValueError(f"unknown opcode: {op}")
        return OPCODE[op]

    @staticmethod
    def encode_rr(op: str, ra: int, rb: int) -> int:
        return (ISA._op(op) << 11) | ((ra & 0x7) << 8) | ((rb & 0x7) << 5)

    @staticmethod
    def encode_ri8(op: str, ra: int, imm8: int) -> int:
        return (ISA._op(op) << 11) | ((ra & 0x7) << 8) | (imm8 & 0xFF)

    @staticmethod
    def encode_j11(op: str, addr11: int) -> int:
        return (ISA._op(op) << 11) | (addr11 & 0x7FF)

    @staticmethod
    def encode_r(op: str, ra: int) -> int:
        return (ISA._op(op) << 11) | ((ra & 0x7) << 8)

    @staticmethod
    def encode_0(op: str) -> int:
        return ISA._op(op) << 11

    @staticmethod
    def decode(word: int) -> Tuple[str, int, int, int]:
        opv = (word >> 11) & 0x1F
        op = OPNAME.get(opv, f"OP{opv}")
        ra = (word >> 8) & 0x7
        rb = (word >> 5) & 0x7
        imm8 = word & 0xFF
        addr11 = word & 0x7FF
        return op, ra, rb, imm8 if op in {"MOVI", "LOAD", "STORE"} else addr11
