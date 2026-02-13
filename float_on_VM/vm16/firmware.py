from __future__ import annotations

from typing import List

from .assembler import AsmRow
from .core import VM16


class Firmware:
    @staticmethod
    def _append_prefixed_subroutine(rows: List[AsmRow], base_rows: List[AsmRow], prefix: str, entry_label: str) -> None:
        label_names = {str(r[1]) for r in base_rows if r and r[0] == "LABEL"}
        branch_ops = {"JMP", "JZ", "JNZ", "JN", "JP", "JC", "JNC", "CALL"}
        rows.append(("LABEL", entry_label))
        for row in base_rows:
            op = str(row[0])
            if op == "HALT":
                rows.append(("RET",))
                continue
            if op == "LABEL":
                rows.append(("LABEL", f"{prefix}{row[1]}"))
                continue
            if op in branch_ops and len(row) > 1 and isinstance(row[1], str) and row[1] in label_names:
                if len(row) == 2:
                    rows.append((op, f"{prefix}{row[1]}"))
                else:
                    rows.append((op, f"{prefix}{row[1]}", row[2]))
                continue
            rows.append(row)

    @staticmethod
    def rows_sum_1_to_n_call() -> List[AsmRow]:
        return [
            ("MOVI", 0, 0),
            ("MOVI", 1, 1),
            ("MOVI", 3, 1),
            ("LOAD", 2, 0),
            ("LABEL", "loop"),
            ("CMP", 1, 2),
            ("JP", "done"),
            ("CALL", "add_i"),
            ("CALL", "inc_i"),
            ("JMP", "loop"),
            ("LABEL", "done"),
            ("STORE", 0, 2),
            ("HALT",),
            ("LABEL", "add_i"),
            ("ADD", 0, 1),
            ("RET",),
            ("LABEL", "inc_i"),
            ("ADD", 1, 3),
            ("RET",),
        ]

    @staticmethod
    def rows_expr_abcd() -> List[AsmRow]:
        return [
            ("LOAD", 0, 0),
            ("LOAD", 1, 2),
            ("LOAD", 2, 4),
            ("LOAD", 3, 6),
            ("ADD", 0, 1),
            ("MUL", 0, 2),
            ("DIV", 0, 3),
            ("STORE", 0, 8),
            ("HALT",),
        ]

    @staticmethod
    def rows_u64_add() -> List[AsmRow]:
        # little-endian 16-bit limbs
        # A at [0..7], B at [8..15], OUT at [16..23]
        return [
            ("LOAD", 0, 0),
            ("LOAD", 1, 8),
            ("ADD", 0, 1),
            ("STORE", 0, 16),

            ("LOAD", 0, 2),
            ("LOAD", 1, 10),
            ("ADDC", 0, 1),
            ("STORE", 0, 18),

            ("LOAD", 0, 4),
            ("LOAD", 1, 12),
            ("ADDC", 0, 1),
            ("STORE", 0, 20),

            ("LOAD", 0, 6),
            ("LOAD", 1, 14),
            ("ADDC", 0, 1),
            ("STORE", 0, 22),
            ("HALT",),
        ]

    @staticmethod
    def rows_u64_sub() -> List[AsmRow]:
        # OUT = A - B (mod 2^64), little-endian 16-bit limbs
        # A at [0..7], B at [8..15], OUT at [16..23]
        return [
            ("LOAD", 0, 0),
            ("LOAD", 1, 8),
            ("SUB", 0, 1),
            ("STORE", 0, 16),

            ("LOAD", 0, 2),
            ("LOAD", 1, 10),
            ("SUBB", 0, 1),
            ("STORE", 0, 18),

            ("LOAD", 0, 4),
            ("LOAD", 1, 12),
            ("SUBB", 0, 1),
            ("STORE", 0, 20),

            ("LOAD", 0, 6),
            ("LOAD", 1, 14),
            ("SUBB", 0, 1),
            ("STORE", 0, 22),
            ("HALT",),
        ]

    @staticmethod
    def rows_u64_mul_low() -> List[AsmRow]:
        # OUT = (A * B) mod 2^64, unsigned.
        # A at [0..7], B at [8..15], OUT at [16..23] (little-endian 16-bit limbs).
        rows: List[AsmRow] = [
            ("MOVI", 6, 1),
            ("MOVI", 5, 0),
            ("STORE", 5, 16),
            ("STORE", 5, 18),
            ("STORE", 5, 20),
            ("STORE", 5, 22),
        ]

        for i in range(4):
            ai = 2 * i
            for j in range(4):
                k = i + j
                if k >= 4:
                    continue
                bj = 8 + 2 * j
                out0 = 16 + 2 * k
                base = f"mul_{i}_{j}"

                rows.extend([
                    ("LOAD", 0, ai),
                    ("LOAD", 1, bj),
                    ("UMUL", 0, 1),  # r0=lo16, r1=hi16

                    ("LOAD", 2, out0),
                    ("ADD", 2, 0),
                    ("STORE", 2, out0),
                ])

                if k + 1 < 4:
                    out1 = 16 + 2 * (k + 1)
                    rows.extend([
                        ("LOAD", 3, out1),
                        ("ADDC", 3, 1),
                        ("STORE", 3, out1),
                    ])

                    for t in range(k + 2, 4):
                        outt = 16 + 2 * t
                        noc = f"{base}_noc_{t}"
                        rows.extend([
                            ("JNC", noc),
                            ("LOAD", 4, outt),
                            ("ADD", 4, 6),
                            ("STORE", 4, outt),
                            ("LABEL", noc),
                        ])

        rows.append(("HALT",))
        return rows

    @staticmethod
    def rows_u64_mul_full() -> List[AsmRow]:
        # OUT128 = A * B (unsigned), full 128-bit result.
        # A at [0..7], B at [8..15], OUT at [16..31], little-endian 16-bit limbs.
        rows: List[AsmRow] = [
            ("MOVI", 6, 1),   # one
            ("MOVI", 5, 0),   # zero
            ("STORE", 5, 16),
            ("STORE", 5, 18),
            ("STORE", 5, 20),
            ("STORE", 5, 22),
            ("STORE", 5, 24),
            ("STORE", 5, 26),
            ("STORE", 5, 28),
            ("STORE", 5, 30),
        ]

        for i in range(4):
            ai = 2 * i
            for j in range(4):
                bj = 8 + 2 * j
                k = i + j
                out0 = 16 + 2 * k
                out1 = out0 + 2
                base = f"mulf_{i}_{j}"

                rows.extend([
                    ("LOAD", 0, ai),
                    ("LOAD", 1, bj),
                    ("UMUL", 0, 1),      # r0=lo16, r1=hi16

                    ("LOAD", 2, out0),
                    ("ADD", 2, 0),
                    ("STORE", 2, out0),

                    ("LOAD", 3, out1),
                    ("ADDC", 3, 1),
                    ("STORE", 3, out1),
                ])

                for t in range(k + 2, 8):
                    outt = 16 + 2 * t
                    noc = f"{base}_noc_{t}"
                    rows.extend([
                        ("JNC", noc),
                        ("LOAD", 4, outt),
                        ("ADD", 4, 6),
                        ("STORE", 4, outt),
                        ("LABEL", noc),
                    ])

        rows.append(("HALT",))
        return rows

    @staticmethod
    def rows_u128_div_u64() -> List[AsmRow]:
        # Unsigned long division core:
        # Q = floor(N / D), R = N % D
        #
        # N (128-bit) at [0..15]  (8 words)
        # D (64-bit)  at [16..23] (4 words), D!=0
        # Q (128-bit) at [24..39] (8 words)
        # R (64-bit)  at [40..47] (4 words)
        # NCUR copy   at [48..63] (8 words)
        # CNT at [64], CMP code at [66], bit tmp at [68]
        return [
            ("MOVI", 6, 1),   # one
            ("MOVI", 5, 0),   # zero
            # r4 = 0x8000
            ("MOVI", 4, 128),
            ("SHL", 4, 6), ("SHL", 4, 6), ("SHL", 4, 6), ("SHL", 4, 6),
            ("SHL", 4, 6), ("SHL", 4, 6), ("SHL", 4, 6), ("SHL", 4, 6),

            # NCUR = N
            ("LOAD", 0, 0), ("STORE", 0, 48),
            ("LOAD", 0, 2), ("STORE", 0, 50),
            ("LOAD", 0, 4), ("STORE", 0, 52),
            ("LOAD", 0, 6), ("STORE", 0, 54),
            ("LOAD", 0, 8), ("STORE", 0, 56),
            ("LOAD", 0, 10), ("STORE", 0, 58),
            ("LOAD", 0, 12), ("STORE", 0, 60),
            ("LOAD", 0, 14), ("STORE", 0, 62),

            # Q=0, R=0
            ("STORE", 5, 24), ("STORE", 5, 26), ("STORE", 5, 28), ("STORE", 5, 30),
            ("STORE", 5, 32), ("STORE", 5, 34), ("STORE", 5, 36), ("STORE", 5, 38),
            ("STORE", 5, 40), ("STORE", 5, 42), ("STORE", 5, 44), ("STORE", 5, 46),

            ("MOVI", 0, 128), ("STORE", 0, 64),

            ("LABEL", "d_loop"),
            ("LOAD", 0, 64), ("CMP", 0, 5), ("JZ", "d_done"),

            # bit = msb(NCUR)
            ("LOAD", 0, 62), ("AND", 0, 4), ("CMP", 0, 5), ("JZ", "d_bit0"),
            ("MOVI", 2, 1), ("JMP", "d_bit_done"),
            ("LABEL", "d_bit0"),
            ("MOVI", 2, 0),
            ("LABEL", "d_bit_done"),

            ("STORE", 2, 68),
            ("CALL", "d_shl1_ncur"),
            ("CALL", "d_shl1_r"),
            ("LOAD", 2, 68),
            ("CMP", 2, 5), ("JZ", "d_no_r_bit"), ("LOAD", 0, 40), ("OR", 0, 6), ("STORE", 0, 40),
            ("LABEL", "d_no_r_bit"),
            ("CALL", "d_shl1_q"),

            ("CALL", "d_cmp_r_d"),
            ("LOAD", 0, 66), ("CMP", 0, 5), ("JZ", "d_no_sub"),
            ("CALL", "d_sub_r_d"),
            ("LOAD", 0, 24), ("OR", 0, 6), ("STORE", 0, 24),
            ("LABEL", "d_no_sub"),

            ("LOAD", 0, 64), ("SUB", 0, 6), ("STORE", 0, 64),
            ("JMP", "d_loop"),

            ("LABEL", "d_done"),
            ("HALT",),

            # NCUR <<= 1 (8 words at 48..62)
            ("LABEL", "d_shl1_ncur"),
            ("LOAD", 0, 48), ("SHL", 0, 6), ("STORE", 0, 48), ("JC", "dn_c01"), ("MOVI", 2, 0), ("JMP", "dn_d01"), ("LABEL", "dn_c01"), ("MOVI", 2, 1), ("LABEL", "dn_d01"),
            ("LOAD", 0, 50), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "dn_o1"), ("OR", 0, 6), ("LABEL", "dn_o1"), ("STORE", 0, 50), ("JC", "dn_c12"), ("MOVI", 2, 0), ("JMP", "dn_d12"), ("LABEL", "dn_c12"), ("MOVI", 2, 1), ("LABEL", "dn_d12"),
            ("LOAD", 0, 52), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "dn_o2"), ("OR", 0, 6), ("LABEL", "dn_o2"), ("STORE", 0, 52), ("JC", "dn_c23"), ("MOVI", 2, 0), ("JMP", "dn_d23"), ("LABEL", "dn_c23"), ("MOVI", 2, 1), ("LABEL", "dn_d23"),
            ("LOAD", 0, 54), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "dn_o3"), ("OR", 0, 6), ("LABEL", "dn_o3"), ("STORE", 0, 54), ("JC", "dn_c34"), ("MOVI", 2, 0), ("JMP", "dn_d34"), ("LABEL", "dn_c34"), ("MOVI", 2, 1), ("LABEL", "dn_d34"),
            ("LOAD", 0, 56), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "dn_o4"), ("OR", 0, 6), ("LABEL", "dn_o4"), ("STORE", 0, 56), ("JC", "dn_c45"), ("MOVI", 2, 0), ("JMP", "dn_d45"), ("LABEL", "dn_c45"), ("MOVI", 2, 1), ("LABEL", "dn_d45"),
            ("LOAD", 0, 58), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "dn_o5"), ("OR", 0, 6), ("LABEL", "dn_o5"), ("STORE", 0, 58), ("JC", "dn_c56"), ("MOVI", 2, 0), ("JMP", "dn_d56"), ("LABEL", "dn_c56"), ("MOVI", 2, 1), ("LABEL", "dn_d56"),
            ("LOAD", 0, 60), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "dn_o6"), ("OR", 0, 6), ("LABEL", "dn_o6"), ("STORE", 0, 60), ("JC", "dn_c67"), ("MOVI", 2, 0), ("JMP", "dn_d67"), ("LABEL", "dn_c67"), ("MOVI", 2, 1), ("LABEL", "dn_d67"),
            ("LOAD", 0, 62), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "dn_o7"), ("OR", 0, 6), ("LABEL", "dn_o7"), ("STORE", 0, 62),
            ("RET",),

            # R <<= 1 (4 words at 40..46)
            ("LABEL", "d_shl1_r"),
            ("LOAD", 0, 40), ("SHL", 0, 6), ("STORE", 0, 40), ("JC", "dr_c01"), ("MOVI", 2, 0), ("JMP", "dr_d01"), ("LABEL", "dr_c01"), ("MOVI", 2, 1), ("LABEL", "dr_d01"),
            ("LOAD", 0, 42), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "dr_o1"), ("OR", 0, 6), ("LABEL", "dr_o1"), ("STORE", 0, 42), ("JC", "dr_c12"), ("MOVI", 2, 0), ("JMP", "dr_d12"), ("LABEL", "dr_c12"), ("MOVI", 2, 1), ("LABEL", "dr_d12"),
            ("LOAD", 0, 44), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "dr_o2"), ("OR", 0, 6), ("LABEL", "dr_o2"), ("STORE", 0, 44), ("JC", "dr_c23"), ("MOVI", 2, 0), ("JMP", "dr_d23"), ("LABEL", "dr_c23"), ("MOVI", 2, 1), ("LABEL", "dr_d23"),
            ("LOAD", 0, 46), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "dr_o3"), ("OR", 0, 6), ("LABEL", "dr_o3"), ("STORE", 0, 46),
            ("RET",),

            # Q <<= 1 (8 words at 24..38)
            ("LABEL", "d_shl1_q"),
            ("LOAD", 0, 24), ("SHL", 0, 6), ("STORE", 0, 24), ("JC", "dq_c01"), ("MOVI", 2, 0), ("JMP", "dq_d01"), ("LABEL", "dq_c01"), ("MOVI", 2, 1), ("LABEL", "dq_d01"),
            ("LOAD", 0, 26), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "dq_o1"), ("OR", 0, 6), ("LABEL", "dq_o1"), ("STORE", 0, 26), ("JC", "dq_c12"), ("MOVI", 2, 0), ("JMP", "dq_d12"), ("LABEL", "dq_c12"), ("MOVI", 2, 1), ("LABEL", "dq_d12"),
            ("LOAD", 0, 28), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "dq_o2"), ("OR", 0, 6), ("LABEL", "dq_o2"), ("STORE", 0, 28), ("JC", "dq_c23"), ("MOVI", 2, 0), ("JMP", "dq_d23"), ("LABEL", "dq_c23"), ("MOVI", 2, 1), ("LABEL", "dq_d23"),
            ("LOAD", 0, 30), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "dq_o3"), ("OR", 0, 6), ("LABEL", "dq_o3"), ("STORE", 0, 30), ("JC", "dq_c34"), ("MOVI", 2, 0), ("JMP", "dq_d34"), ("LABEL", "dq_c34"), ("MOVI", 2, 1), ("LABEL", "dq_d34"),
            ("LOAD", 0, 32), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "dq_o4"), ("OR", 0, 6), ("LABEL", "dq_o4"), ("STORE", 0, 32), ("JC", "dq_c45"), ("MOVI", 2, 0), ("JMP", "dq_d45"), ("LABEL", "dq_c45"), ("MOVI", 2, 1), ("LABEL", "dq_d45"),
            ("LOAD", 0, 34), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "dq_o5"), ("OR", 0, 6), ("LABEL", "dq_o5"), ("STORE", 0, 34), ("JC", "dq_c56"), ("MOVI", 2, 0), ("JMP", "dq_d56"), ("LABEL", "dq_c56"), ("MOVI", 2, 1), ("LABEL", "dq_d56"),
            ("LOAD", 0, 36), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "dq_o6"), ("OR", 0, 6), ("LABEL", "dq_o6"), ("STORE", 0, 36), ("JC", "dq_c67"), ("MOVI", 2, 0), ("JMP", "dq_d67"), ("LABEL", "dq_c67"), ("MOVI", 2, 1), ("LABEL", "dq_d67"),
            ("LOAD", 0, 38), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "dq_o7"), ("OR", 0, 6), ("LABEL", "dq_o7"), ("STORE", 0, 38),
            ("RET",),

            # compare R vs D unsigned -> [66]: 0 less, 1 eq, 2 greater
            ("LABEL", "d_cmp_r_d"),
            ("LOAD", 0, 46), ("LOAD", 1, 22), ("MOV", 2, 0), ("SUB", 2, 1), ("JNC", "d_less"), ("CMP", 2, 5), ("JNZ", "d_greater"),
            ("LOAD", 0, 44), ("LOAD", 1, 20), ("MOV", 2, 0), ("SUB", 2, 1), ("JNC", "d_less"), ("CMP", 2, 5), ("JNZ", "d_greater"),
            ("LOAD", 0, 42), ("LOAD", 1, 18), ("MOV", 2, 0), ("SUB", 2, 1), ("JNC", "d_less"), ("CMP", 2, 5), ("JNZ", "d_greater"),
            ("LOAD", 0, 40), ("LOAD", 1, 16), ("MOV", 2, 0), ("SUB", 2, 1), ("JNC", "d_less"), ("CMP", 2, 5), ("JNZ", "d_greater"),
            ("STORE", 6, 66), ("RET",),
            ("LABEL", "d_less"),
            ("STORE", 5, 66), ("RET",),
            ("LABEL", "d_greater"),
            ("MOVI", 0, 2), ("STORE", 0, 66), ("RET",),

            # R = R - D
            ("LABEL", "d_sub_r_d"),
            ("LOAD", 0, 40), ("LOAD", 1, 16), ("SUB", 0, 1), ("STORE", 0, 40),
            ("LOAD", 0, 42), ("LOAD", 1, 18), ("SUBB", 0, 1), ("STORE", 0, 42),
            ("LOAD", 0, 44), ("LOAD", 1, 20), ("SUBB", 0, 1), ("STORE", 0, 44),
            ("LOAD", 0, 46), ("LOAD", 1, 22), ("SUBB", 0, 1), ("STORE", 0, 46),
            ("RET",),
        ]

    @staticmethod
    def rows_u64_shr1() -> List[AsmRow]:
        # OUT = A >> 1, little-endian 16-bit limbs
        # A at [0..7], OUT at [8..15]
        # Uses R6=1, R5=0x8000
        return [
            ("MOVI", 6, 1),
            ("MOVI", 5, 128),     # 0x0080
            ("SHL", 5, 6),        # 0x0100
            ("SHL", 5, 6),        # 0x0200
            ("SHL", 5, 6),        # 0x0400
            ("SHL", 5, 6),        # 0x0800
            ("SHL", 5, 6),        # 0x1000
            ("SHL", 5, 6),        # 0x2000
            ("SHL", 5, 6),        # 0x4000
            ("SHL", 5, 6),        # 0x8000

            ("LOAD", 0, 6),       # w3
            ("SHR", 0, 6),
            ("STORE", 0, 14),

            ("LOAD", 0, 4),       # w2
            ("SHR", 0, 6),
            ("LOAD", 1, 6),       # old w3
            ("AND", 1, 6),        # keep lsb(old w3)
            ("SHL", 1, 6), ("SHL", 1, 6), ("SHL", 1, 6), ("SHL", 1, 6),
            ("SHL", 1, 6), ("SHL", 1, 6), ("SHL", 1, 6), ("SHL", 1, 6),
            ("SHL", 1, 6), ("SHL", 1, 6), ("SHL", 1, 6), ("SHL", 1, 6),
            ("SHL", 1, 6), ("SHL", 1, 6), ("SHL", 1, 6),
            ("OR", 0, 1),
            ("STORE", 0, 12),

            ("LOAD", 0, 2),       # w1
            ("SHR", 0, 6),
            ("LOAD", 1, 4),       # old w2
            ("AND", 1, 6),
            ("SHL", 1, 6), ("SHL", 1, 6), ("SHL", 1, 6), ("SHL", 1, 6),
            ("SHL", 1, 6), ("SHL", 1, 6), ("SHL", 1, 6), ("SHL", 1, 6),
            ("SHL", 1, 6), ("SHL", 1, 6), ("SHL", 1, 6), ("SHL", 1, 6),
            ("SHL", 1, 6), ("SHL", 1, 6), ("SHL", 1, 6),
            ("OR", 0, 1),
            ("STORE", 0, 10),

            ("LOAD", 0, 0),       # w0
            ("SHR", 0, 6),
            ("LOAD", 1, 2),       # old w1
            ("AND", 1, 6),
            ("SHL", 1, 6), ("SHL", 1, 6), ("SHL", 1, 6), ("SHL", 1, 6),
            ("SHL", 1, 6), ("SHL", 1, 6), ("SHL", 1, 6), ("SHL", 1, 6),
            ("SHL", 1, 6), ("SHL", 1, 6), ("SHL", 1, 6), ("SHL", 1, 6),
            ("SHL", 1, 6), ("SHL", 1, 6), ("SHL", 1, 6),
            ("OR", 0, 1),
            ("STORE", 0, 8),
            ("HALT",),
        ]

    @staticmethod
    def rows_sqrt_q46_newton(iterations: int = 16) -> List[AsmRow]:
        # Q46 sqrt firmware (single program, control flow in VM).
        # Input:  x_q46 u64 at [100..107]
        # Output: sqrt(x)_q46 u64 at [8..15]
        # Work:
        #   g at [108..115], iter at [116], one_q46 at [120..127]
        #   divider scratch follows rows_u128_div_u64 layout:
        #     N [0..15], D [16..23], Q [24..39], R [40..47], NCUR [48..63], CNT [64], CMP [66], bit [68]
        if iterations < 1:
            raise ValueError("iterations must be >= 1")

        rows: List[AsmRow] = [
            ("MOVI", 6, 1),   # one
            ("MOVI", 5, 0),   # zero
            # r3 = 0x8000
            ("MOVI", 3, 128),
            ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6),
            ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6),

            # one_q46 = 1<<46 at [120..127] => words [0,0,0x4000,0]
            ("STORE", 5, 120), ("STORE", 5, 122),
            ("MOVI", 0, 64),
            ("SHL", 0, 6), ("SHL", 0, 6), ("SHL", 0, 6), ("SHL", 0, 6),
            ("SHL", 0, 6), ("SHL", 0, 6), ("SHL", 0, 6), ("SHL", 0, 6),
            ("STORE", 0, 124),
            ("STORE", 5, 126),

            # x==0 => out=0
            ("LOAD", 0, 100), ("LOAD", 1, 102), ("OR", 0, 1), ("LOAD", 1, 104), ("OR", 0, 1), ("LOAD", 1, 106), ("OR", 0, 1),
            ("CMP", 0, 5), ("JNZ", "sq_nonzero"),
            ("STORE", 5, 8), ("STORE", 5, 10), ("STORE", 5, 12), ("STORE", 5, 14),
            ("HALT",),

            ("LABEL", "sq_nonzero"),
            # g = x
            ("LOAD", 0, 100), ("STORE", 0, 108),
            ("LOAD", 0, 102), ("STORE", 0, 110),
            ("LOAD", 0, 104), ("STORE", 0, 112),
            ("LOAD", 0, 106), ("STORE", 0, 114),

            # if x < 1.0 => g = 1.0
            ("LOAD", 0, 106), ("CMP", 0, 5), ("JN", "sq_g_set_one"), ("JP", "sq_g_ok"),
            ("LOAD", 0, 104), ("LOAD", 1, 124), ("CMP", 0, 1), ("JN", "sq_g_set_one"), ("JP", "sq_g_ok"),
            ("LOAD", 0, 102), ("CMP", 0, 5), ("JN", "sq_g_set_one"), ("JP", "sq_g_ok"),
            ("LOAD", 0, 100), ("CMP", 0, 5), ("JN", "sq_g_set_one"),
            ("JMP", "sq_g_ok"),
            ("LABEL", "sq_g_set_one"),
            ("LOAD", 0, 120), ("STORE", 0, 108),
            ("LOAD", 0, 122), ("STORE", 0, 110),
            ("LOAD", 0, 124), ("STORE", 0, 112),
            ("LOAD", 0, 126), ("STORE", 0, 114),
            ("LABEL", "sq_g_ok"),

            ("MOVI", 0, iterations), ("STORE", 0, 116),

            ("LABEL", "sq_loop"),
            ("LOAD", 0, 116), ("CMP", 0, 5), ("JZ", "sq_done"),

            # N = x
            ("LOAD", 0, 100), ("STORE", 0, 0),
            ("LOAD", 0, 102), ("STORE", 0, 2),
            ("LOAD", 0, 104), ("STORE", 0, 4),
            ("LOAD", 0, 106), ("STORE", 0, 6),
            ("STORE", 5, 8), ("STORE", 5, 10), ("STORE", 5, 12), ("STORE", 5, 14),

            # N <<= 46
            ("MOVI", 0, 46), ("STORE", 0, 130),
            ("LABEL", "sq_nshl_loop"),
            ("LOAD", 0, 130), ("CMP", 0, 5), ("JZ", "sq_nshl_done"),
            ("CALL", "sq_shl1_n"),
            ("LOAD", 0, 130), ("SUB", 0, 6), ("STORE", 0, 130),
            ("JMP", "sq_nshl_loop"),
            ("LABEL", "sq_nshl_done"),

            # D = g
            ("LOAD", 0, 108), ("STORE", 0, 16),
            ("LOAD", 0, 110), ("STORE", 0, 18),
            ("LOAD", 0, 112), ("STORE", 0, 20),
            ("LOAD", 0, 114), ("STORE", 0, 22),

            # Q = (x<<46)/g
            ("CALL", "sq_div_entry"),

            # g = (g + q_low64) >> 1
            ("LOAD", 0, 108), ("LOAD", 1, 24), ("ADD", 0, 1), ("STORE", 0, 108),
            ("LOAD", 0, 110), ("LOAD", 1, 26), ("ADDC", 0, 1), ("STORE", 0, 110),
            ("LOAD", 0, 112), ("LOAD", 1, 28), ("ADDC", 0, 1), ("STORE", 0, 112),
            ("LOAD", 0, 114), ("LOAD", 1, 30), ("ADDC", 0, 1), ("STORE", 0, 114),
            ("CALL", "sq_shr1_g"),

            ("LOAD", 0, 116), ("SUB", 0, 6), ("STORE", 0, 116),
            ("JMP", "sq_loop"),

            ("LABEL", "sq_done"),
            ("LOAD", 0, 108), ("STORE", 0, 8),
            ("LOAD", 0, 110), ("STORE", 0, 10),
            ("LOAD", 0, 112), ("STORE", 0, 12),
            ("LOAD", 0, 114), ("STORE", 0, 14),
            ("HALT",),

            # N(0..14) <<= 1
            ("LABEL", "sq_shl1_n"),
            ("LOAD", 0, 0), ("SHL", 0, 6), ("STORE", 0, 0), ("JC", "sqn_c01"), ("MOVI", 2, 0), ("JMP", "sqn_d01"), ("LABEL", "sqn_c01"), ("MOVI", 2, 1), ("LABEL", "sqn_d01"),
            ("LOAD", 0, 2), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "sqn_o1"), ("OR", 0, 6), ("LABEL", "sqn_o1"), ("STORE", 0, 2), ("JC", "sqn_c12"), ("MOVI", 2, 0), ("JMP", "sqn_d12"), ("LABEL", "sqn_c12"), ("MOVI", 2, 1), ("LABEL", "sqn_d12"),
            ("LOAD", 0, 4), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "sqn_o2"), ("OR", 0, 6), ("LABEL", "sqn_o2"), ("STORE", 0, 4), ("JC", "sqn_c23"), ("MOVI", 2, 0), ("JMP", "sqn_d23"), ("LABEL", "sqn_c23"), ("MOVI", 2, 1), ("LABEL", "sqn_d23"),
            ("LOAD", 0, 6), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "sqn_o3"), ("OR", 0, 6), ("LABEL", "sqn_o3"), ("STORE", 0, 6), ("JC", "sqn_c34"), ("MOVI", 2, 0), ("JMP", "sqn_d34"), ("LABEL", "sqn_c34"), ("MOVI", 2, 1), ("LABEL", "sqn_d34"),
            ("LOAD", 0, 8), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "sqn_o4"), ("OR", 0, 6), ("LABEL", "sqn_o4"), ("STORE", 0, 8), ("JC", "sqn_c45"), ("MOVI", 2, 0), ("JMP", "sqn_d45"), ("LABEL", "sqn_c45"), ("MOVI", 2, 1), ("LABEL", "sqn_d45"),
            ("LOAD", 0, 10), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "sqn_o5"), ("OR", 0, 6), ("LABEL", "sqn_o5"), ("STORE", 0, 10), ("JC", "sqn_c56"), ("MOVI", 2, 0), ("JMP", "sqn_d56"), ("LABEL", "sqn_c56"), ("MOVI", 2, 1), ("LABEL", "sqn_d56"),
            ("LOAD", 0, 12), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "sqn_o6"), ("OR", 0, 6), ("LABEL", "sqn_o6"), ("STORE", 0, 12), ("JC", "sqn_c67"), ("MOVI", 2, 0), ("JMP", "sqn_d67"), ("LABEL", "sqn_c67"), ("MOVI", 2, 1), ("LABEL", "sqn_d67"),
            ("LOAD", 0, 14), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "sqn_o7"), ("OR", 0, 6), ("LABEL", "sqn_o7"), ("STORE", 0, 14),
            ("RET",),

            # g(108..114) >>= 1
            ("LABEL", "sq_shr1_g"),
            ("LOAD", 0, 114), ("MOV", 1, 0), ("SHR", 0, 6), ("STORE", 0, 114), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 112), ("MOV", 1, 0), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "sqg_or2"), ("OR", 0, 3), ("LABEL", "sqg_or2"), ("STORE", 0, 112), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 110), ("MOV", 1, 0), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "sqg_or1"), ("OR", 0, 3), ("LABEL", "sqg_or1"), ("STORE", 0, 110), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 108), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "sqg_or0"), ("OR", 0, 3), ("LABEL", "sqg_or0"), ("STORE", 0, 108),
            ("RET",),
        ]

        # append divider as subroutine (prefixed labels, HALT->RET)
        base = Firmware.rows_u128_div_u64()
        label_names = set()
        for row in base:
            if row[0] == "LABEL":
                label_names.add(str(row[1]))
        branch_ops = {"JMP", "JZ", "JNZ", "JN", "JP", "JC", "JNC", "CALL"}
        rows.append(("LABEL", "sq_div_entry"))
        for row in base:
            op = str(row[0])
            if op == "HALT":
                rows.append(("RET",))
                continue
            if op == "LABEL":
                rows.append(("LABEL", f"sqd_{row[1]}"))
                continue
            if op in branch_ops and len(row) > 1 and isinstance(row[1], str) and row[1] in label_names:
                if len(row) == 2:
                    rows.append((op, f"sqd_{row[1]}"))
                else:
                    rows.append((op, f"sqd_{row[1]}", row[2]))
                continue
            rows.append(row)
        return rows

    @staticmethod
    def rows_asin_q46(iterations: int = 32, angle_base: int = 32, sqrt_iterations: int = 16) -> List[AsmRow]:
        # Input:  x_q46 signed at [100..107]
        # Output: asin(x)_q46 signed at [8..15]
        # Requires angle table preloaded by host at [angle_base .. angle_base + iterations*6 - 1].
        rows: List[AsmRow] = [
            ("MOVI", 6, 1),
            ("MOVI", 5, 0),
            ("MOVI", 3, 128),
            ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6),
            ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6),

            # save original x48 for atan y input
            ("LOAD", 0, 100), ("STORE", 0, 224),
            ("LOAD", 0, 102), ("STORE", 0, 226),
            ("LOAD", 0, 104), ("STORE", 0, 228),

            # one_q46 at [120..127]
            ("STORE", 5, 120), ("STORE", 5, 122),
            ("MOVI", 0, 64),
            ("SHL", 0, 6), ("SHL", 0, 6), ("SHL", 0, 6), ("SHL", 0, 6),
            ("SHL", 0, 6), ("SHL", 0, 6), ("SHL", 0, 6), ("SHL", 0, 6),
            ("STORE", 0, 124), ("STORE", 5, 126),

            # xabs = x
            ("LOAD", 0, 100), ("STORE", 0, 132),
            ("LOAD", 0, 102), ("STORE", 0, 134),
            ("LOAD", 0, 104), ("STORE", 0, 136),
            ("LOAD", 0, 106), ("STORE", 0, 138),

            # if negative -> two's complement xabs
            ("LOAD", 0, 106), ("AND", 0, 3), ("CMP", 0, 5), ("JZ", "asin_abs_done"),
            ("MOVI", 4, 0), ("SUB", 4, 6),  # r4=0xFFFF
            ("LOAD", 0, 132), ("XOR", 0, 4), ("STORE", 0, 132),
            ("LOAD", 0, 134), ("XOR", 0, 4), ("STORE", 0, 134),
            ("LOAD", 0, 136), ("XOR", 0, 4), ("STORE", 0, 136),
            ("LOAD", 0, 138), ("XOR", 0, 4), ("STORE", 0, 138),
            ("LOAD", 0, 132), ("ADD", 0, 6), ("STORE", 0, 132),
            ("LOAD", 0, 134), ("ADDC", 0, 5), ("STORE", 0, 134),
            ("LOAD", 0, 136), ("ADDC", 0, 5), ("STORE", 0, 136),
            ("LOAD", 0, 138), ("ADDC", 0, 5), ("STORE", 0, 138),
            ("LABEL", "asin_abs_done"),

            # prepare mul inputs A/B from xabs
            ("LOAD", 0, 132), ("STORE", 0, 0), ("STORE", 0, 8),
            ("LOAD", 0, 134), ("STORE", 0, 2), ("STORE", 0, 10),
            ("LOAD", 0, 136), ("STORE", 0, 4), ("STORE", 0, 12),
            ("LOAD", 0, 138), ("STORE", 0, 6), ("STORE", 0, 14),
            ("CALL", "asin_mul_entry"),

            # (x*x)>>46 in [16..31] via 46 times shr1
            ("MOVI", 0, 46), ("STORE", 0, 130),
            ("LABEL", "asin_shr_loop"),
            ("LOAD", 0, 130), ("CMP", 0, 5), ("JZ", "asin_shr_done"),
            ("CALL", "asin_shr1_128"),
            ("LOAD", 0, 130), ("SUB", 0, 6), ("STORE", 0, 130),
            ("JMP", "asin_shr_loop"),
            ("LABEL", "asin_shr_done"),

            # t = 1 - x2_q46 => [100..107]
            ("LOAD", 0, 120), ("LOAD", 1, 16), ("SUB", 0, 1), ("STORE", 0, 100),
            ("LOAD", 0, 122), ("LOAD", 1, 18), ("SUBB", 0, 1), ("STORE", 0, 102),
            ("LOAD", 0, 124), ("LOAD", 1, 20), ("SUBB", 0, 1), ("STORE", 0, 104),
            ("LOAD", 0, 126), ("LOAD", 1, 22), ("SUBB", 0, 1), ("STORE", 0, 106),

            # root = sqrt(t) -> [8..15]
            ("CALL", "asin_sqrt_entry"),

            # atan input: x=root48, y=orig_x48, z=0
            ("LOAD", 0, 8), ("STORE", 0, 0),
            ("LOAD", 0, 10), ("STORE", 0, 2),
            ("LOAD", 0, 12), ("STORE", 0, 4),
            ("LOAD", 0, 224), ("STORE", 0, 6),
            ("LOAD", 0, 226), ("STORE", 0, 8),
            ("LOAD", 0, 228), ("STORE", 0, 10),
            ("STORE", 5, 12), ("STORE", 5, 14), ("STORE", 5, 16),
            ("CALL", "asin_atan_entry"),

            # out64 sign-extend from z48
            ("LOAD", 0, 12), ("STORE", 0, 8),
            ("LOAD", 0, 14), ("STORE", 0, 10),
            ("LOAD", 0, 16), ("STORE", 0, 12),
            ("LOAD", 0, 16), ("AND", 0, 3), ("CMP", 0, 5), ("JZ", "asin_out_pos"),
            ("MOVI", 0, 0), ("SUB", 0, 6), ("STORE", 0, 14),
            ("JMP", "asin_done"),
            ("LABEL", "asin_out_pos"),
            ("STORE", 5, 14),
            ("LABEL", "asin_done"),
            ("HALT",),

            # shr1 unsigned 128-bit at [16..31]
            ("LABEL", "asin_shr1_128"),
            ("LOAD", 0, 30), ("MOV", 1, 0), ("SHR", 0, 6), ("STORE", 0, 30), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 28), ("MOV", 1, 0), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "as_shr_o6"), ("OR", 0, 3), ("LABEL", "as_shr_o6"), ("STORE", 0, 28), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 26), ("MOV", 1, 0), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "as_shr_o5"), ("OR", 0, 3), ("LABEL", "as_shr_o5"), ("STORE", 0, 26), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 24), ("MOV", 1, 0), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "as_shr_o4"), ("OR", 0, 3), ("LABEL", "as_shr_o4"), ("STORE", 0, 24), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 22), ("MOV", 1, 0), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "as_shr_o3"), ("OR", 0, 3), ("LABEL", "as_shr_o3"), ("STORE", 0, 22), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 20), ("MOV", 1, 0), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "as_shr_o2"), ("OR", 0, 3), ("LABEL", "as_shr_o2"), ("STORE", 0, 20), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 18), ("MOV", 1, 0), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "as_shr_o1"), ("OR", 0, 3), ("LABEL", "as_shr_o1"), ("STORE", 0, 18), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 16), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "as_shr_o0"), ("OR", 0, 3), ("LABEL", "as_shr_o0"), ("STORE", 0, 16),
            ("RET",),
        ]

        Firmware._append_prefixed_subroutine(rows, Firmware.rows_u64_mul_full(), "asin_m_", "asin_mul_entry")
        Firmware._append_prefixed_subroutine(rows, Firmware.rows_sqrt_q46_newton(iterations=sqrt_iterations), "asin_s_", "asin_sqrt_entry")
        Firmware._append_prefixed_subroutine(rows, Firmware.rows_atan_cordic48(iterations=iterations, angle_base=angle_base), "asin_a_", "asin_atan_entry")
        return rows

    @staticmethod
    def rows_acos_q46(iterations: int = 32, angle_base: int = 32, sqrt_iterations: int = 16) -> List[AsmRow]:
        # Input: x_q46 signed at [100..107]
        # Output: acos(x)_q46 signed at [8..15]
        # Requires half_pi_q46 preloaded by host at [240..247]
        rows: List[AsmRow] = [
            ("MOVI", 6, 1),
            ("MOVI", 5, 0),
            ("CALL", "acos_asin_entry"),
            # out = halfpi - asin
            ("LOAD", 0, 240), ("LOAD", 1, 8), ("SUB", 0, 1), ("STORE", 0, 8),
            ("LOAD", 0, 242), ("LOAD", 1, 10), ("SUBB", 0, 1), ("STORE", 0, 10),
            ("LOAD", 0, 244), ("LOAD", 1, 12), ("SUBB", 0, 1), ("STORE", 0, 12),
            ("LOAD", 0, 246), ("LOAD", 1, 14), ("SUBB", 0, 1), ("STORE", 0, 14),
            ("HALT",),
        ]
        Firmware._append_prefixed_subroutine(rows, Firmware.rows_asin_q46(iterations=iterations, angle_base=angle_base, sqrt_iterations=sqrt_iterations), "acos_in_", "acos_asin_entry")
        return rows

    @staticmethod
    def rows_u64_shl1() -> List[AsmRow]:
        # OUT = A << 1, little-endian 16-bit limbs
        # A at [0..7], OUT at [8..15]
        # Uses R6=1, R5=0
        return [
            ("MOVI", 6, 1),
            ("MOVI", 5, 0),

            # word0
            ("LOAD", 0, 0),
            ("SHL", 0, 6),
            ("STORE", 0, 8),
            ("JC", "c01_1"),
            ("MOVI", 2, 0),
            ("JMP", "c01_done"),
            ("LABEL", "c01_1"),
            ("MOVI", 2, 1),
            ("LABEL", "c01_done"),

            # word1
            ("LOAD", 0, 2),
            ("SHL", 0, 6),
            ("CMP", 2, 5),
            ("JZ", "skip_or1"),
            ("OR", 0, 6),
            ("LABEL", "skip_or1"),
            ("STORE", 0, 10),
            ("JC", "c12_1"),
            ("MOVI", 2, 0),
            ("JMP", "c12_done"),
            ("LABEL", "c12_1"),
            ("MOVI", 2, 1),
            ("LABEL", "c12_done"),

            # word2
            ("LOAD", 0, 4),
            ("SHL", 0, 6),
            ("CMP", 2, 5),
            ("JZ", "skip_or2"),
            ("OR", 0, 6),
            ("LABEL", "skip_or2"),
            ("STORE", 0, 12),
            ("JC", "c23_1"),
            ("MOVI", 2, 0),
            ("JMP", "c23_done"),
            ("LABEL", "c23_1"),
            ("MOVI", 2, 1),
            ("LABEL", "c23_done"),

            # word3
            ("LOAD", 0, 6),
            ("SHL", 0, 6),
            ("CMP", 2, 5),
            ("JZ", "skip_or3"),
            ("OR", 0, 6),
            ("LABEL", "skip_or3"),
            ("STORE", 0, 14),
            ("HALT",),
        ]

    @staticmethod
    def rows_u64_cmp() -> List[AsmRow]:
        # unsigned compare
        # A at [0..7], B at [8..15], OUT code at [16..17]:
        # 0 => A<B, 1 => A==B, 2 => A>B
        return [
            ("MOVI", 5, 0),
            ("MOVI", 6, 1),

            # compare w3
            ("LOAD", 0, 6),
            ("LOAD", 1, 14),
            ("MOV", 2, 0),
            ("SUB", 2, 1),
            ("JNC", "less"),
            ("CMP", 2, 5),
            ("JNZ", "greater"),

            # compare w2
            ("LOAD", 0, 4),
            ("LOAD", 1, 12),
            ("MOV", 2, 0),
            ("SUB", 2, 1),
            ("JNC", "less"),
            ("CMP", 2, 5),
            ("JNZ", "greater"),

            # compare w1
            ("LOAD", 0, 2),
            ("LOAD", 1, 10),
            ("MOV", 2, 0),
            ("SUB", 2, 1),
            ("JNC", "less"),
            ("CMP", 2, 5),
            ("JNZ", "greater"),

            # compare w0
            ("LOAD", 0, 0),
            ("LOAD", 1, 8),
            ("MOV", 2, 0),
            ("SUB", 2, 1),
            ("JNC", "less"),
            ("CMP", 2, 5),
            ("JNZ", "greater"),

            ("STORE", 6, 16),    # equal => 1
            ("HALT",),

            ("LABEL", "less"),
            ("STORE", 5, 16),    # less => 0
            ("HALT",),

            ("LABEL", "greater"),
            ("MOVI", 0, 2),
            ("STORE", 0, 16),    # greater => 2
            ("HALT",),
        ]

    @staticmethod
    def rows_u64_clz() -> List[AsmRow]:
        # count leading zeros for 64-bit unsigned
        # A at [0..7], OUT count at [8..9] in [0,64]
        return [
            ("MOVI", 0, 0),       # count
            ("MOVI", 5, 0),       # zero
            ("MOVI", 6, 1),       # one
            ("MOVI", 7, 16),      # sixteen
            ("MOVI", 4, 128),     # build mask 0x8000 in R4
            ("SHL", 4, 6), ("SHL", 4, 6), ("SHL", 4, 6), ("SHL", 4, 6),
            ("SHL", 4, 6), ("SHL", 4, 6), ("SHL", 4, 6), ("SHL", 4, 6),

            ("LOAD", 1, 6),
            ("CMP", 1, 5),
            ("JNZ", "scan_bits"),
            ("ADD", 0, 7),

            ("LOAD", 1, 4),
            ("CMP", 1, 5),
            ("JNZ", "scan_bits"),
            ("ADD", 0, 7),

            ("LOAD", 1, 2),
            ("CMP", 1, 5),
            ("JNZ", "scan_bits"),
            ("ADD", 0, 7),

            ("LOAD", 1, 0),
            ("CMP", 1, 5),
            ("JNZ", "scan_bits"),
            ("ADD", 0, 7),        # now 64
            ("STORE", 0, 8),
            ("HALT",),

            ("LABEL", "scan_bits"),
            ("MOV", 2, 1),        # temp = word
            ("AND", 2, 4),        # temp &= 0x8000
            ("CMP", 2, 5),
            ("JNZ", "done"),
            ("SHL", 1, 6),        # word <<= 1
            ("ADD", 0, 6),        # count++
            ("JMP", "scan_bits"),

            ("LABEL", "done"),
            ("STORE", 0, 8),
            ("HALT",),
        ]

    @staticmethod
    def rows_soft64_add_sameexp_pos() -> List[AsmRow]:
        # Improved prototype: positive normal add with arbitrary exponent gap.
        # Still not full IEEE (no negative/subnormal/NaN/Inf handling in firmware).
        #
        # Layout:
        #   A bits [0..7], B bits [8..15], OUT bits [16..23]
        #   MA [24..31], MB [32..39], SUM [40..47], exp_big [60..61]
        return [
            ("MOVI", 6, 1),   # one
            ("MOVI", 5, 0),   # zero
            ("MOVI", 4, 15),  # mask low4

            # r3 = 0x8000 mask
            ("MOVI", 3, 128),
            ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6),
            ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6),

            # ea -> r0
            ("LOAD", 0, 6),
            ("SHR", 0, 6), ("SHR", 0, 6), ("SHR", 0, 6), ("SHR", 0, 6),
            # eb -> r1
            ("LOAD", 1, 14),
            ("SHR", 1, 6), ("SHR", 1, 6), ("SHR", 1, 6), ("SHR", 1, 6),

            # build MA from A
            ("LOAD", 2, 0), ("STORE", 2, 24),
            ("LOAD", 2, 2), ("STORE", 2, 26),
            ("LOAD", 2, 4), ("STORE", 2, 28),
            ("LOAD", 2, 6), ("AND", 2, 4),
            ("MOVI", 1, 16), ("OR", 2, 1), ("STORE", 2, 30),

            # build MB from B
            ("LOAD", 2, 8), ("STORE", 2, 32),
            ("LOAD", 2, 10), ("STORE", 2, 34),
            ("LOAD", 2, 12), ("STORE", 2, 36),
            ("LOAD", 2, 14), ("AND", 2, 4),
            ("MOVI", 1, 16), ("OR", 2, 1), ("STORE", 2, 38),

            # reload exponents (r0=ea, r1=eb) after mantissa build clobbers
            ("LOAD", 0, 6),
            ("SHR", 0, 6), ("SHR", 0, 6), ("SHR", 0, 6), ("SHR", 0, 6),
            ("LOAD", 1, 14),
            ("SHR", 1, 6), ("SHR", 1, 6), ("SHR", 1, 6), ("SHR", 1, 6),

            # if ea < eb => swap (exp + mantissas)
            ("CMP", 0, 1),
            ("JN", "add_swap"),
            ("JMP", "add_after_swap"),

            ("LABEL", "add_swap"),
            ("MOV", 2, 0), ("MOV", 0, 1), ("MOV", 1, 2),
            ("LOAD", 2, 24), ("LOAD", 6, 32), ("STORE", 6, 24), ("STORE", 2, 32),
            ("LOAD", 2, 26), ("LOAD", 6, 34), ("STORE", 6, 26), ("STORE", 2, 34),
            ("LOAD", 2, 28), ("LOAD", 6, 36), ("STORE", 6, 28), ("STORE", 2, 36),
            ("LOAD", 2, 30), ("LOAD", 6, 38), ("STORE", 6, 30), ("STORE", 2, 38),
            ("MOVI", 6, 1),

            ("LABEL", "add_after_swap"),
            ("STORE", 0, 60),             # exp_big
            ("MOV", 2, 0), ("SUB", 2, 1), # diff = ea-eb
            ("MOVI", 1, 63), ("CMP", 2, 1), ("JN", "add_diff_ok"), ("MOV", 2, 1), ("LABEL", "add_diff_ok"),
            ("STORE", 2, 62),
            ("STORE", 5, 64),             # guard bit from alignment
            ("STORE", 5, 66),             # sticky bits from alignment

            # align MB by diff
            ("LABEL", "add_align_loop"),
            ("LOAD", 2, 62),
            ("CMP", 2, 5),
            ("JZ", "add_align_done"),
            ("LOAD", 0, 32), ("AND", 0, 6),               # dropped bit this shift
            ("LOAD", 1, 62), ("CMP", 1, 6), ("JZ", "add_cap_guard"),
            ("LOAD", 1, 66), ("OR", 1, 0), ("STORE", 1, 66), ("JMP", "add_cap_done"),
            ("LABEL", "add_cap_guard"),
            ("STORE", 0, 64),
            ("LABEL", "add_cap_done"),
            ("CALL", "add_shr1_mb"),
            ("LOAD", 2, 62), ("SUB", 2, 6), ("STORE", 2, 62),
            ("JMP", "add_align_loop"),
            ("LABEL", "add_align_done"),

            # SUM = MA + MB
            ("LOAD", 0, 24), ("LOAD", 1, 32), ("ADD", 0, 1),  ("STORE", 0, 40),
            ("LOAD", 0, 26), ("LOAD", 1, 34), ("ADDC", 0, 1), ("STORE", 0, 42),
            ("LOAD", 0, 28), ("LOAD", 1, 36), ("ADDC", 0, 1), ("STORE", 0, 44),
            ("LOAD", 0, 30), ("LOAD", 1, 38), ("ADDC", 0, 1), ("STORE", 0, 46),

            # overflow bit53? (word3 bit5 set)
            ("MOVI", 2, 32),
            ("LOAD", 1, 46), ("AND", 1, 2),
            ("CMP", 1, 5),
            ("JZ", "add_no_ovf"),
            ("LOAD", 0, 40), ("AND", 0, 6), ("STORE", 0, 68),  # bit dropped by sum>>1 becomes new guard
            ("CALL", "add_shr1_sum"),
            ("LOAD", 0, 66), ("LOAD", 1, 64), ("OR", 0, 1), ("STORE", 0, 66),  # sticky |= old guard
            ("LOAD", 0, 68), ("STORE", 0, 64),                                     # guard = dropped sum bit
            ("LOAD", 0, 60), ("ADD", 0, 6), ("STORE", 0, 60),  # exp_big++

            ("LABEL", "add_no_ovf"),
            # round-to-nearest-even using guard/sticky and current mantissa LSB.
            ("LOAD", 0, 64), ("CMP", 0, 5), ("JZ", "add_after_round"),
            ("LOAD", 0, 66), ("CMP", 0, 5), ("JNZ", "add_round_inc"),
            ("LOAD", 0, 40), ("AND", 0, 6), ("CMP", 0, 5), ("JZ", "add_after_round"),
            ("LABEL", "add_round_inc"),
            ("LOAD", 0, 40), ("ADD", 0, 6), ("STORE", 0, 40),
            ("LOAD", 0, 42), ("ADDC", 0, 5), ("STORE", 0, 42),
            ("LOAD", 0, 44), ("ADDC", 0, 5), ("STORE", 0, 44),
            ("LOAD", 0, 46), ("ADDC", 0, 5), ("STORE", 0, 46),
            # if rounding caused bit53 overflow, renormalize and bump exponent.
            ("MOVI", 2, 32),
            ("LOAD", 1, 46), ("AND", 1, 2), ("CMP", 1, 5), ("JZ", "add_after_round"),
            ("CALL", "add_shr1_sum"),
            ("LOAD", 0, 60), ("ADD", 0, 6), ("STORE", 0, 60),
            ("LABEL", "add_after_round"),
            # exp overflow to 2047 => +Inf
            ("LOAD", 0, 60),
            ("MOVI", 2, 127), ("SHL", 2, 6), ("SHL", 2, 6), ("SHL", 2, 6), ("SHL", 2, 6), ("OR", 2, 4),
            ("CMP", 0, 2), ("JNZ", "add_compose"),
            ("STORE", 5, 16), ("STORE", 5, 18), ("STORE", 5, 20),
            ("SHL", 2, 6), ("SHL", 2, 6), ("SHL", 2, 6), ("SHL", 2, 6),
            ("STORE", 2, 22),
            ("HALT",),

            ("LABEL", "add_compose"),
            # compose output positive sign
            ("LOAD", 0, 60),
            ("SHL", 0, 6), ("SHL", 0, 6), ("SHL", 0, 6), ("SHL", 0, 6),  # exp<<4
            ("LOAD", 1, 46), ("AND", 1, 4), ("OR", 0, 1),
            ("LOAD", 1, 40), ("STORE", 1, 16),
            ("LOAD", 1, 42), ("STORE", 1, 18),
            ("LOAD", 1, 44), ("STORE", 1, 20),
            ("STORE", 0, 22),
            ("HALT",),

            # subroutine: MB >>= 1  (addr 32..38)
            ("LABEL", "add_shr1_mb"),
            ("LOAD", 0, 38), ("MOV", 1, 0), ("SHR", 0, 6), ("STORE", 0, 38), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 36), ("MOV", 1, 0), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "mb_skip_or2"), ("OR", 0, 3), ("LABEL", "mb_skip_or2"), ("STORE", 0, 36), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 34), ("MOV", 1, 0), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "mb_skip_or1"), ("OR", 0, 3), ("LABEL", "mb_skip_or1"), ("STORE", 0, 34), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 32), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "mb_skip_or0"), ("OR", 0, 3), ("LABEL", "mb_skip_or0"), ("STORE", 0, 32),
            ("RET",),

            # subroutine: SUM >>= 1 (addr 40..46)
            ("LABEL", "add_shr1_sum"),
            ("LOAD", 0, 46), ("MOV", 1, 0), ("SHR", 0, 6), ("STORE", 0, 46), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 44), ("MOV", 1, 0), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "sm_skip_or2"), ("OR", 0, 3), ("LABEL", "sm_skip_or2"), ("STORE", 0, 44), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 42), ("MOV", 1, 0), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "sm_skip_or1"), ("OR", 0, 3), ("LABEL", "sm_skip_or1"), ("STORE", 0, 42), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 40), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "sm_skip_or0"), ("OR", 0, 3), ("LABEL", "sm_skip_or0"), ("STORE", 0, 40),
            ("RET",),
        ]

    @staticmethod
    def rows_soft64_sub_align_pos() -> List[AsmRow]:
        # Magnitude subtraction core for opposite-sign add:
        # OUT = |A| - |B| with alignment and normalization.
        # Input A/B are positive normal bit patterns in [0..7]/[8..15].
        # Output magnitude in [16..23].
        return [
            ("MOVI", 6, 1),   # one
            ("MOVI", 5, 0),   # zero
            ("MOVI", 4, 15),  # mask low4
            ("MOVI", 3, 128), # build 0x8000 in r3
            ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6),
            ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6),

            # ea/eb
            ("LOAD", 0, 6),
            ("SHR", 0, 6), ("SHR", 0, 6), ("SHR", 0, 6), ("SHR", 0, 6),
            ("LOAD", 1, 14),
            ("SHR", 1, 6), ("SHR", 1, 6), ("SHR", 1, 6), ("SHR", 1, 6),

            # MA
            ("LOAD", 2, 0), ("STORE", 2, 24),
            ("LOAD", 2, 2), ("STORE", 2, 26),
            ("LOAD", 2, 4), ("STORE", 2, 28),
            ("LOAD", 2, 6), ("AND", 2, 4), ("MOVI", 1, 16), ("OR", 2, 1), ("STORE", 2, 30),
            # MB
            ("LOAD", 2, 8), ("STORE", 2, 32),
            ("LOAD", 2, 10), ("STORE", 2, 34),
            ("LOAD", 2, 12), ("STORE", 2, 36),
            ("LOAD", 2, 14), ("AND", 2, 4), ("MOVI", 1, 16), ("OR", 2, 1), ("STORE", 2, 38),

            # reload ea/eb
            ("LOAD", 0, 6),
            ("SHR", 0, 6), ("SHR", 0, 6), ("SHR", 0, 6), ("SHR", 0, 6),
            ("LOAD", 1, 14),
            ("SHR", 1, 6), ("SHR", 1, 6), ("SHR", 1, 6), ("SHR", 1, 6),

            # ensure ea>=eb by swap if needed
            ("CMP", 0, 1),
            ("JN", "sub_swap_exp"),
            ("JMP", "sub_after_exp"),
            ("LABEL", "sub_swap_exp"),
            ("MOV", 2, 0), ("MOV", 0, 1), ("MOV", 1, 2),
            ("LOAD", 2, 24), ("LOAD", 6, 32), ("STORE", 6, 24), ("STORE", 2, 32),
            ("LOAD", 2, 26), ("LOAD", 6, 34), ("STORE", 6, 26), ("STORE", 2, 34),
            ("LOAD", 2, 28), ("LOAD", 6, 36), ("STORE", 6, 28), ("STORE", 2, 36),
            ("LOAD", 2, 30), ("LOAD", 6, 38), ("STORE", 6, 30), ("STORE", 2, 38),
            ("MOVI", 6, 1),
            ("LABEL", "sub_after_exp"),

            ("STORE", 0, 60),             # exp_big
            ("STORE", 0, 68),             # exp_orig
            ("MOV", 2, 0), ("SUB", 2, 1), # diff
            ("MOVI", 1, 63), ("CMP", 2, 1), ("JN", "sub_diff_ok"), ("MOV", 2, 1), ("LABEL", "sub_diff_ok"),
            ("STORE", 2, 62),
            ("STORE", 5, 64),             # guard bit from alignment
            ("STORE", 5, 66),             # sticky bits from alignment
            ("STORE", 5, 70),             # pre-round applied flag

            # align MB by diff
            ("LABEL", "sub_align_loop"),
            ("LOAD", 2, 62), ("CMP", 2, 5), ("JZ", "sub_align_done"),
            ("LOAD", 0, 32), ("AND", 0, 6),               # dropped bit this shift
            ("LOAD", 1, 62), ("CMP", 1, 6), ("JZ", "sub_cap_guard"),
            ("LOAD", 1, 66), ("OR", 1, 0), ("STORE", 1, 66), ("JMP", "sub_cap_done"),
            ("LABEL", "sub_cap_guard"),
            ("STORE", 0, 64),
            ("LABEL", "sub_cap_done"),
            ("CALL", "sub_shr1_mb"),
            ("LOAD", 2, 62), ("SUB", 2, 6), ("STORE", 2, 62),
            ("JMP", "sub_align_loop"),
            ("LABEL", "sub_align_done"),

            # if MA < MB then swap magnitudes
            ("CALL", "sub_cmp_ma_mb"),
            ("LOAD", 2, 58), ("CMP", 2, 5), ("JZ", "sub_no_swap_mag"),
            ("LOAD", 2, 24), ("LOAD", 1, 32), ("STORE", 1, 24), ("STORE", 2, 32),
            ("LOAD", 2, 26), ("LOAD", 1, 34), ("STORE", 1, 26), ("STORE", 2, 34),
            ("LOAD", 2, 28), ("LOAD", 1, 36), ("STORE", 1, 28), ("STORE", 2, 36),
            ("LOAD", 2, 30), ("LOAD", 1, 38), ("STORE", 1, 30), ("STORE", 2, 38),
            ("LABEL", "sub_no_swap_mag"),

            # SUM = MA - MB
            ("LOAD", 0, 24), ("LOAD", 1, 32), ("SUB", 0, 1),  ("STORE", 0, 40),
            ("LOAD", 0, 26), ("LOAD", 1, 34), ("SUBB", 0, 1), ("STORE", 0, 42),
            ("LOAD", 0, 28), ("LOAD", 1, 36), ("SUBB", 0, 1), ("STORE", 0, 44),
            ("LOAD", 0, 30), ("LOAD", 1, 38), ("SUBB", 0, 1), ("STORE", 0, 46),

            # if SUM==0 => output zero
            ("LOAD", 0, 40), ("LOAD", 1, 42), ("OR", 0, 1), ("LOAD", 1, 44), ("OR", 0, 1), ("LOAD", 1, 46), ("OR", 0, 1),
            ("CMP", 0, 5), ("JNZ", "sub_norm_start"),
            ("STORE", 5, 16), ("STORE", 5, 18), ("STORE", 5, 20), ("STORE", 5, 22),
            ("HALT",),

            # RN-even correction at exp_big scale before normalization:
            # alignment truncation makes SUM too large; subtract 1 ulp when needed.
            ("LABEL", "sub_norm_start"),
            ("LOAD", 0, 64), ("CMP", 0, 5), ("JZ", "sub_norm_loop"),
            ("LOAD", 0, 66), ("CMP", 0, 5), ("JNZ", "sub_round_dec"),
            ("LOAD", 0, 40), ("AND", 0, 6), ("CMP", 0, 5), ("JZ", "sub_norm_loop"),
            ("LABEL", "sub_round_dec"),
            ("STORE", 6, 70),             # pre-round applied
            ("LOAD", 0, 40), ("SUB", 0, 6), ("STORE", 0, 40),
            ("LOAD", 0, 42), ("SUBB", 0, 5), ("STORE", 0, 42),
            ("LOAD", 0, 44), ("SUBB", 0, 5), ("STORE", 0, 44),
            ("LOAD", 0, 46), ("SUBB", 0, 5), ("STORE", 0, 46),
            ("LOAD", 0, 40), ("LOAD", 1, 42), ("OR", 0, 1), ("LOAD", 1, 44), ("OR", 0, 1), ("LOAD", 1, 46), ("OR", 0, 1),
            ("CMP", 0, 5), ("JNZ", "sub_norm_loop"),
            ("STORE", 5, 16), ("STORE", 5, 18), ("STORE", 5, 20), ("STORE", 5, 22),
            ("HALT",),

            # normalize left until bit52 present in word3 bit4
            ("LABEL", "sub_norm_loop"),
            ("MOVI", 2, 16),          # bit52 mask in top word (R2 may be clobbered by subroutine)
            ("LOAD", 1, 46), ("AND", 1, 2), ("CMP", 1, 5), ("JNZ", "sub_norm_done"),
            ("LOAD", 0, 60), ("CMP", 0, 6), ("JZ", "sub_norm_done"),
            ("CALL", "sub_shl1_sum"),
            ("LOAD", 0, 60), ("SUB", 0, 6), ("STORE", 0, 60),
            ("JMP", "sub_norm_loop"),
            ("LABEL", "sub_norm_done"),

            # Tie corner for one-bit cancellation shift (k==1) that pre-round at exp_big misses.
            ("LOAD", 0, 70), ("CMP", 0, 5), ("JNZ", "sub_post_round_done"),
            ("LOAD", 0, 64), ("CMP", 0, 6), ("JNZ", "sub_post_round_done"),  # guard==1
            ("LOAD", 0, 66), ("CMP", 0, 5), ("JNZ", "sub_post_round_done"),  # sticky==0
            ("LOAD", 0, 68), ("LOAD", 1, 60), ("SUB", 0, 1), ("CMP", 0, 6), ("JNZ", "sub_post_round_done"),  # k==1
            ("LOAD", 0, 40), ("SUB", 0, 6), ("STORE", 0, 40),
            ("LOAD", 0, 42), ("SUBB", 0, 5), ("STORE", 0, 42),
            ("LOAD", 0, 44), ("SUBB", 0, 5), ("STORE", 0, 44),
            ("LOAD", 0, 46), ("SUBB", 0, 5), ("STORE", 0, 46),
            ("LABEL", "sub_post_renorm"),
            ("MOVI", 2, 16),
            ("LOAD", 1, 46), ("AND", 1, 2), ("CMP", 1, 5), ("JNZ", "sub_post_round_done"),
            ("LOAD", 0, 60), ("CMP", 0, 6), ("JZ", "sub_post_round_done"),
            ("CALL", "sub_shl1_sum"),
            ("LOAD", 0, 60), ("SUB", 0, 6), ("STORE", 0, 60),
            ("JMP", "sub_post_renorm"),
            ("LABEL", "sub_post_round_done"),

            # compose output
            ("LOAD", 0, 60),
            ("SHL", 0, 6), ("SHL", 0, 6), ("SHL", 0, 6), ("SHL", 0, 6),
            ("LOAD", 1, 46), ("AND", 1, 4), ("OR", 0, 1),
            ("LOAD", 1, 40), ("STORE", 1, 16),
            ("LOAD", 1, 42), ("STORE", 1, 18),
            ("LOAD", 1, 44), ("STORE", 1, 20),
            ("STORE", 0, 22),
            ("HALT",),

            # MB >>= 1 (32..38)
            ("LABEL", "sub_shr1_mb"),
            ("LOAD", 0, 38), ("MOV", 1, 0), ("SHR", 0, 6), ("STORE", 0, 38), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 36), ("MOV", 1, 0), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "smb_or2"), ("OR", 0, 3), ("LABEL", "smb_or2"), ("STORE", 0, 36), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 34), ("MOV", 1, 0), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "smb_or1"), ("OR", 0, 3), ("LABEL", "smb_or1"), ("STORE", 0, 34), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 32), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "smb_or0"), ("OR", 0, 3), ("LABEL", "smb_or0"), ("STORE", 0, 32),
            ("RET",),

            # SUM <<= 1 (40..46)
            ("LABEL", "sub_shl1_sum"),
            ("LOAD", 0, 40), ("SHL", 0, 6), ("STORE", 0, 40), ("JC", "ss_c01_1"), ("MOVI", 2, 0), ("JMP", "ss_c01_d"), ("LABEL", "ss_c01_1"), ("MOVI", 2, 1), ("LABEL", "ss_c01_d"),
            ("LOAD", 0, 42), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "ss_or1"), ("OR", 0, 6), ("LABEL", "ss_or1"), ("STORE", 0, 42), ("JC", "ss_c12_1"), ("MOVI", 2, 0), ("JMP", "ss_c12_d"), ("LABEL", "ss_c12_1"), ("MOVI", 2, 1), ("LABEL", "ss_c12_d"),
            ("LOAD", 0, 44), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "ss_or2"), ("OR", 0, 6), ("LABEL", "ss_or2"), ("STORE", 0, 44), ("JC", "ss_c23_1"), ("MOVI", 2, 0), ("JMP", "ss_c23_d"), ("LABEL", "ss_c23_1"), ("MOVI", 2, 1), ("LABEL", "ss_c23_d"),
            ("LOAD", 0, 46), ("SHL", 0, 6), ("CMP", 2, 5), ("JZ", "ss_or3"), ("OR", 0, 6), ("LABEL", "ss_or3"), ("STORE", 0, 46),
            ("RET",),

            # compare MA vs MB unsigned (sets RAM58: 0 if MA>=MB, 1 if MA<MB)
            ("LABEL", "sub_cmp_ma_mb"),
            ("LOAD", 0, 30), ("LOAD", 1, 38), ("MOV", 2, 0), ("SUB", 2, 1), ("JNC", "sc_less"), ("CMP", 2, 5), ("JNZ", "sc_ge"),
            ("LOAD", 0, 28), ("LOAD", 1, 36), ("MOV", 2, 0), ("SUB", 2, 1), ("JNC", "sc_less"), ("CMP", 2, 5), ("JNZ", "sc_ge"),
            ("LOAD", 0, 26), ("LOAD", 1, 34), ("MOV", 2, 0), ("SUB", 2, 1), ("JNC", "sc_less"), ("CMP", 2, 5), ("JNZ", "sc_ge"),
            ("LOAD", 0, 24), ("LOAD", 1, 32), ("MOV", 2, 0), ("SUB", 2, 1), ("JNC", "sc_less"),
            ("LABEL", "sc_ge"),
            ("STORE", 5, 58), ("RET",),
            ("LABEL", "sc_less"),
            ("STORE", 6, 58), ("RET",),
        ]

    @staticmethod
    def rows_fp64_add_norm() -> List[AsmRow]:
        # Limited fp64 add in pure VM firmware:
        # - handles finite normal inputs (all sign combinations)
        # - no NaN/Inf/subnormal handling
        # Layout:
        #   A bits [0..7], B bits [8..15], OUT [16..23]
        #   sa [90], sb [92], cmp_code [94] (0<,1==,2>)
        rows: List[AsmRow] = [
            ("MOVI", 5, 0),   # zero
            ("MOVI", 6, 1),   # one
            # r3=0x8000 (sign bit)
            ("MOVI", 3, 128),
            ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6),
            ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6),
            # r4=0x7FFF
            ("MOVI", 4, 0), ("SUB", 4, 6), ("XOR", 4, 3),

            # sa
            ("LOAD", 0, 6), ("AND", 0, 3), ("CMP", 0, 5), ("JZ", "fadd_sa0"),
            ("STORE", 6, 90), ("JMP", "fadd_sa_done"),
            ("LABEL", "fadd_sa0"), ("STORE", 5, 90),
            ("LABEL", "fadd_sa_done"),
            # sb
            ("LOAD", 0, 14), ("AND", 0, 3), ("CMP", 0, 5), ("JZ", "fadd_sb0"),
            ("STORE", 6, 92), ("JMP", "fadd_sb_done"),
            ("LABEL", "fadd_sb0"), ("STORE", 5, 92),
            ("LABEL", "fadd_sb_done"),

            # clear sign bits in A/B (magnitude only)
            ("LOAD", 0, 6), ("AND", 0, 4), ("STORE", 0, 6),
            ("LOAD", 0, 14), ("AND", 0, 4), ("STORE", 0, 14),

            # if sa == sb -> magnitude add
            ("LOAD", 0, 90), ("LOAD", 1, 92), ("CMP", 0, 1), ("JNZ", "fadd_diffsign"),
            ("CALL", "fadd_add_entry"),
            ("LOAD", 0, 90), ("CMP", 0, 5), ("JZ", "fadd_done"),
            ("LOAD", 0, 22), ("OR", 0, 3), ("STORE", 0, 22),
            ("JMP", "fadd_done"),

            ("LABEL", "fadd_diffsign"),
            # compare |A| vs |B| unsigned into [94]
            ("CALL", "fadd_cmp_entry"),
            ("LOAD", 0, 16), ("STORE", 0, 94),
            # OUT = |A|-|B| (abs, subroutine handles ordering)
            ("CALL", "fadd_sub_entry"),

            # if equal, force +0
            ("LOAD", 0, 94), ("MOVI", 1, 2), ("CMP", 0, 1), ("JZ", "fadd_gt"),
            ("CMP", 0, 5), ("JZ", "fadd_lt"),
            ("LOAD", 0, 16), ("LOAD", 1, 18), ("OR", 0, 1), ("LOAD", 1, 20), ("OR", 0, 1), ("LOAD", 1, 22), ("OR", 0, 1),
            ("CMP", 0, 5), ("JNZ", "fadd_done"),
            ("STORE", 5, 22), ("JMP", "fadd_done"),

            ("LABEL", "fadd_gt"),
            ("LOAD", 0, 90), ("CMP", 0, 5), ("JZ", "fadd_done"),
            ("LOAD", 0, 22), ("OR", 0, 3), ("STORE", 0, 22),
            ("JMP", "fadd_done"),

            ("LABEL", "fadd_lt"),
            ("LOAD", 0, 92), ("CMP", 0, 5), ("JZ", "fadd_done"),
            ("LOAD", 0, 22), ("OR", 0, 3), ("STORE", 0, 22),

            ("LABEL", "fadd_done"),
            ("HALT",),
        ]

        Firmware._append_prefixed_subroutine(rows, Firmware.rows_u64_cmp(), "fadd_c_", "fadd_cmp_entry")
        Firmware._append_prefixed_subroutine(rows, Firmware.rows_soft64_add_sameexp_pos(), "fadd_a_", "fadd_add_entry")
        Firmware._append_prefixed_subroutine(rows, Firmware.rows_soft64_sub_align_pos(), "fadd_s_", "fadd_sub_entry")
        return rows

    @staticmethod
    def rows_sin_cordic48(iterations: int = 32, angle_base: int = 64) -> List[AsmRow]:
        # 48-bit signed fixed-point CORDIC rotation mode for sin(x).
        # Memory layout (little-endian 16-bit words):
        # x: [0,2,4], y: [6,8,10], z: [12,14,16], temp: [18,20,22], x_old: [24,26,28]
        # angles table starts at angle_base, each angle is 3 words (48-bit signed).
        if iterations < 1:
            raise ValueError("iterations must be >= 1")
        if angle_base < 0 or angle_base + iterations * 6 > VM16.RAM_SIZE:
            raise ValueError(f"angles table does not fit in {VM16.RAM_SIZE}B RAM")

        rows: List[AsmRow] = [
            ("MOVI", 6, 1),   # one
            ("MOVI", 5, 0),   # zero
            # r3 = 0x8000
            ("MOVI", 3, 128),
            ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6),
            ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6),
        ]

        for i in range(iterations):
            a0 = angle_base + i * 6
            a1 = a0 + 2
            a2 = a0 + 4
            pos_label = f"sin_pos_{i}"
            end_label = f"sin_end_{i}"

            rows.extend([
                # if z >= 0 => positive branch
                ("LOAD", 0, 16),
                ("AND", 0, 3),
                ("CMP", 0, 5),
                ("JZ", pos_label),

                # negative branch: x += (y>>i); y -= (x_old>>i); z += angle[i]
                ("LOAD", 0, 0), ("STORE", 0, 24),
                ("LOAD", 0, 2), ("STORE", 0, 26),
                ("LOAD", 0, 4), ("STORE", 0, 28),

                ("LOAD", 0, 6), ("STORE", 0, 18),
                ("LOAD", 0, 8), ("STORE", 0, 20),
                ("LOAD", 0, 10), ("STORE", 0, 22),
            ])
            for _ in range(i):
                rows.append(("CALL", "sin_sar1_temp"))
            rows.extend([
                ("LOAD", 0, 0), ("LOAD", 1, 18), ("ADD", 0, 1), ("STORE", 0, 0),
                ("LOAD", 0, 2), ("LOAD", 1, 20), ("ADDC", 0, 1), ("STORE", 0, 2),
                ("LOAD", 0, 4), ("LOAD", 1, 22), ("ADDC", 0, 1), ("STORE", 0, 4),

                ("LOAD", 0, 24), ("STORE", 0, 18),
                ("LOAD", 0, 26), ("STORE", 0, 20),
                ("LOAD", 0, 28), ("STORE", 0, 22),
            ])
            for _ in range(i):
                rows.append(("CALL", "sin_sar1_temp"))
            rows.extend([
                ("LOAD", 0, 6), ("LOAD", 1, 18), ("SUB", 0, 1), ("STORE", 0, 6),
                ("LOAD", 0, 8), ("LOAD", 1, 20), ("SUBB", 0, 1), ("STORE", 0, 8),
                ("LOAD", 0, 10), ("LOAD", 1, 22), ("SUBB", 0, 1), ("STORE", 0, 10),

                ("LOAD", 0, 12), ("LOAD", 1, a0), ("ADD", 0, 1), ("STORE", 0, 12),
                ("LOAD", 0, 14), ("LOAD", 1, a1), ("ADDC", 0, 1), ("STORE", 0, 14),
                ("LOAD", 0, 16), ("LOAD", 1, a2), ("ADDC", 0, 1), ("STORE", 0, 16),
                ("JMP", end_label),

                ("LABEL", pos_label),
                # positive branch: x -= (y>>i); y += (x_old>>i); z -= angle[i]
                ("LOAD", 0, 0), ("STORE", 0, 24),
                ("LOAD", 0, 2), ("STORE", 0, 26),
                ("LOAD", 0, 4), ("STORE", 0, 28),

                ("LOAD", 0, 6), ("STORE", 0, 18),
                ("LOAD", 0, 8), ("STORE", 0, 20),
                ("LOAD", 0, 10), ("STORE", 0, 22),
            ])
            for _ in range(i):
                rows.append(("CALL", "sin_sar1_temp"))
            rows.extend([
                ("LOAD", 0, 0), ("LOAD", 1, 18), ("SUB", 0, 1), ("STORE", 0, 0),
                ("LOAD", 0, 2), ("LOAD", 1, 20), ("SUBB", 0, 1), ("STORE", 0, 2),
                ("LOAD", 0, 4), ("LOAD", 1, 22), ("SUBB", 0, 1), ("STORE", 0, 4),

                ("LOAD", 0, 24), ("STORE", 0, 18),
                ("LOAD", 0, 26), ("STORE", 0, 20),
                ("LOAD", 0, 28), ("STORE", 0, 22),
            ])
            for _ in range(i):
                rows.append(("CALL", "sin_sar1_temp"))
            rows.extend([
                ("LOAD", 0, 6), ("LOAD", 1, 18), ("ADD", 0, 1), ("STORE", 0, 6),
                ("LOAD", 0, 8), ("LOAD", 1, 20), ("ADDC", 0, 1), ("STORE", 0, 8),
                ("LOAD", 0, 10), ("LOAD", 1, 22), ("ADDC", 0, 1), ("STORE", 0, 10),

                ("LOAD", 0, 12), ("LOAD", 1, a0), ("SUB", 0, 1), ("STORE", 0, 12),
                ("LOAD", 0, 14), ("LOAD", 1, a1), ("SUBB", 0, 1), ("STORE", 0, 14),
                ("LOAD", 0, 16), ("LOAD", 1, a2), ("SUBB", 0, 1), ("STORE", 0, 16),

                ("LABEL", end_label),
            ])

        rows.extend([
            ("HALT",),

            # temp (18..22) arithmetic right shift by 1 (signed 48-bit)
            ("LABEL", "sin_sar1_temp"),
            ("LOAD", 0, 22), ("MOV", 1, 0), ("SAR", 0, 6), ("STORE", 0, 22), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 20), ("MOV", 1, 0), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "sin_sar1_no_or1"), ("OR", 0, 3), ("LABEL", "sin_sar1_no_or1"), ("STORE", 0, 20), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 18), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "sin_sar1_no_or0"), ("OR", 0, 3), ("LABEL", "sin_sar1_no_or0"), ("STORE", 0, 18),
            ("RET",),
        ])
        return rows

    @staticmethod
    def rows_atan_cordic48(iterations: int = 32, angle_base: int = 64) -> List[AsmRow]:
        # 48-bit signed fixed-point CORDIC vectoring mode for atan(y/x).
        # Input layout:
        #   x: [0,2,4], y: [6,8,10], z: [12,14,16]
        #   z should be initialized to 0.
        # Temp:
        #   temp [18,20,22], x_old [24,26,28]
        if iterations < 1:
            raise ValueError("iterations must be >= 1")
        if angle_base < 0 or angle_base + iterations * 6 > VM16.RAM_SIZE:
            raise ValueError(f"angles table does not fit in {VM16.RAM_SIZE}B RAM")

        rows: List[AsmRow] = [
            ("MOVI", 6, 1),   # one
            ("MOVI", 5, 0),   # zero
            # r3 = 0x8000
            ("MOVI", 3, 128),
            ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6),
            ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6), ("SHL", 3, 6),
        ]

        for i in range(iterations):
            a0 = angle_base + i * 6
            a1 = a0 + 2
            a2 = a0 + 4
            pos_label = f"atan_pos_{i}"
            end_label = f"atan_end_{i}"

            rows.extend([
                # if y == 0, done early to avoid quantized tail-angle bias
                ("LOAD", 0, 6), ("LOAD", 1, 8), ("OR", 0, 1), ("LOAD", 1, 10), ("OR", 0, 1), ("CMP", 0, 5), ("JZ", "atan_done"),

                # if y >= 0 => positive branch
                ("LOAD", 0, 10),
                ("AND", 0, 3),
                ("CMP", 0, 5),
                ("JZ", pos_label),

                # negative branch: x -= (y>>i); y += (x_old>>i); z -= angle[i]
                ("LOAD", 0, 0), ("STORE", 0, 24),
                ("LOAD", 0, 2), ("STORE", 0, 26),
                ("LOAD", 0, 4), ("STORE", 0, 28),

                ("LOAD", 0, 6), ("STORE", 0, 18),
                ("LOAD", 0, 8), ("STORE", 0, 20),
                ("LOAD", 0, 10), ("STORE", 0, 22),
            ])
            for _ in range(i):
                rows.append(("CALL", "atan_sar1_temp"))
            rows.extend([
                ("LOAD", 0, 0), ("LOAD", 1, 18), ("SUB", 0, 1), ("STORE", 0, 0),
                ("LOAD", 0, 2), ("LOAD", 1, 20), ("SUBB", 0, 1), ("STORE", 0, 2),
                ("LOAD", 0, 4), ("LOAD", 1, 22), ("SUBB", 0, 1), ("STORE", 0, 4),

                ("LOAD", 0, 24), ("STORE", 0, 18),
                ("LOAD", 0, 26), ("STORE", 0, 20),
                ("LOAD", 0, 28), ("STORE", 0, 22),
            ])
            for _ in range(i):
                rows.append(("CALL", "atan_sar1_temp"))
            rows.extend([
                ("LOAD", 0, 6), ("LOAD", 1, 18), ("ADD", 0, 1), ("STORE", 0, 6),
                ("LOAD", 0, 8), ("LOAD", 1, 20), ("ADDC", 0, 1), ("STORE", 0, 8),
                ("LOAD", 0, 10), ("LOAD", 1, 22), ("ADDC", 0, 1), ("STORE", 0, 10),

                ("LOAD", 0, 12), ("LOAD", 1, a0), ("SUB", 0, 1), ("STORE", 0, 12),
                ("LOAD", 0, 14), ("LOAD", 1, a1), ("SUBB", 0, 1), ("STORE", 0, 14),
                ("LOAD", 0, 16), ("LOAD", 1, a2), ("SUBB", 0, 1), ("STORE", 0, 16),
                ("JMP", end_label),

                ("LABEL", pos_label),
                # positive branch: x += (y>>i); y -= (x_old>>i); z += angle[i]
                ("LOAD", 0, 0), ("STORE", 0, 24),
                ("LOAD", 0, 2), ("STORE", 0, 26),
                ("LOAD", 0, 4), ("STORE", 0, 28),

                ("LOAD", 0, 6), ("STORE", 0, 18),
                ("LOAD", 0, 8), ("STORE", 0, 20),
                ("LOAD", 0, 10), ("STORE", 0, 22),
            ])
            for _ in range(i):
                rows.append(("CALL", "atan_sar1_temp"))
            rows.extend([
                ("LOAD", 0, 0), ("LOAD", 1, 18), ("ADD", 0, 1), ("STORE", 0, 0),
                ("LOAD", 0, 2), ("LOAD", 1, 20), ("ADDC", 0, 1), ("STORE", 0, 2),
                ("LOAD", 0, 4), ("LOAD", 1, 22), ("ADDC", 0, 1), ("STORE", 0, 4),

                ("LOAD", 0, 24), ("STORE", 0, 18),
                ("LOAD", 0, 26), ("STORE", 0, 20),
                ("LOAD", 0, 28), ("STORE", 0, 22),
            ])
            for _ in range(i):
                rows.append(("CALL", "atan_sar1_temp"))
            rows.extend([
                ("LOAD", 0, 6), ("LOAD", 1, 18), ("SUB", 0, 1), ("STORE", 0, 6),
                ("LOAD", 0, 8), ("LOAD", 1, 20), ("SUBB", 0, 1), ("STORE", 0, 8),
                ("LOAD", 0, 10), ("LOAD", 1, 22), ("SUBB", 0, 1), ("STORE", 0, 10),

                ("LOAD", 0, 12), ("LOAD", 1, a0), ("ADD", 0, 1), ("STORE", 0, 12),
                ("LOAD", 0, 14), ("LOAD", 1, a1), ("ADDC", 0, 1), ("STORE", 0, 14),
                ("LOAD", 0, 16), ("LOAD", 1, a2), ("ADDC", 0, 1), ("STORE", 0, 16),

                ("LABEL", end_label),
            ])

        rows.extend([
            ("LABEL", "atan_done"),
            ("HALT",),

            # temp (18..22) arithmetic right shift by 1 (signed 48-bit)
            ("LABEL", "atan_sar1_temp"),
            ("LOAD", 0, 22), ("MOV", 1, 0), ("SAR", 0, 6), ("STORE", 0, 22), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 20), ("MOV", 1, 0), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "atan_sar1_no_or1"), ("OR", 0, 3), ("LABEL", "atan_sar1_no_or1"), ("STORE", 0, 20), ("MOV", 2, 1), ("AND", 2, 6),
            ("LOAD", 0, 18), ("SHR", 0, 6), ("CMP", 2, 5), ("JZ", "atan_sar1_no_or0"), ("OR", 0, 3), ("LABEL", "atan_sar1_no_or0"), ("STORE", 0, 18),
            ("RET",),
        ])
        return rows
