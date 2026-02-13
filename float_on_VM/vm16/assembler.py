from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

from .isa import ISA

AsmRow = Tuple[object, ...]


class Assembler:
    @staticmethod
    def assemble_rows(rows: Sequence[AsmRow]) -> List[int]:
        labels: Dict[str, int] = {}
        pc = 0
        for row in rows:
            if row[0] == "LABEL":
                labels[str(row[1])] = pc
            else:
                op = str(row[0])
                if op in {"JMP", "JZ", "JNZ", "JN", "JP", "JC", "JNC", "CALL"}:
                    pc += 2  # extended jump/call: opcode word + 16-bit address word
                else:
                    pc += 1

        words: List[int] = []
        for row in rows:
            op = str(row[0])
            if op == "LABEL":
                continue
            a = row[1] if len(row) > 1 else None
            b = row[2] if len(row) > 2 else None
            if isinstance(a, str) and a in labels:
                a = labels[a]
            if isinstance(b, str) and b in labels:
                b = labels[b]

            if op in {"MOV", "ADD", "ADDC", "SUB", "SUBB", "MUL", "DIV", "CMP", "AND", "OR", "XOR", "SHL", "SHR", "SAR", "UMUL"}:
                words.append(ISA.encode_rr(op, int(a), int(b)))
            elif op in {"MOVI", "LOAD", "STORE"}:
                words.append(ISA.encode_ri8(op, int(a), int(b)))
            elif op in {"JMP", "JZ", "JNZ", "JN", "JP", "JC", "JNC", "CALL"}:
                # Extended form to support larger ROM spaces.
                words.append(ISA.encode_j11(op, 0))
                words.append(int(a) & 0xFFFF)
            elif op in {"PUSH", "POP"}:
                words.append(ISA.encode_r(op, int(a)))
            elif op in {"HALT", "RET"}:
                words.append(ISA.encode_0(op))
            else:
                raise ValueError(f"unsupported row op: {op}")
        return words

    @staticmethod
    def words_to_bin(words: Sequence[int]) -> bytes:
        out = bytearray()
        for w in words:
            out.append(w & 0xFF)
            out.append((w >> 8) & 0xFF)
        return bytes(out)

    @staticmethod
    def bin_to_words(data: bytes) -> List[int]:
        if len(data) % 2 != 0:
            raise ValueError("ROM binary length must be even")
        out: List[int] = []
        for i in range(0, len(data), 2):
            out.append(data[i] | (data[i + 1] << 8))
        return out

    @staticmethod
    def parse_asm_text(text: str) -> List[AsmRow]:
        rows: List[AsmRow] = []
        for raw in text.splitlines():
            line = raw.split(";", 1)[0].strip()
            if not line:
                continue
            if line.endswith(":"):
                rows.append(("LABEL", line[:-1].strip()))
                continue

            line = line.replace(",", " ")
            parts = [p for p in line.split() if p]
            op = parts[0].upper()

            def parse_operand(tok: str) -> object:
                t = tok.strip()
                if t.upper().startswith("R") and t[1:].isdigit():
                    return int(t[1:])
                if t.lstrip("-").isdigit():
                    return int(t)
                return t

            if len(parts) == 1:
                rows.append((op,))
            elif len(parts) == 2:
                rows.append((op, parse_operand(parts[1])))
            elif len(parts) == 3:
                rows.append((op, parse_operand(parts[1]), parse_operand(parts[2])))
            else:
                raise ValueError(f"invalid asm line: {raw}")
        return rows

    @staticmethod
    def disassemble_words(words: Sequence[int]) -> List[str]:
        out: List[str] = []
        i = 0
        while i < len(words):
            w = words[i]
            op, ra, rb, imm = ISA.decode(w)
            if op in {"MOV", "ADD", "ADDC", "SUB", "SUBB", "MUL", "DIV", "CMP", "AND", "OR", "XOR", "SHL", "SHR", "SAR", "UMUL"}:
                out.append(f"{op} R{ra}, R{rb}")
                i += 1
            elif op in {"MOVI", "LOAD", "STORE"}:
                out.append(f"{op} R{ra}, {imm}")
                i += 1
            elif op in {"JMP", "JZ", "JNZ", "JN", "JP", "JC", "JNC", "CALL"}:
                if i + 1 >= len(words):
                    out.append(f"{op} <missing_addr>")
                    i += 1
                else:
                    out.append(f"{op} {words[i+1]}")
                    i += 2
            elif op in {"PUSH", "POP"}:
                out.append(f"{op} R{ra}")
                i += 1
            elif op in {"HALT", "RET"}:
                out.append(op)
                i += 1
            else:
                out.append(f".word 0x{w:04X}")
                i += 1
        return out


def compile_firmware(rows: Sequence[AsmRow]) -> bytes:
    words = Assembler.assemble_rows(rows)
    return Assembler.words_to_bin(words)
