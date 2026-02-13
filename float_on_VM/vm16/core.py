from __future__ import annotations

from typing import List

from .assembler import Assembler
from .isa import ISA
from .utils import i16, u16


class VM16:
    REG_COUNT = 8
    SP_REG = 7
    RAM_SIZE = 512
    ROM_MAX_WORDS = 16384
    STACK_TOP = 0x01FE

    def __init__(self) -> None:
        self.reg: List[int] = [0] * self.REG_COUNT
        self.ram = bytearray(self.RAM_SIZE)
        self.rom_words: List[int] = []
        self.pc = 0
        self.halted = False
        self.steps = 0
        self.flags = {"Z": False, "N": False, "P": False, "C": False}
        self.reg[self.SP_REG] = self.STACK_TOP

    def reset(self) -> None:
        self.reg = [0] * self.REG_COUNT
        self.ram = bytearray(self.RAM_SIZE)
        self.rom_words = []
        self.pc = 0
        self.halted = False
        self.steps = 0
        self.flags = {"Z": False, "N": False, "P": False, "C": False}
        self.reg[self.SP_REG] = self.STACK_TOP

    def load_rom_bin(self, rom_bin: bytes) -> None:
        self.rom_words = Assembler.bin_to_words(rom_bin)
        if len(self.rom_words) > self.ROM_MAX_WORDS:
            raise ValueError(f"ROM too large: {len(self.rom_words)} words (max {self.ROM_MAX_WORDS})")
        self.pc = 0
        self.halted = False
        self.steps = 0
        self.reg[self.SP_REG] = self.STACK_TOP

    def _read_word(self, addr: int) -> int:
        a = addr % self.RAM_SIZE
        b = (a + 1) % self.RAM_SIZE
        return self.ram[a] | (self.ram[b] << 8)

    def _write_word(self, addr: int, value: int) -> None:
        a = addr % self.RAM_SIZE
        b = (a + 1) % self.RAM_SIZE
        v = u16(value)
        self.ram[a] = v & 0xFF
        self.ram[b] = (v >> 8) & 0xFF

    def _stack_push_word(self, value: int) -> None:
        sp = u16(self.reg[self.SP_REG] - 2)
        if sp < 0x40:
            raise RuntimeError("stack overflow")
        self._write_word(sp, value)
        self.reg[self.SP_REG] = sp

    def _stack_pop_word(self) -> int:
        sp = self.reg[self.SP_REG]
        if sp >= self.STACK_TOP:
            raise RuntimeError("stack underflow")
        out = self._read_word(sp)
        self.reg[self.SP_REG] = u16(sp + 2)
        return out

    def mem_set_i16(self, addr: int, value: int) -> None:
        self._write_word(addr, value)

    def mem_get_i16(self, addr: int) -> int:
        return i16(self._read_word(addr))

    def _set_cmp_flags(self, diff_i16: int) -> None:
        self.flags["Z"] = diff_i16 == 0
        self.flags["N"] = diff_i16 < 0
        self.flags["P"] = diff_i16 > 0

    def _r(self, idx: int) -> int:
        return self.reg[idx]

    def _set_r(self, idx: int, value: int) -> None:
        self.reg[idx] = u16(value)

    def step(self) -> None:
        if self.halted:
            return
        if self.pc < 0 or self.pc >= len(self.rom_words):
            raise RuntimeError(f"PC out of ROM range: {self.pc}")
        word = self.rom_words[self.pc]
        self.pc += 1
        self.steps += 1
        op, ra, rb, imm = ISA.decode(word)

        def fetch_ext_addr() -> int:
            if self.pc < 0 or self.pc >= len(self.rom_words):
                raise RuntimeError("missing extended address word")
            addr = self.rom_words[self.pc] & 0xFFFF
            self.pc += 1
            return addr

        if op == "HALT":
            self.halted = True
        elif op == "MOVI":
            self._set_r(ra, imm)
        elif op == "MOV":
            self._set_r(ra, self._r(rb))
        elif op == "LOAD":
            self._set_r(ra, self._read_word(imm))
        elif op == "STORE":
            self._write_word(imm, self._r(ra))
        elif op == "ADD":
            a = self._r(ra)
            b = self._r(rb)
            s = a + b
            self.flags["C"] = s > 0xFFFF
            self._set_r(ra, s)
        elif op == "ADDC":
            a = self._r(ra)
            b = self._r(rb)
            c = 1 if self.flags["C"] else 0
            s = a + b + c
            self.flags["C"] = s > 0xFFFF
            self._set_r(ra, s)
        elif op == "SUB":
            a = self._r(ra)
            b = self._r(rb)
            self.flags["C"] = a >= b
            self._set_r(ra, a - b)
        elif op == "SUBB":
            a = self._r(ra)
            b = self._r(rb)
            borrow = 0 if self.flags["C"] else 1
            d = a - b - borrow
            self.flags["C"] = d >= 0
            self._set_r(ra, d)
        elif op == "MUL":
            self._set_r(ra, i16(self._r(ra)) * i16(self._r(rb)))
        elif op == "UMUL":
            a = self._r(ra)
            b = self._r(rb)
            p = a * b
            self._set_r(ra, p & 0xFFFF)
            self._set_r(rb, (p >> 16) & 0xFFFF)
        elif op == "DIV":
            rhs = i16(self._r(rb))
            if rhs == 0:
                raise ZeroDivisionError("VM DIV by zero")
            lhs = i16(self._r(ra))
            self._set_r(ra, int(lhs / rhs))
        elif op == "CMP":
            self._set_cmp_flags(i16(self._r(ra)) - i16(self._r(rb)))
        elif op == "JMP":
            self.pc = fetch_ext_addr()
        elif op == "JZ":
            addr = fetch_ext_addr()
            if self.flags["Z"]:
                self.pc = addr
        elif op == "JNZ":
            addr = fetch_ext_addr()
            if not self.flags["Z"]:
                self.pc = addr
        elif op == "JN":
            addr = fetch_ext_addr()
            if self.flags["N"]:
                self.pc = addr
        elif op == "JP":
            addr = fetch_ext_addr()
            if self.flags["P"]:
                self.pc = addr
        elif op == "JC":
            addr = fetch_ext_addr()
            if self.flags["C"]:
                self.pc = addr
        elif op == "JNC":
            addr = fetch_ext_addr()
            if not self.flags["C"]:
                self.pc = addr
        elif op == "PUSH":
            self._stack_push_word(self._r(ra))
        elif op == "POP":
            self._set_r(ra, self._stack_pop_word())
        elif op == "CALL":
            addr = fetch_ext_addr()
            self._stack_push_word(self.pc)
            self.pc = addr
        elif op == "RET":
            self.pc = self._stack_pop_word()
        elif op == "AND":
            self._set_r(ra, self._r(ra) & self._r(rb))
        elif op == "OR":
            self._set_r(ra, self._r(ra) | self._r(rb))
        elif op == "XOR":
            self._set_r(ra, self._r(ra) ^ self._r(rb))
        elif op == "SHL":
            cnt = self._r(rb) & 0xF
            a = self._r(ra)
            if cnt == 0:
                self.flags["C"] = False
                self._set_r(ra, a)
            else:
                self.flags["C"] = ((a >> (16 - cnt)) & 0x1) != 0
                self._set_r(ra, a << cnt)
        elif op == "SHR":
            cnt = self._r(rb) & 0xF
            a = self._r(ra)
            if cnt == 0:
                self.flags["C"] = False
                self._set_r(ra, a)
            else:
                self.flags["C"] = ((a >> (cnt - 1)) & 0x1) != 0
                self._set_r(ra, a >> cnt)
        elif op == "SAR":
            cnt = self._r(rb) & 0xF
            a = i16(self._r(ra))
            if cnt == 0:
                self.flags["C"] = False
                self._set_r(ra, a)
            else:
                raw = self._r(ra)
                self.flags["C"] = ((raw >> (cnt - 1)) & 0x1) != 0
                self._set_r(ra, a >> cnt)
        else:
            raise ValueError(f"unknown decoded op: {op}")

    def run(self, max_steps: int = 100_000) -> None:
        while not self.halted and self.steps < max_steps:
            self.step()
        if not self.halted:
            raise RuntimeError(f"Program did not halt within {max_steps} steps")
