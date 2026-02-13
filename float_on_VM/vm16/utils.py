from __future__ import annotations


def u16(v: int) -> int:
    return v & 0xFFFF


def i16(v: int) -> int:
    v &= 0xFFFF
    return v - 0x10000 if v & 0x8000 else v

