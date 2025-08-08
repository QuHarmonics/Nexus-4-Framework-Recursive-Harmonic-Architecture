#!/usr/bin/env python3
"""
cycle_full_analysis.py

Comprehensive analysis of a 13‐state byte cycle:

 1) Alignment with SHA‐256 round constants at moduli:
    256, 128, 64, 32, 16, 8, and 2^16 (65536).
 2) Nibble‐chain decomposition:
    • Split each byte into high/low 4‐bit nibbles.
    • Sum each nibble pair.
    • Group nibble‐sums in adjacent pairs → digital root.
 3) 16‐bit “echo” check: treat every two bytes as one 16‐bit word,
    apply the same nibble logic to confirm echoing.

Usage:
    python cycle_full_analysis.py
"""

from typing import List, Dict, Tuple
from collections import Counter
import sys

# Your 13‐state cycle (8‐bit values)
STATES: List[int] = [
    0x14, 0x15, 0x1A, 0x0F, 0x2E, 0x39, 0x4B,
    0x5C, 0x6D, 0x7E, 0x83, 0x94, 0xA5
]

# SHA‐256 round constants K_j (32‐bit)
K: List[int] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]

def compute_min_distances(states: List[int], M: int) -> Dict[int, Tuple[int,int]]:
    """
    For each state s, find (j, distance) where
    distance = min((s - (K_j mod M)) mod M, (K_j mod M) - s mod M).
    """
    K_mod = [k % M for k in K]
    result: Dict[int, Tuple[int,int]] = {}
    for s in states:
        best_j, best_d = None, M
        for j, kj in enumerate(K_mod):
            delta = (s - kj) % M
            dist  = min(delta, M - delta)
            if dist < best_d:
                best_d, best_j = dist, j
        result[s] = (best_j, best_d)
    return result

def nibble_chain(states: List[int]) -> None:
    """
    Decompose each byte into nibbles, sum, pairwise digital‐root,
    and 16‐bit echo check.
    """
    def split_nibbles(b: int) -> Tuple[int,int]:
        return (b >> 4) & 0xF, b & 0xF

    def digital_root(n: int) -> int:
        while n >= 10:
            n = sum(int(d) for d in str(n))
        return n

    nib_sums = [sum(split_nibbles(b)) for b in states]
    chained  = [digital_root(nib_sums[i] + nib_sums[i+1]) 
                if i+1 < len(nib_sums) else digital_root(nib_sums[i])
                for i in range(0, len(nib_sums), 2)]

    print("\n--- Nibble‐Chain Analysis ---")
    print("States (hex):      ", [f"0x{b:02X}" for b in states])
    print("Nibble sums:       ", nib_sums)
    print("Pairwise droots:   ", chained)

    # 16‐bit echo check
    print("\n16-bit word echo check:")
    for i in range(0, len(states), 2):
        hi, lo = states[i], states[i+1] if i+1 < len(states) else 0
        hi_hi, hi_lo = split_nibbles(hi)
        lo_hi, lo_lo = split_nibbles(lo)
        sum1 = hi_hi + hi_lo
        sum2 = lo_hi + lo_lo
        print(f" Word {i//2}: [0x{hi:02X},0x{lo:02X}] → "
              f"nib sums [{sum1},{sum2}] → "
              f"droots [{digital_root(sum1)},{digital_root(sum2)}]")

def full_alignment_report():
    moduli = [256, 128, 64, 32, 16, 8, 2**16]
    print("--- Alignment with SHA-256 Constants ---")
    for M in moduli:
        results = compute_min_distances(STATES, M)
        exact = sum(1 for (_,d) in results.values() if d == 0)
        top5 = Counter(j for (j,_) in results.values()).most_common(5)
        print(f"\nModulus {M}: {exact}/{len(STATES)} exact matches")
        print(" Top-5 closest K_j:", top5)

if __name__ == "__main__":
    full_alignment_report()
    nibble_chain(STATES)
