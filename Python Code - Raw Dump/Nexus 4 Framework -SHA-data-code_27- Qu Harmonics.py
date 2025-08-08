#!/usr/bin/env python3
"""
summarize_cycle_alignment.py

Processes cycle‐vs‐round‐constant distance & Hamming tables across multiple moduli,
and produces a concise summary of alignment patterns.

1) For each modulus M:
   • Count how many states have exact matches (distance == 0).
   • Tally which K_j indices appear most frequently as closest.
2) Identify states that are exact matches at all moduli.
3) Identify K_j indices that “dominate” (closest most often) across all moduli.
4) Print a human‐readable report.

Usage:
    python summarize_cycle_alignment.py
"""

from collections import Counter, defaultdict
from typing import List, Dict, Tuple

# Your cycle states
STATES: List[int] = [21, 8, 154, 65, 253, 155, 94, 203, 29, 205, 209, 153, 88]

# Full list of SHA-256 round constants K_j (32-bit)
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

def compute_closest(states: List[int], M: int) -> Dict[int, Tuple[int,int]]:
    """
    For each state s, find (j, distance) where K_j mod M is closest.
    Returns dict s -> (j, dist).
    """
    K_mod = [k % M for k in K]
    result = {}
    for s in states:
        best_j, best_d = None, M
        for j, kj in enumerate(K_mod):
            delta = (s - kj) % M
            dist = min(delta, M - delta)
            if dist < best_d:
                best_d, best_j = dist, j
        result[s] = (best_j, best_d)
    return result

def summarize():
    moduli = [256, 128, 64, 32, 16, 8]
    # store per-modulus info
    exact_counts = {}
    j_tallies = {}
    exact_states_per_mod = {}

    for M in moduli:
        closest = compute_closest(STATES, M)
        # count exact matches
        exact_states = [s for s,(j,d) in closest.items() if d == 0]
        exact_counts[M] = len(exact_states)
        exact_states_per_mod[M] = set(exact_states)
        # tally which j's occur
        tally = Counter(j for j,(j,d) in closest.items())
        j_tallies[M] = tally

    # states exact in all moduli
    always_exact = set(STATES)
    for sset in exact_states_per_mod.values():
        always_exact &= sset

    # overall j frequency across moduli
    overall_j = Counter()
    for tally in j_tallies.values():
        overall_j.update(dict(tally))

    # produce report
    print("\n=== Exact‐Match Counts per Modulus ===")
    for M in moduli:
        print(f" M = {M:3d} → {exact_counts[M]:2d} / {len(STATES)} states exact")

    print("\n=== States Exact at Every Modulus ===")
    if always_exact:
        print(" ", sorted(always_exact))
    else:
        print("  (none)")

    print("\n=== Top 5 Round Constants by Total Closest Frequency ===")
    for j, freq in overall_j.most_common(5):
        print(f"  K_{j:02d}: {freq} times")

    print("\n=== Most Frequent Closest j per Modulus ===")
    for M in moduli:
        most_j, count = j_tallies[M].most_common(1)[0]
        print(f" M = {M:3d} → K_{most_j:02d} ({count} states)")

if __name__ == "__main__":
    summarize()
