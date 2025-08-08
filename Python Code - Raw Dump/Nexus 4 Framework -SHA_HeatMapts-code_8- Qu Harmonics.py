#!/usr/bin/env python3
"""
msb_flip_constant_correlation.py

For each flip in bit7 (and bit6), record the residue before the flip,
and list which SHA-256 constants (mod 64) match that residue.

Usage:
    python msb_flip_constant_correlation.py
"""

import numpy as np
from collections import Counter

# --- Replace these with your actual data arrays ---
# residue: shape (T,) integers in [0,255] or [0,63]
# bit_array: shape (T,8), bits of residue (bit_array[:,0]=MSB)
# K: list of 64 SHA-256 constants

# For demo purposes we re-use the random example from before:
np.random.seed(0)
T = 1000
residue = np.random.randint(0, 64, size=T)
bit_array = ((residue[:,None] & (1 << np.arange(7,-1,-1))) > 0).astype(int)
K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]

# Precompute K mod 64
K_mod64 = [k & 0x3F for k in K]  # mask lower 6 bits

# Find bit 7 flips
flip_times = np.where(bit_array[:-1,0] != bit_array[1:,0])[0]  # 0 = MSB

# Correlate residues to K indices
correlations = []
for t in flip_times:
    r = residue[t]
    # find all j where K_mod64[j] == r
    js = [j for j, km in enumerate(K_mod64) if km == r]
    correlations.extend(js)

# Summarize most common constants at MSB flip times
counter = Counter(correlations)
print("MSB Flip → Top SHA-256 K_j (mod 64) correlations:")
for j, cnt in counter.most_common(10):
    print(f"  K_{j:02d} (0x{K[j]:08x}) → {cnt} flips")

# (Optional) Repeat for bit6 by using bit_array[:,1]
