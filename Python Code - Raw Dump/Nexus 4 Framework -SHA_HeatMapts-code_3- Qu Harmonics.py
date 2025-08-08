#!/usr/bin/env python3
"""
bit_level_driver_analysis.py

Compute per-bit Shannon entropy for “folded” vs. “unfolded” SHA-derived residue states,
and plot a temporal bit-flip heatmap to identify transition spikes.

Requirements:
    - numpy
    - matplotlib

Usage:
    python bit_level_driver_analysis.py
"""

import numpy as np
import matplotlib.pyplot as plt

# Example data: replace these with your actual arrays
# residue: sequence of integer residue values (e.g. mod 64) over time
# folded_mask: boolean array same length as residue, True when state is in “folded” regime
# bit_array: 2D array shape (T, 8) of the 8 bits of residue[t], from most significant (bit7) to LSB (bit0)
np.random.seed(0)
T = 1000
residue = np.random.randint(0, 64, size=T)
# define “folded” as residue < 32 for example
folded_mask = residue < 32

# Build bit_array
bit_array = ((residue[:, None] & (1 << np.arange(7, -1, -1))) > 0).astype(int)

def shannon_entropy(p: np.ndarray) -> float:
    """Compute Shannon entropy (base-2) for a Bernoulli distribution with probability p."""
    if p == 0 or p == 1:
        return 0.0
    return -p*np.log2(p) - (1-p)*np.log2(1-p)

# 1) Compute per-bit probabilities and entropy in folded vs. unfolded regimes
folded_bits   = bit_array[folded_mask]
unfolded_bits = bit_array[~folded_mask]
entropies = {'folded': [], 'unfolded': []}

for regime, bits in [('folded', folded_bits), ('unfolded', unfolded_bits)]:
    Hs = []
    # for each of the 8 bit positions
    for b in range(8):
        p1 = bits[:, b].mean()
        Hs.append(shannon_entropy(p1))
    entropies[regime] = Hs

# Plot per-bit entropy comparison
bits = np.arange(7, -1, -1)  # bit indices MSB to LSB
plt.figure()
plt.plot(bits, entropies['folded'], marker='o', label='Folded')
plt.plot(bits, entropies['unfolded'], marker='s', label='Unfolded')
plt.xlabel('Bit position (7=MSB, 0=LSB)')
plt.ylabel('Shannon Entropy (bits)')
plt.title('Per-Bit Entropy: Folded vs. Unfolded States')
plt.legend()
plt.xticks(bits)
plt.gca().invert_xaxis()
plt.tight_layout()
plt.show()

# 2) Temporal bit-flip heatmap
# Compute flips: 1 where bit changes between t and t+1
flip_matrix = np.abs(np.diff(bit_array, axis=0))

plt.figure()
plt.imshow(flip_matrix.T, aspect='auto', origin='lower')
plt.colorbar(label='Flip (1=yes,0=no)')
plt.xlabel('Time step')
plt.ylabel('Bit position (7→0)')
plt.yticks(np.arange(8), [f'bit{b}' for b in range(7, -1, -1)])
plt.title('Temporal Bit-Flip Heatmap')
plt.tight_layout()
plt.show()
