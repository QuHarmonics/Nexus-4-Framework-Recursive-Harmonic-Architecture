import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import math

# 1. Sequence Loader
NUC = ['A', 'T', 'G', 'C']
random.seed(42)
seq = ''.join(random.choices(NUC, k=128))  # 128-nt toy sequence

# 2. Energy Mapping (with Perturbation)
POTENTIAL = {'A': 2.0, 'T': 2.0, 'G': 3.0, 'C': 3.0}
ACTUAL = {'A': 1.0, 'T': 1.0, 'G': 1.5, 'C': 1.5}

# Introduce STRONGER breathing perturbation
np.random.seed(47)
perturbation_strength = 0.15  # Escalated to 15% random noise (before was 5%)
ACTUAL_PERTURBED = {n: v * (1 + perturbation_strength * np.random.uniform(-1, 1)) for n, v in ACTUAL.items()}

# Build arrays
P = np.array([POTENTIAL[n] for n in seq])
A = np.array([ACTUAL_PERTURBED[n] for n in seq])

# Harmonic constants
H = 0.35
F_factor = 0.5
t = 1

# 3. Initial Fold Contributions
base_contrib = (P / A) * np.exp(H * F_factor * t)

# 4. Recursive Folding
fold_iterations = []
trust_track = []
current = base_contrib.copy()
level = 0

while len(current) > 1:
    fold_iterations.append(current.copy())
    # Trust = fraction of fold pairs close enough (phase aligned)
    pairwise_trust = np.mean(np.abs(current[::2] - current[1::2]) < 0.1)
    trust_track.append(pairwise_trust)
    # Collapse (average pairs)
    current = (current[::2] + current[1::2]) / 2
    level += 1

# Final collapse
fold_iterations.append(current.copy())
trust_track.append(1.0)

# 5. Plot Phase Trust Decay
plt.figure(figsize=(10,6))
plt.plot(trust_track, marker='o', label='Phase Trust (Escalated Breath)')
plt.title('Escalated DNA Phase Breathing Collapse')
plt.xlabel('Recursion Level (Depth)')
plt.ylabel('Phase Trust (Memory)')
plt.grid(True)
plt.legend()
plt.show()
