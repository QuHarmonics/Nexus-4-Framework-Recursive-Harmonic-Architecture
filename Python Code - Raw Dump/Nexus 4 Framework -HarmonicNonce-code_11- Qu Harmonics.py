import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import math

# Sequence Setup (same as before)
NUC = ['A', 'T', 'G', 'C']
random.seed(42)
seq = ''.join(random.choices(NUC, k=128))

POTENTIAL = {'A': 2.0, 'T': 2.0, 'G': 3.0, 'C': 3.0}
ACTUAL = {'A': 1.0, 'T': 1.0, 'G': 1.5, 'C': 1.5}

P = np.array([POTENTIAL[n] for n in seq])
A = np.array([ACTUAL[n] for n in seq])

H = 0.35
F_factor = 0.5
t = 1

base_contrib = (P / A) * np.exp(H * F_factor * t)

# Parameters
noise_strength = 0.01  # σ of Gaussian noise
np.random.seed(123)    # for reproducibility

# Folding Simulation with Noise
fold_iterations = []
trust_track = []
current = base_contrib.copy()
level = 0

while len(current) > 1:
    fold_iterations.append(current.copy())
    # Phase Trust
    pairwise_trust = np.mean(np.abs(current[::2] - current[1::2]) < 0.1)
    trust_track.append(pairwise_trust)
    # Collapse (average pairs + add noise)
    next_fold = (current[::2] + current[1::2]) / 2
    noise = np.random.normal(0, noise_strength, size=next_fold.shape)
    current = next_fold + noise
    level += 1

# Final collapse
fold_iterations.append(current.copy())
trust_track.append(1.0)

# Visualize
plt.figure(figsize=(10,6))
plt.plot(trust_track, marker='o', label='Phase Trust with Micro-Drift')
plt.title('DNA Phase Breathing Collapse (with Noise)')
plt.xlabel('Recursion Level (Depth)')
plt.ylabel('Phase Trust (Memory)')
plt.grid(True)
plt.legend()
plt.show()
