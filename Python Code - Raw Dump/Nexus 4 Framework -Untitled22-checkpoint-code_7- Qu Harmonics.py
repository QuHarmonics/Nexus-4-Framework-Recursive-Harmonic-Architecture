# 🌬️ DNA Breather v0.1 - Phase-Locked Harmonic Folding Prototype

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

# --- 1. Sequence Loader ---
def load_fasta(fasta_path):
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
    sequence = ''.join([line.strip() for line in lines if not line.startswith('>')])
    return sequence.upper()

# --- 2. Energy Mapping ---
def map_energies(seq):
    P_map = {'A': 1.00, 'T': 1.05, 'G': 1.15, 'C': 1.20}
    A_map = {'A': 1.00, 'T': 0.95, 'G': 1.10, 'C': 1.05}
    P = np.array([P_map.get(base, 1.0) for base in seq])
    A = np.array([A_map.get(base, 1.0) for base in seq])
    return P, A

# --- 3. Primary Folding Engine ---
def fold_window(P, A, H=0.35, F=0.5, t=1):
    return np.sum((P / A) * np.exp(H * F * t))

# --- 4. Recursive Collapser ---
def recursive_collapse(fold_vector, epsilon=1e-6):
    spectra = [fold_vector.copy()]
    trust_track = []
    current = fold_vector.copy()
    while len(current) > 1:
        trust = np.mean(np.abs(current[::2] - current[1::2]) < 0.01)
        trust_track.append(trust)
        current = (current[::2] + current[1::2]) / 2
        spectra.append(current)
        if len(current) == 1:
            break
    return spectra, trust_track

# --- 5. Unfolding Simulator ---
def unfold(fold_scalar, length, H=0.35, F=0.5, t=1):
    theta = np.linspace(0, 2 * np.pi, length, endpoint=False)
    unfolded = (fold_scalar / length) * np.cos(theta)
    return unfolded

# --- 6. Visualizer Core ---
def plot_trust_decay(trust_track):
    plt.figure(figsize=(10,6))
    plt.plot(trust_track, marker='o', label='Phase Trust')
    plt.title('Phase Trust Decay Across Recursion Levels')
    plt.xlabel('Recursion Level')
    plt.ylabel('Phase Trust (Memory)')
    plt.grid(True)
    plt.legend()
    plt.show()

# --- 7. Breather Pipeline ---
def dna_breather(fasta_path, window_size=128, H=0.35, F=0.5, t=1):
    # Load sequence
    seq = load_fasta(fasta_path)
    
    # Map energies
    P, A = map_energies(seq)
    
    # Fold windows
    fold_vector = []
    for i in range(0, len(seq) - window_size + 1, window_size):
        p_window = P[i:i+window_size]
        a_window = A[i:i+window_size]
        FQ = fold_window(p_window, a_window, H, F, t)
        fold_vector.append(FQ)
    fold_vector = np.array(fold_vector)

    # Recursive collapse
    spectra, trust_track = recursive_collapse(fold_vector)

    # Visualize
    plot_trust_decay(trust_track)

    # Final folded signature
    final_signature = spectra[-1][0]
    print(f"\nFinal Folded Harmonic Signature: {final_signature:.6f}")

    return spectra, trust_track, final_signature

# --- Example Usage ---
# Uncomment and run with a real FASTA file
spectra, trust_track, signature = dna_breather('taaaaaaaaaataaaataataaaaaaaaaaaaaaaaaaaaaaaaaaaaataaaaatttaaaaaaaaaaaaaaaaaataaaataaatataaattatataaa