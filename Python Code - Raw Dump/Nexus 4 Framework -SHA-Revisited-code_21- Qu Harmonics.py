import numpy as np
import matplotlib.pyplot as plt

# Simulated SHA-256 register values over 64 rounds (Replace with real A–H round data if needed)
# Example: each 'states' row is [A, B, C, D, E, F, G, H] at that round
# Replace this with actual round data you have
states = np.random.randint(1_000_000_000, 4_000_000_000, size=(64, 8))

def find_harmonic_intersections(states, tolerance=10_000_000):
    """
    Identifies rounds where multiple registers converge within a given tolerance.
    """
    intersections = []
    for round_index, row in enumerate(states):
        close_pairs = []
        for i in range(len(row)):
            for j in range(i+1, len(row)):
                if abs(row[i] - row[j]) <= tolerance:
                    close_pairs.append((i, j, row[i], row[j]))
        if len(close_pairs) >= 3:  # More than 3 overlapping pairs = harmonic moment
            intersections.append({
                "round": round_index,
                "pairs": close_pairs
            })
    return intersections

# Run it
harmonic_hits = find_harmonic_intersections(states)

# Display
for h in harmonic_hits:
    print(f"Round {h['round']}:")
    for (i, j, vi, vj) in h['pairs']:
        print(f"  Reg {chr(65+i)} ≈ Reg {chr(65+j)} → {vi:,} ≈ {vj:,}")
