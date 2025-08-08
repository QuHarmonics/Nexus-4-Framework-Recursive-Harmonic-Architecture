import numpy as np
import matplotlib.pyplot as plt

# === Generate Simulated SHA-256 Register Values (Use real data if available) ===
rng = np.random.default_rng()
states = rng.integers(1_000_000_000, 4_000_000_000, size=(64, 8), dtype=np.int64)
# shape: (rounds, registers) → 64 rounds × 8 registers (A–H)

# === Harmonic Convergence Finder ===
def find_harmonic_intersections(states, tolerance=10_000_000):
    """
    Identifies rounds where multiple registers converge within a tolerance.

    Parameters:
        states (ndarray): Shape (rounds, 8), register values per round.
        tolerance (int): Maximum allowed value difference to consider harmonic.

    Returns:
        List of dicts: Each with round, converging indices, and mean value.
    """
    intersections = []
    for i, row in enumerate(states):
        diffs = np.abs(row[:, None] - row)
        np.fill_diagonal(diffs, np.inf)  # ignore self-comparison
        harmonics = np.where(diffs < tolerance)
        grouped = set()
        for a, b in zip(*harmonics):
            if a < b:  # avoid duplicates
                grouped.add(a)
                grouped.add(b)
        if grouped:
            values = [row[j] for j in grouped]
            intersections.append({
                "round": i,
                "registers": list(grouped),
                "mean_value": int(np.mean(values))
            })
    return intersections

# === Run Harmonic Finder ===
harmonics = find_harmonic_intersections(states, tolerance=10_000_000)

# === Plot Round States and Mark Harmonics ===
plt.figure(figsize=(14, 6))
for i in range(8):
    plt.plot(states[:, i], label=chr(65 + i))  # A–H

for h in harmonics:
    plt.axvline(h["round"], color='gray', linestyle='--', alpha=0.5)

plt.title("SHA-256 Register Round States with Harmonic Intersections")
plt.xlabel("Round")
plt.ylabel("Register Value")
plt.legend(title="Registers (A–H)")
plt.grid(True)
plt.tight_layout()
plt.show()

# === Preview Harmonic Events ===
for h in harmonics[:5]:
    print(f"Round {h['round']}: Registers {h['registers']} aligned near {h['mean_value']}")
