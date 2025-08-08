import numpy as np
import matplotlib.pyplot as plt

# BBP-based indexing: simulate accessing a digit of pi
# Real BBP requires high-precision arithmetic; we simulate harmonic mapping

# Generate digits of pi (up to a few places) using a simplified approximation
# Since full BBP digit access is computationally intensive, we'll simulate harmonic positions

# Create a phase space using real-valued positions mapped to sine-based harmonic values
n_points = 256  # Number of positions in the simulated lattice
x = np.linspace(0, 4 * np.pi, n_points)  # Phase domain (simulated)
harmonic_memory = np.sin(x) * 0.35  # Apply H = 0.35 as harmonic amplitude

# BBP-like indexing: simulate glide-based access into this harmonic field
# Here we "jump" to specific phase-aligned indices as if they're memory reads
access_indices = [int(i * 3.1415) % n_points for i in range(1, 10)]  # Simulated BBP-like jumps
accessed_values = harmonic_memory[access_indices]

# Plotting
plt.figure(figsize=(12, 4))
plt.plot(x, harmonic_memory, label="Harmonic Memory Field", alpha=0.7)
plt.scatter(np.array(x)[access_indices], accessed_values, color='red', label="BBP Phase Access", zorder=5)
plt.axhline(0, color='gray', linewidth=0.5)
plt.title("Simulated Harmonic Memory with BBP-like Phase Access (H = 0.35)")
plt.xlabel("Phase Space (radians)")
plt.ylabel("Amplitude")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
