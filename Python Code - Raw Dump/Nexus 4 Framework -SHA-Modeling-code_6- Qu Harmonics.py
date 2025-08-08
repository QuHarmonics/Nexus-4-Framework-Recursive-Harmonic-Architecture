import numpy as np
import matplotlib.pyplot as plt

# -------------------------------
# Quantum Folding and Unfolding Simulation
# -------------------------------

# Given sample data from earlier (Universal Harmonic Resonance)
P = np.array([5.02095126, 8.37535711, 9.47325124, 4.65076768, 2.20552366])
A = np.array([6.57021632, 6.21162321, 6.00069682, 6.09087803, 7.18640782])
n = len(P)

# Previously computed H from Mark 1 simulation:
H_value = 0.9271994949802294

# Set parameters for quantum folding:
F_factor = 0.5   # Folding factor
t = 1.0          # Time (or recursive depth) parameter

# Compute the folded contributions for each segment:
FQ_components = (P / A) * np.exp(H_value * F_factor * t)
# Sum them to get the overall folded quantum state:
FQ = np.sum(FQ_components)

print("Quantum Folding Simulation:")
print("P (Potential energies):", P)
print("A (Actual energies):", A)
print("Folding factor F =", F_factor, ", Time t =", t)
print("Computed F(Q)_i for each segment:", FQ_components)
print("Overall Folded Quantum State F(Q):", FQ)
print()

# -------------------------------
# Quantum Unfolding Simulation
# -------------------------------
# Choose random phase angles for each segment between 0 and 2π:
theta = np.random.uniform(0, 2*np.pi, size=n)
# Set a small residual term zeta:
zeta = 0.1

# Compute the unfolded quantum state:
UQ = np.sum(FQ_components * np.cos(theta)) + zeta

print("Quantum Unfolding Simulation:")
print("Random phase angles theta (radians):", theta)
print("Computed Unfolded Quantum State U(Q):", UQ)
print()

# Optional: Visualizing the folded contributions and their cosine adjustments
plt.figure(figsize=(8,5))
plt.bar(np.arange(n), FQ_components, alpha=0.7, label="Folded Components F(Q)_i")
plt.plot(np.arange(n), FQ_components * np.cos(theta), 'ro', label="Cosine-adjusted components")
plt.title("Quantum Folding Components and Unfolding Adjustments")
plt.xlabel("Segment Index")
plt.ylabel("Value")
plt.legend()
plt.grid(True)
plt.show()
