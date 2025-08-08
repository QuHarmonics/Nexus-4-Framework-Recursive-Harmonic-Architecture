import numpy as np
import matplotlib.pyplot as plt

# -------------------------------
# 1. Universal Harmonic Resonance (Mark 1)
# -------------------------------
n = 5  # number of segments

# Randomly generate potential energies (P) and actual energies (A)
P = np.random.uniform(1, 10, size=n)
A = np.random.uniform(1, 10, size=n)

# Compute harmonic resonance H
H_value = np.sum(P) / np.sum(A)

print("Universal Harmonic Resonance (Mark 1)")
print("Potential energies P:", P)
print("Actual energies A:", A)
print("Computed H =", H_value)
print()

# -------------------------------
# 2. Kulik Recursive Reflection (KRR)
# -------------------------------
# Formula: R(t) = R0 * exp(H * F * t)
R0 = 1.0         # initial reflection state
F_factor = 0.5   # folding factor (example value)

# Define time values (e.g., t from 0 to 10)
t_values = np.linspace(0, 10, 100)
R_t = R0 * np.exp(H_value * F_factor * t_values)

# Plot R(t) over time
plt.figure(figsize=(8, 5))
plt.plot(t_values, R_t, label=r'$R(t) = R_0 e^{(H\cdot F\cdot t)}$')
plt.title("Kulik Recursive Reflection (KRR)")
plt.xlabel("Time t")
plt.ylabel("Reflection State R(t)")
plt.legend()
plt.grid(True)
plt.show()

# -------------------------------
# 3. Dynamic Noise Filtering (DNF)
# -------------------------------
# Simulate n random noise differences ΔN_i
k_constant = 0.1  # noise sensitivity factor
Delta_N = np.random.uniform(-5, 5, size=n)

# Compute the filtered noise output N(t)
N_value = np.sum(Delta_N / (1 + k_constant * np.abs(Delta_N)))

print("Dynamic Noise Filtering (DNF)")
print("Noise differences ΔN:", Delta_N)
print("Filtered Noise Output N(t) =", N_value)
