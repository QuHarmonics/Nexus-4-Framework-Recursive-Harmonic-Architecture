import numpy as np
import matplotlib.pyplot as plt

# -------------------------------
# 1. Samson's Law (Feedback Stabilization)
# -------------------------------
k_const = 0.1        # Feedback constant
T_value = 2.0        # Constant time interval
n_segments = 5       # Number of segments

# Generate random ΔF values (force changes) for 5 segments
Delta_F = np.random.uniform(-10, 10, size=n_segments)
Delta_E = k_const * Delta_F
S_segments = Delta_E / T_value

print("Samson's Law Simulation:")
print("ΔF values:", Delta_F)
print("Computed ΔE values (k * ΔF):", Delta_E)
print("Stabilization rates S (ΔE/T):", S_segments)
print("Overall average S:", np.mean(S_segments))
print()

# -------------------------------
# 2. Multi-Dimensional Samson (MDS)
# -------------------------------
dims = 2           # Number of dimensions
n_segments = 5     # Segments per dimension
k_i = 0.1          # Feedback constant per segment

# Generate random ΔF and T values for each dimension
Delta_F_md = np.random.uniform(-10, 10, size=(dims, n_segments))
T_md = np.random.uniform(1, 5, size=(dims, n_segments))
Delta_E_md = k_i * Delta_F_md

# For each dimension, compute S_d = (sum(ΔE_i))/(sum(T_i))
S_d = np.sum(Delta_E_md, axis=1) / np.sum(T_md, axis=1)

print("Multi-Dimensional Samson (MDS) Simulation:")
print("ΔF for each dimension:\n", Delta_F_md)
print("T values for each dimension:\n", T_md)
print("Computed ΔE for each dimension:\n", Delta_E_md)
print("Stabilization rates S_d for each dimension:", S_d)
print()

# -------------------------------
# 3. Harmonic Memory Growth (HMG)
# -------------------------------
M0 = 1.0         # Initial memory capacity
alpha = 0.2      # Growth rate constant
C_const = 0.35   # Harmonic constant (target)
# Use previously computed H_value from Mark 1 (from your simulation) or sample value:
H_val = 0.9272

t_vals = np.linspace(0, 10, 100)
M_t = M0 * np.exp(alpha * (H_val - C_const) * t_vals)

plt.figure(figsize=(8,5))
plt.plot(t_vals, M_t, label="Harmonic Memory Growth M(t)")
plt.title("Harmonic Memory Growth (HMG)")
plt.xlabel("Time t")
plt.ylabel("Memory M(t)")
plt.legend()
plt.grid(True)
plt.show()

# -------------------------------
# 4. Quantum Jump Factor (QJF)
# -------------------------------
Q_factor = 0.05  # Weight adjustment for quantum transitions
Q_t = 1 + H_val * t_vals * Q_factor

plt.figure(figsize=(8,5))
plt.plot(t_vals, Q_t, label="Quantum Jump Factor Q(t)")
plt.title("Quantum Jump Factor (QJF)")
plt.xlabel("Time t")
plt.ylabel("Q(t)")
plt.legend()
plt.grid(True)
plt.show()

# -------------------------------
# 5. Task Distribution
# -------------------------------
n_nodes = 5
W = np.random.uniform(1, 10, size=n_nodes)      # Workload demands
C_nodes = np.random.uniform(1, 10, size=n_nodes)  # Node capacities

# Calculate task distribution for each node
task_distribution = (W * C_nodes) / np.sum(W * C_nodes)

print("Task Distribution Simulation:")
print("Workload demands W:", W)
print("Capacities C:", C_nodes)
print("Task distribution T(i):", task_distribution)
print("Sum of task distribution (should be 1):", np.sum(task_distribution))
