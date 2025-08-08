import numpy as np
import matplotlib.pyplot as plt

# Set random seed for reproducibility
np.random.seed(42)

# -------------------------------
# Simulation Parameters
# -------------------------------
time_steps = 20            # Number of discrete time steps
num_segments = 5           # Number of segments for potential/actual energies
k_noise = 0.1              # Noise filtering constant
k_feedback = 0.1           # Samson’s Law feedback constant
alpha_mem = 0.2            # Memory growth rate
C_const = 0.35             # Harmonic constant
fold_factor = 0.5          # Folding factor for quantum folding
zeta = 0.1                 # Residual term for unfolding

# Time array for reference (e.g., 0..19)
t_vals = np.arange(time_steps)

# -------------------------------
# 1. Generate Time-Varying Potential & Actual Energies
# -------------------------------
# Each time step has random potential (P) & actual (A) energies
# We'll store them in arrays: P[t, i], A[t, i]
P_data = np.random.uniform(1, 10, size=(time_steps, num_segments))
A_data = np.random.uniform(1, 10, size=(time_steps, num_segments))

# -------------------------------
# 2. Compute Mark 1 Resonance H(t) for each time step
#    H(t) = sum(P_i) / sum(A_i)
# -------------------------------
H_values = []
for t in range(time_steps):
    P_sum = np.sum(P_data[t])
    A_sum = np.sum(A_data[t])
    H_t = P_sum / A_sum
    H_values.append(H_t)
H_values = np.array(H_values)

# -------------------------------
# 3. Dynamic Noise Filtering (DNF)
#    We'll simulate random noise differences ΔN for each time step, then filter them.
# -------------------------------
Delta_N = np.random.uniform(-5, 5, size=time_steps)  # random noise differences
N_filtered = []
for t in range(time_steps):
    # Single-step DNF formula:
    # N(t) = ΔN / (1 + k * |ΔN|)
    # But let's accumulate over time or just store individually
    dn = Delta_N[t] / (1 + k_noise * abs(Delta_N[t]))
    N_filtered.append(dn)
N_filtered = np.array(N_filtered)

# -------------------------------
# 4. Samson’s Law (basic version)
#    ΔE = k * ΔF,  S = ΔE / T
#    We'll simulate ΔF for each time step, then compute S(t).
# -------------------------------
Delta_F = np.random.uniform(-10, 10, size=time_steps)
T_value = 2.0  # fixed time interval
Delta_E = k_feedback * Delta_F
S_values = Delta_E / T_value  # stabilization rate

# -------------------------------
# 5. Harmonic Memory Growth (HMG)
#    M(t) = M0 * exp(alpha * (H - C) * t)
#    We'll treat t as discrete steps, or do a cumulative approach.
# -------------------------------
M0 = 1.0
M_array = []
for t in range(time_steps):
    # Use the average H up to time t to represent the system's current resonance
    # Or you could just use H_values[t].
    H_avg = np.mean(H_values[: t + 1])
    M_t = M0 * np.exp(alpha_mem * (H_avg - C_const) * t)
    M_array.append(M_t)
M_array = np.array(M_array)

# -------------------------------
# 6. Quantum Folding/Unfolding at final time step
#    We'll demonstrate with the final P, A at t = time_steps - 1
# -------------------------------
final_t = time_steps - 1
P_final = P_data[final_t]
A_final = A_data[final_t]
H_final = H_values[final_t]

# Compute F(Q)_i for each segment:
FQ_components = (P_final / A_final) * np.exp(H_final * fold_factor * 1.0)
FQ_total = np.sum(FQ_components)

# Unfolding with random phase angles
theta = np.random.uniform(0, 2 * np.pi, size=num_segments)
UQ_value = np.sum(FQ_components * np.cos(theta)) + zeta

# -------------------------------
# Plot Results
# -------------------------------
fig, axs = plt.subplots(3, 2, figsize=(12, 10))
axs = axs.flatten()

# 1) H(t) over time
axs[0].plot(t_vals, H_values, marker='o')
axs[0].set_title("Mark 1 Resonance H(t)")
axs[0].set_xlabel("Time step")
axs[0].set_ylabel("H(t)")
axs[0].grid(True)

# 2) Dynamic Noise Filtering
axs[1].plot(t_vals, Delta_N, 'r--', label="Raw ΔN")
axs[1].plot(t_vals, N_filtered, 'b-', label="Filtered Noise")
axs[1].set_title("Dynamic Noise Filtering (DNF)")
axs[1].set_xlabel("Time step")
axs[1].set_ylabel("Noise Value")
axs[1].legend()
axs[1].grid(True)

# 3) Samson's Law stabilization rates
axs[2].plot(t_vals, Delta_F, 'r--', label="ΔF")
axs[2].plot(t_vals, S_values, 'b-', label="S = ΔE / T")
axs[2].set_title("Samson’s Law Over Time")
axs[2].set_xlabel("Time step")
axs[2].set_ylabel("Value")
axs[2].legend()
axs[2].grid(True)

# 4) Harmonic Memory Growth
axs[3].plot(t_vals, M_array, 'g-')
axs[3].set_title("Harmonic Memory Growth (HMG)")
axs[3].set_xlabel("Time step")
axs[3].set_ylabel("M(t)")
axs[3].grid(True)

# 5) Final Step: Quantum Folding Components
axs[4].bar(np.arange(num_segments), FQ_components, alpha=0.7, color='orange')
axs[4].set_title(f"Quantum Folding at t={final_t} (FQ_total={FQ_total:.2f})")
axs[4].set_xlabel("Segment Index")
axs[4].set_ylabel("F(Q)_i")
axs[4].grid(True)

# 6) Final Step: Unfolding with phase
cos_adjusted = FQ_components * np.cos(theta)
axs[5].bar(np.arange(num_segments), FQ_components, alpha=0.2, color='orange', label="F(Q)_i")
axs[5].scatter(np.arange(num_segments), cos_adjusted, color='red', label="F(Q)_i * cos(θ_i)")
axs
