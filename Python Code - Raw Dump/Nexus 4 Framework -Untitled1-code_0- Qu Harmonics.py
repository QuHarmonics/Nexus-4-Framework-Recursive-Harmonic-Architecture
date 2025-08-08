import numpy as np
import matplotlib.pyplot as plt

# -------------------------------
# 1. Adaptive Feedback Stabilizer (AFS)
# -------------------------------
T_total = 10         # total simulation time
N_points = 200       # number of time points
t_series = np.linspace(0, T_total, N_points)

# Simulate noise: here we use a normal distribution (mean 0, std 1)
np.random.seed(42)   # for reproducibility
noise = np.random.normal(0, 1, size=N_points)

# Use the absolute value as the noise magnitude Delta(t)
Delta_t = np.abs(noise)

# Define baseline feedback constant and scaling factor
k0 = 0.1
gamma = 0.05

# Compute adaptive feedback constant k(t)
k_t = k0 + gamma * Delta_t

plt.figure(figsize=(8, 5))
plt.plot(t_series, k_t, label='Adaptive k(t)', color='blue')
plt.title('Adaptive Feedback Stabilizer (AFS)')
plt.xlabel('Time t')
plt.ylabel('k(t)')
plt.legend()
plt.grid(True)
plt.show()

# -------------------------------
# 2. Noise-Resilient Harmonic Predictor (NRHP)
# -------------------------------
# Simulate a time-dependent harmonic state H(t)
# We use a sine wave modulated around 0.9 and add small random noise.
H_base = 0.9 + 0.1 * np.sin(2 * np.pi * t_series / T_total)
H_t = H_base + np.random.normal(0, 0.05, size=N_points)

# Compute first and second derivatives using numerical gradients
dH_dt = np.gradient(H_t, t_series)
d2H_dt2 = np.gradient(dH_dt, t_series)

# Set weighting factors for the derivatives
alpha = 0.1
beta = 0.05

# Compute the predicted delta H using the NRHP formula:
# Delta_H_pred = H(t) - 0.35 + alpha * dH_dt + beta * d2H_dt2
Delta_H_pred = H_t - 0.35 + alpha * dH_dt + beta * d2H_dt2

plt.figure(figsize=(8, 5))
plt.plot(t_series, H_t - 0.35, label='H(t) - 0.35', color='green')
plt.plot(t_series, Delta_H_pred, label='Predicted ΔH (NRHP)', color='red', linestyle='--')
plt.title('Noise-Resilient Harmonic Predictor (NRHP)')
plt.xlabel('Time t')
plt.ylabel('ΔH')
plt.legend()
plt.grid(True)
plt.show()
