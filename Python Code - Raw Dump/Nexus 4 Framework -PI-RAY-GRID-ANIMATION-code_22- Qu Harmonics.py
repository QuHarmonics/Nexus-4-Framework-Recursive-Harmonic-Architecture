import numpy as np
import matplotlib.pyplot as plt

# System parameters
m = 1.0  # mass (kg)
k = 100.0  # spring constant (N/m)
c_initial = 2.0  # initial damping coefficient (kg/s)
dt = 0.01  # time step (s)
t_max = 5.0  # simulation time (s)

# Nexus constants
zeta_target = 0.35
kappa = 10.0  # smoothing strength

# Calculate target damping coefficient (Samson V2)
c_target = 2 * np.sqrt(k * m) * zeta_target  # 7.0 kg/s

# Initialize arrays
t = np.arange(0, t_max, dt)
x = np.zeros(len(t))  # position
v = np.zeros(len(t))  # velocity
c = np.zeros(len(t))  # damping coefficient over time
x[0] = 1.0  # initial displacement
c[0] = c_initial

# Simulation loop
for n in range(len(t) - 1):
    # Current damping ratio
    zeta_n = c[n] / (2 * np.sqrt(k * m))
    
    # Mary’s Spirit: Smoothly adjust c towards c_target
    c[n + 1] = c[n] * (1 + np.exp(-kappa * (zeta_n - zeta_target)))
    if abs(zeta_n - zeta_target) < 0.01:  # Clamp when close
        c[n + 1] = c_target
    
    # Equations of motion
    a = (-k * x[n] - c[n] * v[n]) / m
    v[n + 1] = v[n] + a * dt
    x[n + 1] = x[n] + v[n] * dt

# Plot results
plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.plot(t, x, label="Displacement")
plt.title("Nexus-Stabilized Oscillator")
plt.xlabel("Time (s)")
plt.ylabel("Displacement (m)")
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(t, c, label="Damping Coefficient", color="orange")
plt.axhline(y=c_target, color="r", linestyle="--", label="Target c (7.0)")
plt.xlabel("Time (s)")
plt.ylabel("Damping Coefficient (kg/s)")
plt.legend()

plt.tight_layout()
plt.show()