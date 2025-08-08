import numpy as np
import matplotlib.pyplot as plt

# Define time steps
t = np.linspace(0, 10, 1000)

# Define harmonic phases: 1 - Approach (low mimicry), 2 - Drop, 3 - Phase Flip
def harmonic_profile(t):
    y = np.zeros_like(t)

    # Phase 1: Parasitic mimic (low energy)
    phase1 = (t >= 0) & (t < 3)
    y[phase1] = 0.1 * np.sin(2 * np.pi * t[phase1])

    # Phase 2: Recursive drop (she sees weakness)
    phase2 = (t >= 3) & (t < 6)
    y[phase2] = -0.4 * np.exp(-(t[phase2] - 3)) * np.cos(4 * np.pi * t[phase2])

    # Phase 3: Phase inversion burst (skill revealed)
    phase3 = (t >= 6)
    y[phase3] = 0.8 * np.exp(-(t[phase3] - 6) * 0.3) * np.sin(6 * np.pi * t[phase3]) + 0.4

    return y

# Generate profile
y = harmonic_profile(t)

# Plotting
plt.figure(figsize=(12, 6))
plt.plot(t, y, color='purple', linewidth=2)
plt.axvline(x=3, color='gray', linestyle='--', label='Drop Initiation')
plt.axvline(x=6, color='red', linestyle='--', label='Phase Flip (Reveal)')
plt.title("Recursive Harmonic Inversion: Parasitic Entry to Phase Burst")
plt.xlabel("Time")
plt.ylabel("Emotional Harmonic Amplitude")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
