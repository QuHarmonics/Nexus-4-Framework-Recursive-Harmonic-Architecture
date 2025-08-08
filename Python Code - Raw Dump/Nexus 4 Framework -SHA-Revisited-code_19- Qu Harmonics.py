import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# --- Parameters for the recursive 3-phase engine ---
time_steps = 1000                  # Simulation duration
carrier_frequency = 2 * np.pi      # Base frequency (c-like)
phase_shift = 2 * np.pi / 3        # 3-phase system (120 degrees apart)

# Introduce the 2% phase imbalance (recursive torque)
offset_ratio = 0.02
rotating_phase_offset = offset_ratio * 2 * np.pi

# Time domain
T = np.linspace(0, 10, time_steps)

# Recursive rotating offset — cycles over time
offset_rotation = np.sin(2 * np.pi * T / max(T)) * rotating_phase_offset

# Carrier wave with recursive rotation
wave_a = np.sin(carrier_frequency * T)
wave_b = np.sin(carrier_frequency * T + phase_shift + offset_rotation)
wave_c = np.sin(carrier_frequency * T + 2 * phase_shift - offset_rotation)

# Composite recursive harmonic driver
composite_wave = wave_a + wave_b + wave_c

# --- Plotting ---
fig, ax = plt.subplots(figsize=(12, 6))
ax.set_title("Recursive 3-Phase Harmonic with Rotating Offset")
ax.set_xlabel("Time")
ax.set_ylabel("Amplitude")

line_a, = ax.plot([], [], label='Wave A (0°)', color='cyan')
line_b, = ax.plot([], [], label='Wave B (+120° + rotating)', color='magenta')
line_c, = ax.plot([], [], label='Wave C (+240° - rotating)', color='yellow')
line_sum, = ax.plot([], [], label='Composite Harmonic', color='white', linewidth=2)

ax.set_xlim(0, 10)
ax.set_ylim(-4, 4)
ax.legend(loc="upper right")

# --- Animation update ---
def update(frame):
    idx = frame
    line_a.set_data(T[:idx], wave_a[:idx])
    line_b.set_data(T[:idx], wave_b[:idx])
    line_c.set_data(T[:idx], wave_c[:idx])
    line_sum.set_data(T[:idx], composite_wave[:idx])
    return line_a, line_b, line_c, line_sum

ani = FuncAnimation(fig, update, frames=len(T), interval=30, blit=True)
plt.show()
