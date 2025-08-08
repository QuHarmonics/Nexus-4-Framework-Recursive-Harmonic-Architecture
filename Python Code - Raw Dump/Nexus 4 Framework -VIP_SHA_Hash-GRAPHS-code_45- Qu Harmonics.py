# Re-running the adjusted waveforms and alignment history for SHA reconstruction analysis
import numpy as np
import matplotlib.pyplot as plt

# Constants (Normalized) from prior steps
constants = [0.27264203, 0.46389402, 0.74472339, 0.9576116, 0.23494206, 0.36852961, 0.59924109, 0.7011437]

# Example Final Waveform (from uploaded images for context)
final_waveform = np.array([
    [-0.04011812, 0.02865528, 0.14288313, 0.2227326, -0.09226802, -0.02109399, 0.087651, 0.14646431],
    [0.02865528, 0.13158962, 0.2940543, 0.41154043, -0.02462803, 0.06705249, 0.21326083, 0.28720288],
    [0.14288313, 0.2940543, 0.52299908, 0.69299634, 0.0927981, 0.21055904, 0.406464, 0.49874197],
    [0.2227326, 0.41154043, 0.69299634, 0.90416811, 0.17259968, 0.31217014, 0.54846415, 0.65648111],
    [-0.09226802, -0.02462803, 0.0927981, 0.17259968, -0.15813491, -0.07946322, 0.03755358, 0.10393196],
    [-0.02109399, 0.06705249, 0.21055904, 0.31217014, -0.07946322, 0.00668178, 0.14040057, 0.21110628],
    [0.087651, 0.21326083, 0.406464, 0.54846415, 0.03755358, 0.14040057, 0.30898473, 0.3906885],
    [0.14646431, 0.28720288, 0.49874197, 0.65648111, 0.10393196, 0.21110628, 0.3906885, 0.47245402]
])

# Iterative alignment adjustment
iterations = 50
std_dev_history = []
aligned_waveform = final_waveform.copy()

for i in range(iterations):
    # Apply harmonic corrections iteratively
    aligned_waveform += np.sin(aligned_waveform * constants[i % len(constants)])
    std_dev = np.std(aligned_waveform - final_waveform)
    std_dev_history.append(std_dev)

# Plot results
fig, axs = plt.subplots(2, 1, figsize=(8, 12))
# Adjusted Waveform Visualization
ax = axs[0]
ax.imshow(aligned_waveform, cmap='viridis', interpolation='nearest')
ax.set_title("Adjusted Waveform After Iterations")
ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")

# Alignment History
ax = axs[1]
ax.plot(range(1, iterations + 1), std_dev_history, marker='o', label='Std Dev to Harmonic Target')
ax.set_title("Alignment History")
ax.set_xlabel("Iteration")
ax.set_ylabel("Std Dev")
ax.legend()

plt.tight_layout()
plt.show()
