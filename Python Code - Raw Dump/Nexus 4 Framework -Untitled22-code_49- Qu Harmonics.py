import numpy as np
import matplotlib.pyplot as plt

# Constants
num_rounds = 50
input_data = np.random.rand(8, 8)  # Simulating an 8x8 input lattice
padding_data = np.zeros_like(input_data)  # Padding waste lattice (zeros for now)
harmonic_constant = 0.35  # Resonance constant

# Initialize waveforms
waveforms = [input_data.copy()]
padding_waveforms = [padding_data.copy()]

# Function to simulate SHA round (waveform transformation)
def sha_round_transform(data, padding, round_num):
    # Transform input lattice using harmonic and padding interaction
    new_data = np.abs(np.sin(data + harmonic_constant * round_num)) + padding
    new_padding = np.abs(np.cos(padding - harmonic_constant * round_num))
    return new_data, new_padding

# Perform SHA-like transformation over multiple rounds
for round_num in range(1, num_rounds + 1):
    new_data, new_padding = sha_round_transform(waveforms[-1], padding_waveforms[-1], round_num)
    waveforms.append(new_data)
    padding_waveforms.append(new_padding)

# Visualize the final waveform and padding
fig = plt.figure(figsize=(12, 6))

# Final waveform
ax1 = fig.add_subplot(121, projection='3d')
X, Y = np.meshgrid(range(input_data.shape[0]), range(input_data.shape[1]))
ax1.plot_surface(X, Y, waveforms[-1], cmap="viridis")
ax1.set_title("Final Waveform (Aligned with Output)")

# Padding waveform
ax2 = fig.add_subplot(122, projection='3d')
ax2.plot_surface(X, Y, padding_waveforms[-1], cmap="plasma")
ax2.set_title("Final Padding (Waste Alignment)")

plt.show()

# Analyze how the padding shifts results in iterations
padding_std = [np.std(pad) for pad in padding_waveforms]
waveform_std = [np.std(wave) for wave in waveforms]

# Plot entropy trends
plt.figure(figsize=(10, 6))
plt.plot(range(num_rounds + 1), waveform_std, label="Waveform Entropy", marker='o')
plt.plot(range(num_rounds + 1), padding_std, label="Padding Entropy", marker='x')
plt.axhline(harmonic_constant, color='red', linestyle='--', label="Harmonic Constant (Target)")
plt.title("Entropy Trends Across Iterations")
plt.xlabel("Iteration")
plt.ylabel("Entropy (Std Dev)")
plt.legend()
plt.grid(True)
plt.show()

# Return final state for further analysis
final_state = {
    "Waveform": waveforms[-1],
    "Padding": padding_waveforms[-1],
    "Waveform Entropy": waveform_std,
    "Padding Entropy": padding_std
}

final_state
