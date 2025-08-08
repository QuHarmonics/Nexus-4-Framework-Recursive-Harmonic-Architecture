import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft
from mpl_toolkits.mplot3d import Axes3D

# Constants defined by the SHA-256 algorithm
K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
]

# Differences and ratios
differences = [K[i] - K[i - 1] for i in range(1, len(K))]
ratios = [K[i] / K[i - 1] for i in range(1, len(K))]

# Waveform
waveform_y = np.array(K) / max(K)
frequencies = np.fft.fftfreq(len(waveform_y), d=1)
fft_values = fft(waveform_y)

# Dominant harmonics
dominant_frequencies = np.argsort(-np.abs(fft_values))[:5]
dominant_magnitudes = np.abs(fft_values[dominant_frequencies])
dominant_frequencies_values = frequencies[dominant_frequencies]

# Harmonic behavior interpretation
def harmonic_behavior_analysis(frequencies, magnitudes):
    behaviors = []
    for freq, mag in zip(frequencies, magnitudes):
        if freq == 0.0:
            behaviors.append("Baseline stability (central alignment)")
        elif abs(freq) < 0.1:
            behaviors.append("Low-frequency oscillations (structural reinforcement)")
        else:
            behaviors.append("Higher frequency oscillations (fine-tuning or diffusion)")
    return behaviors

harmonic_behaviors = harmonic_behavior_analysis(
    dominant_frequencies_values,
    dominant_magnitudes
)

# ASM simulation
def simulate_advanced_kinetics(asm_instructions, initial_value=0):
    value = initial_value
    kinetic_progress = []
    for instruction in asm_instructions:
        if instruction.startswith("add"):
            delta = int(instruction.split(" ")[1], 16)
            value += delta
        elif instruction.startswith("sub"):
            delta = int(instruction.split(" ")[1], 16)
            value -= delta
        elif instruction.startswith("nop"):
            value = value
        elif instruction.startswith("xor"):
            delta = int(instruction.split(" ")[1], 16)
            value ^= delta
        elif instruction.startswith("shl"):
            shift = int(instruction.split(" ")[1], 16)
            value <<= shift
        kinetic_progress.append(value)
    return kinetic_progress

asm_instructions = [f"add {delta:#x}" if delta > 0 else f"sub {abs(delta):#x}" for delta in differences]
kinetic_progress = simulate_advanced_kinetics(asm_instructions[:64])
kinetic_normalized = np.array(kinetic_progress) / max(abs(np.array(kinetic_progress)))

# Surface matrix
interaction_waveform = np.outer(waveform_y, kinetic_normalized)
Z_interaction = interaction_waveform

# Attractor highlights near 0.35
collapse_attractors = np.isclose(Z_interaction, 0.35, atol=0.05)

# Meshgrid
X = np.linspace(0, Z_interaction.shape[1] - 1, Z_interaction.shape[1])
Y = np.linspace(0, Z_interaction.shape[0] - 1, Z_interaction.shape[0])
X, Y = np.meshgrid(X, Y)

# Plot
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, Z_interaction, cmap='plasma', edgecolor='none')
ax.scatter(X[collapse_attractors], Y[collapse_attractors], Z_interaction[collapse_attractors], color='cyan', s=5)
ax.set_title("Refined SHA Harmonic Surface with Attractor Overlay")
ax.set_xlabel("X Axis (Index)")
ax.set_ylabel("Y Axis (Index)")
ax.set_zlabel("Amplitude")
plt.show()

# Surface area
def calculate_surface_area(Z):
    return np.sum(np.sqrt(1 + np.gradient(Z, axis=0)**2 + np.gradient(Z, axis=1)**2))

surface_area = calculate_surface_area(Z_interaction)

# Output
results = {
    "Surface Area Approximation": surface_area,
    "Dominant Frequencies": dominant_frequencies_values.tolist(),
    "Dominant Magnitudes": dominant_magnitudes.tolist(),
    "Harmonic Behaviors": harmonic_behaviors,
    "Collapse Attractors (≈0.35)": int(np.sum(collapse_attractors))
}

print(results)
