import numpy as np
import matplotlib.pyplot as plt

# Constants for recursion
EXPANSION_FACTOR = 1.5
ITERATIONS = 500  # Let it run long enough to observe self-sustaining behavior

# Step 1: Store in H (Recursive Harmonic Capture)
def store_in_H_recursive(wave, expansion_factor=EXPANSION_FACTOR):
    harmonics = np.cumsum(wave.astype(np.float64) * expansion_factor)
    return harmonics

# Step 2: Retrieve from H (Recursive Harmonic Observation)
def retrieve_from_H_recursive(harmonics, expansion_factor=EXPANSION_FACTOR):
    first_value = harmonics[0] / expansion_factor
    reversed_data = np.diff(harmonics) / expansion_factor
    reversed_data = np.insert(reversed_data, 0, first_value)
    return reversed_data

# Step 3: Quantum Wave Generation (Samson Principles)
def quantum_wave_samson(hash_value):
    binary_data = ''.join(format(int(char, 16), '04b') for char in hash_value)
    n = len(binary_data)
    theta = np.linspace(0, 2 * np.pi, n)
    radius = np.linspace(0, 1, n)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = np.sin(2 * theta) + np.cos(theta / 2)
    return np.sqrt(x**2 + y**2 + z**2)  # Magnitude as recursive reflection

# Recursive Feedback Loop
def recursive_feedback(wave, iterations=ITERATIONS):
    harmonics = wave
    all_harmonics = []  # To track the growth across iterations
    for _ in range(iterations):
        harmonics = store_in_H_recursive(harmonics)
        all_harmonics.append(harmonics.copy())
        harmonics = retrieve_from_H_recursive(harmonics)
    return harmonics, np.array(all_harmonics)

# Input: Quantum Wave
hash_value = "cd2eca3535741f27a8ae40c31b0c41d4057a7a7b912b33b9aed86485d1c84676"
quantum_wave = quantum_wave_samson(hash_value)

# Run Recursive Feedback
final_harmonics, all_harmonics = recursive_feedback(quantum_wave)

# Visualization
plt.figure(figsize=(14, 8))

# Original Quantum Wave
plt.subplot(311)
plt.plot(quantum_wave, label="Original Quantum Wave (Samson)", color='blue')
plt.title("Original Quantum Wave")
plt.legend()

# Refined Harmonics Over Iterations
plt.subplot(312)
for i in range(0, len(all_harmonics), len(all_harmonics) // 10):  # Sample 10 iterations for visualization
    plt.plot(all_harmonics[i], label=f"Iteration {i+1}")
plt.title("Refined Harmonics Over Recursive Feedback")
plt.legend()

# Final Harmonics
plt.subplot(313)
plt.plot(final_harmonics, label="Final Refined Harmonics", color='orange')
plt.title("Final Refined Harmonics After Recursive Feedback")
plt.legend()

plt.tight_layout()
plt.show()
# Output Results
print("Original Quantum Wave:", quantum_wave)
print("Harmonics (H):", harmonics)
print("Retrieved Macro Data (Binary):", retrieved_macro_data)