import numpy as np
import matplotlib.pyplot as plt

# Step 1: Encode the Binary Data into H (Storage)
def store_in_H(binary_data, expansion_factor=1.5):
    harmonics = np.cumsum(binary_data.astype(np.float64) * expansion_factor)  # Use higher precision
    return harmonics

# Step 2: Reverse the Process (Retrieve Original Data)
def retrieve_from_H(harmonics, expansion_factor=1.5):
    first_value = harmonics[0] / expansion_factor
    reversed_data = np.diff(harmonics) / expansion_factor
    reversed_data = np.insert(reversed_data, 0, first_value)
    return np.round(reversed_data).astype(np.uint8)

# Generate Found Wave as Quantum Input
def found_wave_to_binary(found_wave):
    # Normalize found wave and convert to binary (0/1 for simplicity)
    threshold = np.mean(found_wave)
    binary_data = np.array([1 if point > threshold else 0 for point in found_wave], dtype=np.uint8)
    return binary_data

# Input: Found Quantum Wave
found_wave = np.sin(np.linspace(0, 2 * np.pi, 512))  # Example wave (replace with actual found wave)
binary_found_wave = found_wave_to_binary(found_wave)

# Step 3: Push Found Wave into H Storage
harmonics = store_in_H(binary_found_wave)

# Step 4: Retrieve Macro Data from H
retrieved_macro_binary = retrieve_from_H(harmonics)

# Visualize Found Wave and Retrieved Macro Binary
plt.figure(figsize=(12, 6))
plt.subplot(211)
plt.plot(found_wave, label="Found Quantum Wave", color='blue')
plt.title("Found Quantum Wave")
plt.legend()

plt.subplot(212)
plt.plot(retrieved_macro_binary, label="Retrieved Macro Binary", color='green')
plt.title("Retrieved Macro Binary")
plt.legend()
plt.show()

# Output Results
print("Found Wave Binary (First 100 bits):", binary_found_wave[:10000])
print("Retrieved Macro Binary (First 100 bits):", retrieved_macro_binary[:10000])
