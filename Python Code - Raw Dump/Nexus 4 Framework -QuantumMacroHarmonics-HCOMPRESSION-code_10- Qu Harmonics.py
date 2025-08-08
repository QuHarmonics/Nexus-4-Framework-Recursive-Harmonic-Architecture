import numpy as np
from scipy.io.wavfile import read, write
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

# File paths
input_file_path = r'd:\test.mp3'  # Input wave file
compressed_file_path = r'd:\compressed_audio.hc'  # Compressed harmonic file
decoded_file_path = r'd:\decoded_audio.wav'  # Decoded wave file

# Load the wave file as raw data
sample_rate, audio_data = read(input_file_path)

# Ensure no modifications are made to the incoming data
# All channels and raw data will be processed as-is

# Step 1: Encode the Audio Data into H (Storage)
def store_in_H(audio_data, expansion_factor=1.5):
    harmonics = np.cumsum(audio_data.astype(np.float64) * expansion_factor)  # High precision cumulative sum
    return harmonics

# Step 2: Reverse the Process (Retrieve Original Audio Data)
def retrieve_from_H(harmonics, expansion_factor=1.5):
    first_value = harmonics[0] / expansion_factor  # Correct the first value
    reversed_data = np.diff(harmonics) / expansion_factor
    reversed_data = np.insert(reversed_data, 0, first_value)  # Insert corrected first value
    return np.round(reversed_data).astype(audio_data.dtype)  # Ensure original data type

# Encode the audio data into H
harmonics = store_in_H(audio_data)

# Save the compressed harmonic data as .hc
harmonics.astype(np.float64).tofile(compressed_file_path)

# Decode the audio data from H
retrieved_audio = retrieve_from_H(harmonics)

# Save the decoded audio to a new wave file
write(decoded_file_path, sample_rate, retrieved_audio)

# Compare File Sizes
original_size = os.path.getsize(input_file_path)
compressed_size = os.path.getsize(compressed_file_path)
decoded_size = os.path.getsize(decoded_file_path)

print(f"Original Audio File Size: {original_size / 1024:.2f} KB")
print(f"Compressed Harmonics File Size: {compressed_size / 1024:.2f} KB")
print(f"Decoded Audio File Size: {decoded_size / 1024:.2f} KB")

# Visualize H(n) in 3D
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Generate 3D coordinates
x = np.arange(len(harmonics))
y = harmonics
z = np.sin(x / 100.0)  # Add depth with a sinusoidal wave

# Plot the 3D visualization
ax.plot(x, y, z, label="H(n) in 3D", color='b', lw=2)
ax.scatter(x, y, z, color='r', s=5, label="Nodes")

# Labels and title
ax.set_title("3D Visualization of H(n)", fontsize=16)
ax.set_xlabel("Iteration (n)", fontsize=12)
ax.set_ylabel("H(n)", fontsize=12)
ax.set_zlabel("Z-axis Wave", fontsize=12)
ax.legend()

plt.show()

# Visualize the original and retrieved audio waveforms
plt.figure(figsize=(10, 6))
plt.plot(audio_data[:1000], label='Original Audio', alpha=0.7)
plt.plot(retrieved_audio[:1000], label='Decoded Audio', linestyle='dashed', alpha=0.7)
plt.legend()
plt.title("Original vs Decoded Audio Waveform")
plt.show()

print("Original and Decoded audio data match:", np.array_equal(audio_data, retrieved_audio))
