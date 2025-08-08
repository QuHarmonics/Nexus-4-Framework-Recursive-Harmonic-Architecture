import numpy as np
from scipy.signal import chirp
from scipy.fftpack import ifft
import matplotlib.pyplot as plt

# Constants
golden_ratio = 1.618  # Golden ratio for recursive growth
expansion_factor = 1.5  # Key parameter for unfolding chaos
sampling_rate = 1024  # Resolution for waveform reconstruction

# Generate chaos from harmonized hash
def generate_chaos(hash_data):
    """
    Reverse process: Expand a hash (harmonized data) into its chaotic original form.
    """
    # Decode the hash into a numerical sequence
    hash_sequence = np.array([int(hash_data[i:i+2], 16) for i in range(0, len(hash_data), 2)], dtype=np.float64)
    
    # Initialize chaotic expansion
    chaos = np.zeros(sampling_rate)
    time = np.linspace(0, 1, sampling_rate)
    
    # Recursive unfolding using harmonics
    for i, value in enumerate(hash_sequence):
        frequency = golden_ratio * value / max(hash_sequence)  # Frequency scales with the hash value
        amplitude = expansion_factor * value / 255.0  # Normalize amplitude to [0, expansion_factor]
        chaos += amplitude * chirp(time, f0=frequency, f1=frequency * 2, t1=1, method='linear')

    # Normalize the result to prevent overflow
    chaos = chaos / np.max(np.abs(chaos))
    return chaos

# Visualization of the chaotic reconstruction
def visualize_chaos(chaos_waveform):
    plt.figure(figsize=(10, 6))
    plt.plot(chaos_waveform, color='purple', label="Reconstructed Chaos")
    plt.title("Reconstructed Chaos from Harmonized Data", fontsize=16)
    plt.xlabel("Time Steps", fontsize=12)
    plt.ylabel("Amplitude", fontsize=12)
    plt.legend()
    plt.grid()
    plt.show()

# Input SHA hash (e.g., "0753e9bcdebcbe1bbbf19b6bfe6c06bd9606ce89f3c6d57d6d03a6c257f44b3f")
hash_input = "0753e9bcdebcbe1bbbf19b6bfe6c06bd9606ce89f3c6d57d6d03a6c257f44b3f"

# Generate and visualize chaos
chaos_waveform = generate_chaos(hash_input)
visualize_chaos(chaos_waveform)
