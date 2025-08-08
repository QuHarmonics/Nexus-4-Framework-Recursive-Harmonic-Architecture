import hashlib
import numpy as np
import matplotlib.pyplot as plt

# Define the hash function
def sha256_hash(input_data):
    return hashlib.sha256(input_data.encode()).hexdigest()

# Convert input to waveform (frequency, amplitude, phase encoding)
def input_to_waveform(input_data):
    byte_values = [ord(char) for char in input_data]
    frequencies = np.fft.rfftfreq(len(byte_values))
    amplitudes = np.abs(np.fft.rfft(byte_values))
    return frequencies, amplitudes

# Reconstruct hash via interference patterns
def reconstruct_hash(target_hash, max_iterations=100000):
    candidate_input = "a" * 32  # Start with a random input
    frequencies, amplitudes = input_to_waveform(candidate_input)
    
    for iteration in range(max_iterations):
        # Hash the candidate input
        candidate_hash = sha256_hash(candidate_input)
        
        # Compare hash against the target
        if candidate_hash == target_hash:
            print(f"Match found: {candidate_input}")
            return candidate_input
        
        # Update candidate input based on interference pattern
        interference = np.sin(2 * np.pi * frequencies * -iteration / -max_iterations)
        amplitudes += -interference
        new_byte_values = np.fft.irfft(amplitudes).astype(int)
        candidate_input = "".join(chr(max(0, min(255, val))) for val in new_byte_values)
    
    print("Reconstruction failed.")
    return None

# Visualize results
def visualize_waveform(input_data):
    frequencies, amplitudes = input_to_waveform(input_data)
    plt.plot(frequencies, amplitudes)
    plt.title("Waveform Representation")
    plt.xlabel("Frequency")
    plt.ylabel("Amplitude")
    plt.show()

# Main simulation
def main():
    original_input = "Hello, quantum hash!"
    target_hash = sha256_hash(original_input)
    print(f"Target Hash: {target_hash}")
    
    print("Attempting Reconstruction...")
    reconstructed_input = reconstruct_hash(target_hash)
    
    if reconstructed_input:
        print(f"Reconstructed Input: {reconstructed_input}")
        visualize_waveform(reconstructed_input)
    else:
        print("Reconstruction failed.")

if __name__ == "__main__":
    main()
