import hashlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, ifft

# Generate a hash (SHA-256) from a given input
def generate_hash(input_data):
    hash_object = hashlib.sha256(input_data.encode())
    return hash_object.hexdigest()

# Convert the hash into a numerical signal
def hash_to_signal(hash_value):
    # Convert hex characters to integers and normalize
    signal = [int(hash_value[i:i+2], 16) / 255.0 for i in range(0, len(hash_value), 2)]
    return np.array(signal)

# Nyquist oversampling to double the resolution
def nyquist_oversample(signal, factor=2):
    original_length = len(signal)
    oversampled_signal = np.interp(
        np.linspace(0, original_length, original_length * factor, endpoint=False),
        np.arange(original_length),
        signal
    )
    return oversampled_signal

# Perform FFT to analyze harmonics
def analyze_harmonics(signal, sampling_rate):
    N = len(signal)
    fft_values = fft(signal)
    fft_frequencies = fftfreq(N, 1 / sampling_rate)

    # Return only the positive frequencies and normalize
    positive_freqs = fft_frequencies[:N // 2]
    magnitudes = np.abs(fft_values[:N // 2]) / N
    return positive_freqs, magnitudes

# Reconstruct the signal from FFT (optional)
def reconstruct_signal(fft_values, N):
    reconstructed_signal = ifft(fft_values)
    return np.real(reconstructed_signal)

# Main function to run the analysis
def main():
    # Input and hash generation
    input_data = "Hello, Nyquist and FFT!"  # Example input
    hash_value = generate_hash(input_data)
    print(f"Input Data: {input_data}")
    print(f"Hash Value: {hash_value}")

    # Convert hash to numerical signal
    signal = hash_to_signal(hash_value)
    print(f"Original Signal Length: {len(signal)}")

    # Nyquist oversampling
    oversampled_signal = nyquist_oversample(signal, factor=2)
    print(f"Oversampled Signal Length: {len(oversampled_signal)}")

    # Analyze harmonics using FFT
    sampling_rate = len(oversampled_signal)  # Assume one sample per unit
    freqs, magnitudes = analyze_harmonics(oversampled_signal, sampling_rate)

    # Plot the results
    plt.figure(figsize=(15, 5))

    # Plot original signal
    plt.subplot(1, 3, 1)
    plt.plot(signal, label="Original Signal")
    plt.title("Original Signal")
    plt.xlabel("Index")
    plt.ylabel("Amplitude")
    plt.legend()

    # Plot oversampled signal
    plt.subplot(1, 3, 2)
    plt.plot(oversampled_signal, label="Oversampled Signal")
    plt.title("Oversampled Signal")
    plt.xlabel("Index")
    plt.ylabel("Amplitude")
    plt.legend()

    # Plot FFT magnitudes
    plt.subplot(1, 3, 3)
    plt.stem(freqs, magnitudes, label="FFT Magnitudes")  # Fixed line
    plt.title("Harmonic Analysis (FFT)")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Magnitude")
    plt.legend()

    plt.tight_layout()
    plt.show()

# Run the program
if __name__ == "__main__":
    main()
