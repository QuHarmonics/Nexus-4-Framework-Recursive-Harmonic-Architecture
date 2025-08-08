import hashlib
import numpy as np
import matplotlib.pyplot as plt

# SHA-256 Construction Analysis
def sha256_one_round(message):
    """Simulates the first round of SHA-256 for reverse testing."""
    # Simplified: Operates on first block of message (512 bits)
    # Padding message to 512 bits
    padded_message = message.ljust(64, b'\x00')[:64]  # Truncate or pad
    hash_obj = hashlib.sha256()
    hash_obj.update(padded_message)
    return hash_obj.hexdigest()

# Reverse SHA-256 Single Round
def reverse_sha256_single_round(hash_value, guess_input):
    """Attempts to reverse-engineer a single round of SHA-256"""
    guess_hash = sha256_one_round(guess_input)
    if guess_hash == hash_value:
        return f"Recovered Input: {guess_input}"
    else:
        return f"Failed: {guess_input} does not match {hash_value}"

# Nyquist Expansion of Hash
def nyquist_expansion(hash_value):
    """Attempts to expand the hash as a frequency signal."""
    # Convert hash to numeric data
    hash_bytes = bytes.fromhex(hash_value)
    signal = np.array([int(byte) / 255.0 for byte in hash_bytes])  # Normalize

    # Oversample using interpolation
    oversampled_signal = np.interp(
        np.linspace(0, len(signal) - 1, len(signal) * 2),
        np.arange(len(signal)),
        signal,
    )

    # FFT Analysis
    fft = np.fft.fft(oversampled_signal)
    freqs = np.fft.fftfreq(len(oversampled_signal))

    # Plotting
    plt.figure(figsize=(12, 6))
    plt.subplot(1, 3, 1)
    plt.plot(signal, label="Original Signal")
    plt.legend()
    plt.subplot(1, 3, 2)
    plt.plot(oversampled_signal, label="Oversampled Signal")
    plt.legend()
    plt.subplot(1, 3, 3)
    plt.stem(freqs[:len(freqs)//2], np.abs(fft[:len(freqs)//2]), linefmt="r-", markerfmt="o", label="FFT Magnitudes")
    plt.legend()
    plt.show()

    return oversampled_signal

# Main Runner
def main():
    # Hash input
    data = "Hello, Nyquist and SHA-256!".encode()
    hash_value = hashlib.sha256(data).hexdigest()

    print(f"Hash Value: {hash_value}")

    # Reverse SHA-256 Example
    reverse_result = reverse_sha256_single_round(hash_value, data)
    print(f"Reverse SHA-256 Attempt: {reverse_result}")

    # Nyquist Test
    print("Performing Nyquist Expansion...")
    expanded_signal = nyquist_expansion(hash_value)

if __name__ == "__main__":
    main()
