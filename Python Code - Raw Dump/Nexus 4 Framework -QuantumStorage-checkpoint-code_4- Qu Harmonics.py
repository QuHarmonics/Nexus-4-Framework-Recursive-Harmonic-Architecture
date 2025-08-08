import wave
import numpy as np
import os

# Constants
HARMONIC_CONSTANT = 0.35
TOLERANCE = 0.01
MAX_ITERATIONS = 100

# 1. Read the WAV File
def read_wav_file(filename):
    """Reads a .wav file and extracts audio data."""
    with wave.open(filename, 'rb') as wav_file:
        n_channels = wav_file.getnchannels()
        sample_width = wav_file.getsampwidth()
        framerate = wav_file.getframerate()
        n_frames = wav_file.getnframes()

        # Extract raw audio data
        raw_data = wav_file.readframes(n_frames)
        
        # Convert bytes to integer data
        dtype = np.int16 if sample_width == 2 else np.int8
        audio_data = np.frombuffer(raw_data, dtype=dtype)
        
        return audio_data, n_channels, sample_width, framerate

# 2. Save WAV File
def save_wav_file(filename, audio_data, n_channels, sample_width, framerate):
    """Saves audio data back to a .wav file."""
    with wave.open(filename, 'wb') as wav_file:
        wav_file.setnchannels(n_channels)
        wav_file.setsampwidth(sample_width)
        wav_file.setframerate(framerate)

        # Convert integer data to bytes and write
        wav_file.writeframes(audio_data.tobytes())

# 3. Convert Audio Data to Binary
def audio_to_binary(audio_data):
    """Converts signed audio data to binary string."""
    max_bits = 16 if audio_data.dtype == np.int16 else 8
    binary_data = ''.join(
        format(sample & (2**max_bits - 1), f'0{max_bits}b') for sample in audio_data
    )
    return binary_data

# 4. Convert Binary Back to Audio
def binary_to_audio(binary_data, dtype):
    """Converts binary string back to signed audio data."""
    max_bits = 16 if dtype == np.int16 else 8
    num_samples = len(binary_data) // max_bits
    audio_data = []

    for i in range(num_samples):
        value = int(binary_data[i * max_bits:(i + 1) * max_bits], 2)
        # Convert from two's complement if necessary
        if value >= 2**(max_bits - 1):
            value -= 2**max_bits
        audio_data.append(value)

    return np.array(audio_data, dtype=dtype)

# 5. Harmonic Compression/Expansion
def harmonize_data(data, harmonic_constant, compress=True):
    """
    Compress or expand data by aligning it with the harmonic constant.
    If compress=True, align data for compression.
    If compress=False, expand compressed data.
    """
    data = np.copy(data)
    gain = 1.0 if compress else -1.0

    for _ in range(MAX_ITERATIONS):
        delta = np.mean(data) - harmonic_constant
        adjustment = delta * gain
        data -= adjustment
        data = np.clip(data, 0, 1)  # Ensure data stays within valid range
        if abs(delta) < TOLERANCE:
            break

    return data

# 6. Compress and Save Data
def compress_audio(audio_data, filename):
    """Compresses audio data, saves it to an .npz file, and returns compression stats."""
    print("Original Audio Length:", len(audio_data))

    # Convert to binary
    binary_data = audio_to_binary(audio_data)
    print("Binary Length:", len(binary_data))

    # Convert binary data to matrix
    binary_length = len(binary_data)
    matrix_size = int(np.ceil(np.sqrt(binary_length)))
    data_matrix = np.zeros((matrix_size, matrix_size), dtype=np.float64)

    # Fill matrix with binary data
    for i, bit in enumerate(binary_data):
        data_matrix[i // matrix_size, i % matrix_size] = int(bit)

    # Harmonic compression
    compressed_data = harmonize_data(data_matrix, HARMONIC_CONSTANT, compress=True)
    print("Compressed Data Non-Zero Elements:", np.count_nonzero(compressed_data))

    # Save compressed data
    np.savez_compressed(filename, compressed_matrix=compressed_data)
    print(f"Compressed data saved to {filename}")

    return compressed_data

# 7. Expand and Restore Audio
def expand_audio(filename, original_dtype):
    """Expands compressed audio data and restores it."""
    # Load compressed data
    loaded_data = np.load(filename)['compressed_matrix']

    # Harmonic expansion
    expanded_data = harmonize_data(loaded_data, HARMONIC_CONSTANT, compress=False)
    
    # Flatten and convert back to binary
    flattened_binary = ''.join(str(int(round(bit))) for bit in expanded_data.flatten())
    
    # Convert binary back to audio data
    restored_audio = binary_to_audio(flattened_binary, original_dtype)
    print("Restored Audio Length:", len(restored_audio))

    return restored_audio

# Main Framework
if __name__ == "__main__":
    # Input and Output Files
    input_wav_file = "d:\\test.wav"
    compressed_file = "d:\\compressed_audio.npz"
    output_wav_file = "d:\\restored.wav"

    # Read the input WAV file
    audio_data, n_channels, sample_width, framerate = read_wav_file(input_wav_file)
    print("Read WAV File:", input_wav_file)

    # Compress and save the audio
    compressed_data = compress_audio(audio_data, compressed_file)

    # Expand and restore the audio
    restored_audio = expand_audio(compressed_file, audio_data.dtype)

    # Save the restored audio back to a new WAV file
    save_wav_file(output_wav_file, restored_audio, n_channels, sample_width, framerate)
    print("Restored WAV File Saved:", output_wav_file)
