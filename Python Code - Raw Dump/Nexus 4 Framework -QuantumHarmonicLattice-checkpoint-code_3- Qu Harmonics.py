from tqdm import tqdm
import wave
import numpy as np

# Constants
HARMONIC_CONSTANT = 0.35
TOLERANCE = 0.01
MAX_ITERATIONS = 100

# 1. Read WAV File
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

# 3. Convert Audio to Binary with Progress Meter
def audio_to_binary(audio_data):
    """Converts signed audio data to binary string with progress tracking."""
    max_bits = 16 if audio_data.dtype == np.int16 else 8
    binary_data = ''.join(
        format(sample & (2**max_bits - 1), f'0{max_bits}b') 
        for sample in tqdm(audio_data, desc="Converting audio to binary")
    )
    return binary_data

# 4. Convert Binary to Matrix with Progress Meter
def binary_to_matrix(binary_data):
    """Converts binary data to a 2D matrix with progress tracking."""
    binary_length = len(binary_data)
    matrix_size = int(np.ceil(np.sqrt(binary_length)))
    data_matrix = np.zeros((matrix_size, matrix_size), dtype=np.float64)

    for i, bit in tqdm(enumerate(binary_data), total=binary_length, desc="Filling matrix"):
        data_matrix[i // matrix_size, i % matrix_size] = int(bit)

    return data_matrix

# 5. Harmonize Data with Progress Meter
def harmonize_data(data, harmonic_constant, compress=True):
    """Compress or expand data by aligning it with the harmonic constant, with progress tracking."""
    data = np.copy(data)
    gain = 1.0 if compress else -1.0

    for iteration in tqdm(range(MAX_ITERATIONS), desc="Harmonic adjustment"):
        delta = np.mean(data) - harmonic_constant
        adjustment = delta * gain
        data -= adjustment
        data = np.clip(data, 0, 1)  # Ensure binary range
        if abs(delta) < TOLERANCE:
            break

    return data

# 6. Binary to Audio
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

# Main Workflow
if __name__ == "__main__":
    # Input and Output Files
    input_wav_file = "d:\\test.wav"
    compressed_file = "d:\\compressed_audio.npz"
    output_wav_file = "d:\\restored.wav"

    # Step 1: Read the Input WAV File
    audio_data, n_channels, sample_width, framerate = read_wav_file(input_wav_file)
    print(f"Read WAV File: {input_wav_file}")
    print(f"Channels: {n_channels}, Sample Width: {sample_width}, Frame Rate: {framerate}")

    # Step 2: Convert Audio to Binary
    binary_data = audio_to_binary(audio_data)

    # Step 3: Convert Binary to Matrix
    data_matrix = binary_to_matrix(binary_data)

    # Step 4: Apply Harmonic Compression
    compressed_data = harmonize_data(data_matrix, HARMONIC_CONSTANT, compress=True)
    print("Compressed Data Non-Zero Elements:", np.count_nonzero(compressed_data))

    # Step 5: Save Compressed Data
    np.savez_compressed(compressed_file, compressed_matrix=compressed_data)
    print(f"Compressed data saved to {compressed_file}")

    # Step 6: Load Compressed Data and Expand
    loaded_data = np.load(compressed_file)['compressed_matrix']
    expanded_data = harmonize_data(loaded_data, HARMONIC_CONSTANT, compress=False)

    # Step 7: Flatten Expanded Data and Convert Back to Audio
    flattened_binary = ''.join(str(int(round(bit))) for bit in tqdm(expanded_data.flatten(), desc="Flattening matrix"))
    restored_audio = binary_to_audio(flattened_binary, audio_data.dtype)

    # Step 8: Save Restored Audio to a New WAV File
    save_wav_file(output_wav_file, restored_audio, n_channels, sample_width, framerate)
    print(f"Restored WAV File Saved: {output_wav_file}")
