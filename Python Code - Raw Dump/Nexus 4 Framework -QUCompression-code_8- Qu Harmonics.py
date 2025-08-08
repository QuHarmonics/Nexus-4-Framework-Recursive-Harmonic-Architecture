from tqdm import tqdm
import numpy as np
import os

# Constants
HARMONIC_CONSTANT = 0.35
TOLERANCE = 0.01
MAX_ITERATIONS = 100

# Utility Functions
def file_to_base16(filename):
    """
    Read any file as raw bytes and convert it to Base-16 representation with a progress meter.
    """
    file_size = os.path.getsize(filename)
    base16_data = []

    with open(filename, 'rb') as f, tqdm(total=file_size, desc="Encoding file to Base-16") as pbar:
        while chunk := f.read(1024):  # Read in chunks
            base16_data.extend(format(byte, '02x') for byte in chunk)
            pbar.update(len(chunk))
    
    return ''.join(base16_data)

def base16_to_file(base16_data, output_filename):
    """
    Decode Base-16 representation back to the original file format with a progress meter.
    """
    byte_count = len(base16_data) // 2
    with open(output_filename, 'wb') as f, tqdm(total=byte_count, desc="Restoring file from Base-16") as pbar:
        for i in range(0, len(base16_data), 2):
            byte = int(base16_data[i:i+2], 16)
            f.write(byte.to_bytes(1, 'big'))
            pbar.update(1)

def base16_to_matrix(base16_data):
    """
    Convert Base-16 data into a binary matrix for harmonic compression.
    """
    binary_data = ''.join(format(int(char, 16), '04b') for char in base16_data)
    matrix_size = int(np.ceil(np.sqrt(len(binary_data))))
    
    data_matrix = np.zeros((matrix_size, matrix_size), dtype=np.float64)
    for i, bit in tqdm(enumerate(binary_data), total=len(binary_data), desc="Filling binary matrix"):
        row, col = divmod(i, matrix_size)
        data_matrix[row, col] = int(bit)
    
    return data_matrix

def matrix_to_base16(matrix, original_length):
    """
    Convert the harmonic-compressed matrix back to Base-16 representation.
    """
    flattened_data = ''.join(str(int(round(bit))) for bit in matrix.flatten()[:original_length])
    base16_data = ''.join(
        format(int(flattened_data[i:i+4], 2), 'x') for i in range(0, len(flattened_data), 4)
    )
    return base16_data

# Harmonic Compression and Expansion
def harmonize_data(data, harmonic_constant, compress=True):
    """
    Compress or expand data by aligning with the harmonic constant.
    """
    data = np.copy(data)
    gain = 1.0 if compress else -1.0

    for _ in tqdm(range(MAX_ITERATIONS), desc="Harmonic adjustment"):
        delta = np.mean(data) - harmonic_constant
        adjustment = delta * gain
        data -= adjustment
        data = np.clip(data, 0, 1)  # Ensure binary range
        if abs(delta) < TOLERANCE:
            break

    return data

# Full Compression and Decompression Workflow
def hybrid_file_compression(input_file, compressed_file, restored_file):
    """
    Compress and decompress any file using harmonic compression.
    """
    # Step 1: Encode File as Base-16
    base16_data = file_to_base16(input_file)
    print(f"Base-16 Encoded Data Length: {len(base16_data)}")

    # Step 2: Convert Base-16 to Matrix
    data_matrix = base16_to_matrix(base16_data)
    print(f"Matrix Dimensions: {data_matrix.shape}")

    # Step 3: Harmonic Compression
    compressed_matrix = harmonize_data(data_matrix, HARMONIC_CONSTANT, compress=True)
    compressed_size = np.count_nonzero(compressed_matrix)
    print(f"Compressed Data Non-Zero Elements: {compressed_size}")

    # Step 4: Save Compressed Data
    np.savez_compressed(compressed_file, compressed_matrix=compressed_matrix)
    print(f"Compressed data saved to {compressed_file}")

    # Step 5: Load Compressed Data and Expand
    loaded_data = np.load(compressed_file)['compressed_matrix']
    expanded_matrix = harmonize_data(loaded_data, HARMONIC_CONSTANT, compress=False)

    # Step 6: Convert Matrix Back to Base-16
    restored_base16 = matrix_to_base16(expanded_matrix, len(base16_data))

    # Step 7: Decode Base-16 to Original File
    base16_to_file(restored_base16, restored_file)

    # Compression Statistics
    compression_ratio = compressed_size / (len(base16_data) * 4)
    print(f"Compression Ratio: {compression_ratio:.4f}")
    print(f"Compression Successful. Restored file: {restored_file}")

# Main Execution
if __name__ == "__main__":
    # Input and Output Files
    input_file = "d:\\test.wav"           # Input file path (any file)
    compressed_file = "d:\\compressed_data.npz"  # Compressed file
    restored_file = "d:\\restored_test.wav"      # Restored output file

    # Execute compression and decompression
    hybrid_file_compression(input_file, compressed_file, restored_file)
