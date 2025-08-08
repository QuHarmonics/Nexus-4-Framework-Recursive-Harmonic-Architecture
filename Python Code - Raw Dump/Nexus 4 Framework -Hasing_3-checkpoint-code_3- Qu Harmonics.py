# Define the functions as previously shared with large input testing.

import numpy as np
from math import ceil, log2

def zero_pad_to_block(seed, block_size):
    """
    Pad the binary representation of the seed with zeros to fit the block size.
    """
    seed_binary = "".join(f"{ord(c):08b}" for c in seed)
    while len(seed_binary) % block_size != 0:
        seed_binary += "0"
    return [int(bit) for bit in seed_binary], seed_binary

def binary_to_text(binary_data):
    """
    Convert binary data back to text, ensuring only complete 8-bit chunks are processed.
    """
    binary_string = "".join(map(str, binary_data))
    binary_string = binary_string[:len(binary_string) - (len(binary_string) % 8)]  # Ensure valid 8-bit chunks
    chars = [chr(int(binary_string[i:i+8], 2)) for i in range(0, len(binary_string), 8)]
    return "".join(chars).rstrip("\x00")

def binary_to_hex(binary_data):
    """
    Convert binary data to hexadecimal representation.
    """
    binary_string = "".join(map(str, binary_data))
    hex_string = hex(int(binary_string, 2))[2:].upper()
    return hex_string

def compress_large_input(seed, block_size=512):
    """
    Compress input dynamically by padding and processing binary data.
    """
    padded_binary, padded_string = zero_pad_to_block(seed, block_size)
    print(f"Original Binary Representation (Without Padding): {''.join(f'{ord(c):08b}' for c in seed[:64])}...")
    print(f"Original Binary Representation (With Padding): {padded_string[:256]}... (truncated)")

    print("\nPadded Binary Data:")
    print(f"{padded_binary[:64]}... (truncated)")

    # Calculate and display hash
    hash_hex = binary_to_hex(padded_binary)
    print("\nHash (Hexadecimal Representation):")
    print(hash_hex[:64] + "... (truncated)")

    return padded_binary, hash_hex

def decompress_large_input(binary_data):
    """
    Decompress binary data and convert it back to text.
    """
    reconstructed_text = binary_to_text(binary_data)
    print("\nReconstructed Text from Binary (First 64 Characters):")
    print(reconstructed_text[:64] + "... (truncated)")
    return reconstructed_text

# Generate a large seed of 3000 characters for testing
seed = "a" * 3000
print("Original Text (First 64 Characters):", seed[:64] + "... (truncated)")

# Compress the input
compressed_binary, hash_hex = compress_large_input(seed)

# Decompress the binary back to text
decompressed_text = decompress_large_input(compressed_binary)

# Validate the decompression result matches the original input
validation_result = "Success" if decompressed_text == seed else "Failure"
print("\nValidation Result:", validation_result)
