import numpy as np
from math import ceil, log2

def calculate_block_size(seed_length, multiplier=10):
    return 2 ** ceil(log2(seed_length)) * multiplier

def zero_pad_to_block(seed, block_size):
    seed_binary = "".join(f"{ord(c):08b}" for c in seed)
    while len(seed_binary) % block_size != 0:
        seed_binary += "0"
    return [int(bit) for bit in seed_binary], len(seed_binary.rstrip('0'))  # Return both padded and true bit length

def determine_container_size(data_length):
    base = int(ceil(log2(data_length)))
    container_size = 2 ** base
    return base, container_size

def binary_to_base(data, base):
    value = int("".join(str(bit) for bit in data), 2)
    result = []
    while value > 0:
        result.append(value % base)
        value //= base
    return result[::-1]

def base_to_binary(data, base):
    value = sum(d * (base ** i) for i, d in enumerate(reversed(data)))
    binary = bin(value)[2:]
    return [int(bit) for bit in binary.zfill(len(data) * int(log2(base)))]

def binary_to_hex(binary_data):
    binary_string = "".join(map(str, binary_data))
    return hex(int(binary_string, 2))[2:].upper()

def binary_to_text(binary_data):
    binary_string = "".join(map(str, binary_data))
    chars = [chr(int(binary_string[i:i+8], 2)) for i in range(0, len(binary_string), 8)]
    return "".join(chars).rstrip("\x00")

def compress_large_input(seed, multiplier=1, initial_base=16, final_base=18):
    original_seed = seed.strip()

    block_size = calculate_block_size(len(original_seed), multiplier)
    padded_data, true_length = zero_pad_to_block(original_seed, block_size=block_size)

    base, container_size = determine_container_size(len(padded_data))
    container = [0] * container_size
    container[:len(padded_data)] = padded_data[:container_size]

    current_data = container
    for b in range(initial_base, final_base + 1):
        converted_base = binary_to_base(current_data, b)
        current_data = base_to_binary(converted_base, b)
        container = [0] * container_size
        container[:len(current_data)] = current_data[:container_size]

    final_hex = binary_to_hex(container)
    return container, final_hex, original_seed, true_length

def decompress_large_input(compressed_data, true_length, initial_base=18, final_base=16):
    current_data = compressed_data

    for b in range(initial_base, final_base - 1, -1):
        converted_base = binary_to_base(current_data, b)
        current_data = base_to_binary(converted_base, b)

    # ✅ Truncate to original true bit length
    decompressed_data = current_data[:true_length]

    original_text = binary_to_text(decompressed_data)
    return original_text

# 🧪 Example usage
seed = "test"   # <---- Try with bigger inputs too
compressed_data, compressed_hex, adjusted_seed, true_length = compress_large_input(seed, multiplier=2)

print("\nCompressed Final Hex:", compressed_hex)

# Decompress
decompressed_text = decompress_large_input(compressed_data, true_length)

print("\nDecompressed Text:", decompressed_text)

# Verify
