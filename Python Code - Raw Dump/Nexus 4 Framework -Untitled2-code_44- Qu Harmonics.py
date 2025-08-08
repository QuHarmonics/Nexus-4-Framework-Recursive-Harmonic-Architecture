import numpy as np
from math import ceil, log2

def zero_pad_to_block(seed, block_size=1024):
    """
    Pad the binary representation of the seed with zeros to fit the block size.
    """
    # Convert seed to binary representation
    seed_binary = "".join(f"{ord(c):08b}" for c in seed)

    # Pad with zeros to fit the block size
    padded_binary = seed_binary
    while len(padded_binary) % block_size != 0:
        padded_binary += "0"

    return [int(bit) for bit in padded_binary]  # Return as a list of bits

def determine_container_size(data_length):
    """
    Determine the smallest base and container size that fits the data length.
    """
    base = int(ceil(log2(data_length)))
    container_size = 2 ** base
    return base, container_size

def binary_to_base(data, base):
    """
    Convert binary data to a specified base.
    """
    value = int("".join(str(bit) for bit in data), 2)  # Convert binary to integer
    result = []
    while value > 0:
        result.append(value % base)
        value //= base
    return result[::-1]  # Return reversed for proper order

def base_to_binary(data, base):
    """
    Convert data in a specified base back to binary.
    """
    value = sum(d * (base ** i) for i, d in enumerate(reversed(data)))  # Convert base to integer
    binary = bin(value)[2:]  # Convert integer to binary string
    return [int(bit) for bit in binary.zfill(len(data) * int(log2(base)))]  # Ensure proper bit length

def binary_to_hex(binary_data):
    """
    Convert binary data to hexadecimal representation.
    """
    binary_string = "".join(map(str, binary_data))
    hex_string = hex(int(binary_string, 2))[2:].upper()  # Convert to hex and remove "0x" prefix
    return hex_string

def binary_to_text(binary_data):
    """
    Convert binary data back to text.
    """
    # Strip leading zeros
    binary_string = "".join(map(str, binary_data)).lstrip("0")

    # Convert binary string to text
    chars = [chr(int(binary_string[i:i+8], 2)) for i in range(0, len(binary_string), 8)]
    return "".join(chars).rstrip("\x00")  # Strip null characters added during padding

def compress_large_input(seed, initial_base=16, final_base=18):
    """
    Compress large input dynamically by adjusting container size and base.
    """
    # Remove the first zero and comma from the seed
    seed = seed.lstrip("0").lstrip(",").strip()

    # Apply pure zero padding to the seed
    padded_data = zero_pad_to_block(seed, block_size=512)
    data_length = len(padded_data)
    base, container_size = determine_container_size(data_length)

    print(f"Input Length: {data_length}, Determined Base: {base}, Container Size: {container_size}")

    # Initialize the fixed container
    container = [0] * container_size
    container[:len(padded_data)] = padded_data[:container_size]  # Fill initial container

    print("\nInitial Container (Binary):", container)

    # Iteratively perform base conversions
    current_data = container
    for b in range(initial_base, final_base + 1):
        # Step 1: Convert to next base
        converted_base = binary_to_base(current_data, b)
        print(f"\nBase-{b} Representation:", converted_base)

        # Step 2: Convert back to binary
        current_data = base_to_binary(converted_base, b)
        print(f"Binary Representation after Base-{b} Conversion:", current_data)

        # Step 3: Store the converted data in the container
        container = [0] * container_size  # Reset container
        container[:len(current_data)] = current_data[:container_size]  # Ensure it fits
        print(f"Container after Base-{b} Conversion:", container)

    # Convert final binary container to hexadecimal
    final_hex = binary_to_hex(container)
    print("\nFinal Compressed Data in Hexadecimal:", final_hex)

    return container, final_hex, seed

def decompress_large_input(compressed_data, initial_base=18, final_base=16):
    """
    Decompress large input using dynamically adjusted container size and base.
    """
    current_data = compressed_data

    for b in range(initial_base, final_base - 1, -1):
        # Step 1: Convert binary to base
        converted_base = binary_to_base(current_data, b)
        print(f"\nBase-{b} Representation during Decompression:", converted_base)

        # Step 2: Convert base back to binary
        current_data = base_to_binary(converted_base, b)
        print(f"Binary Representation after Base-{b} Conversion:", current_data)

    # Strip trailing zeros to get the original data
    decompressed_data = current_data[:len(current_data) - current_data[::-1].index(1)]
    print("\nDecompressed Data (Binary):", decompressed_data)

    # Convert binary to text
    original_text = binary_to_text(decompressed_data)
    print("\nDecompressed Text:", original_text)

    return original_text

# Example usage
seed = "Watson, I need you"
compressed_data, compressed_hex, adjusted_seed = compress_large_input(seed)

print("\nFinal Hash in Hexadecimal:", compressed_hex)

# Decompress the data
decompressed_text = decompress_large_input(compressed_data)

# Verify that the decompressed text matches the adjusted seed
assert adjusted_seed == decompressed_text, "The decompressed text does not match the adjusted seed!"
