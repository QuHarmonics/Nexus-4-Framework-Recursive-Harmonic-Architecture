import numpy as np

def pad_sha_style(seed):
    """
    Apply SHA-like padding to the input seed.
    """
    # Convert seed to binary representation
    seed_binary = "".join(f"{ord(c):08b}" for c in seed)

    # Add a single '1' bit
    padded_binary = seed_binary + "1"

    # Add '0' bits until length is 448 mod 512
    while len(padded_binary) % 512 != 448:
        padded_binary += "0"

    # Append the original seed length as a 64-bit binary
    seed_length = len(seed_binary)
    padded_binary += f"{seed_length:064b}"

    return [int(bit) for bit in padded_binary]  # Return as a list of bits

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
    return [int(bit) for bit in binary.zfill(len(data) * int(np.log2(base)))]  # Ensure proper bit length

def binary_to_hex(binary_data):
    """
    Convert binary data to hexadecimal representation.
    """
    binary_string = "".join(map(str, binary_data))
    hex_string = hex(int(binary_string, 2))[2:].upper()  # Convert to hex and remove "0x" prefix
    return hex_string

def compress_with_fixed_container(seed, initial_base=16, final_base=18, container_size=265):
    """
    Perform compression using a fixed container and iterative base conversions.
    """
    # Apply SHA-style padding to the seed
    padded_data = pad_sha_style(seed)

    print("Padded Data (Binary):", "".join(map(str, padded_data)))

    # Embed the number of binary digits (len) into the data
    binary_count = len(padded_data)
    count_binary = [int(bit) for bit in bin(binary_count)[2:].zfill(16)]  # 16-bit count
    padded_data = count_binary + padded_data

    # Initialize the fixed container
    container = [0] * container_size
    container[:len(padded_data)] = padded_data[:container_size]  # Fill initial container

    print("\nInitial Container (Binary):", container)

    # Iteratively perform base conversions
    current_data = container
    for base in range(initial_base, final_base + 1):
        # Step 1: Convert to next base
        converted_base = binary_to_base(current_data, base)
        print(f"\nBase-{base} Representation:", converted_base)

        # Step 2: Convert back to binary
        current_data = base_to_binary(converted_base, base)
        print(f"Binary Representation after Base-{base} Conversion:", current_data)

        # Step 3: Store the converted data in the fixed container
        container = [0] * container_size  # Reset container
        container[:len(current_data)] = current_data[:container_size]  # Ensure it fits
        print(f"Container after Base-{base} Conversion:", container)

    # Convert final binary container to hexadecimal
    final_hex = binary_to_hex(container)
    print("\nFinal Compressed Data in Hexadecimal:", final_hex)

    return container, final_hex

def decompress_with_fixed_container(compressed_data, initial_base=18, final_base=16):
    """
    Decompress the data using the binary count indicator.
    """
    current_data = compressed_data

    # Extract the binary count from the first 16 bits
    count_binary = current_data[:16]
    binary_count = int("".join(map(str, count_binary)), 2)
    print("\nExtracted Binary Count:", binary_count)

    # Remove the count bits from the compressed data
    current_data = current_data[16:]

    for base in range(initial_base, final_base - 1, -1):
        # Step 1: Convert binary to base
        converted_base = binary_to_base(current_data, base)
        print(f"\nBase-{base} Representation during Decompression:", converted_base)

        # Step 2: Convert base back to binary
        current_data = base_to_binary(converted_base, base)
        print(f"Binary Representation after Base-{base} Conversion:", current_data)

    # Stop at the indicated binary count
    decompressed_data = current_data[:binary_count]
    print("\nDecompressed Data (Binary):", decompressed_data)

    return decompressed_data

# Example usage
seed = "Watson, I need you"
compressed_data, compressed_hex = compress_with_fixed_container(seed, initial_base=16, final_base=18)

print("\nFinal Hash in Hexadecimal:", compressed_hex)

# Decompress the data
decompressed_data = decompress_with_fixed_container(compressed_data, initial_base=18, final_base=16)
