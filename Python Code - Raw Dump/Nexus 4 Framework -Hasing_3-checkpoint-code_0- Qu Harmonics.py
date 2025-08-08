import numpy as np
from math import ceil, log2

def zero_pad_to_block(seed, block_size):
    """
    Pad the binary representation of the seed with zeros to fit the block size.
    """
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

def compress_large_input(seed, initial_base=16, target_length=256, max_base=128):
    """
    Compress large input dynamically by incrementing the base
    until the hash fits the target length or the maximum base is reached.
    """
    # Convert seed to binary and pad to block size
    padded_data = zero_pad_to_block(seed, block_size=512)
    current_data = padded_data
    current_base = initial_base

    print(f"Input Length (Text): {len(seed)}, Padded Binary Length: {len(padded_data)}")

    # Iteratively increase the base to reduce binary length
    while len(current_data) > target_length:
        if current_base >= max_base:
            print("Maximum base reached without achieving target length.")
            break

        current_base += 1
        print(f"\nConverting to Base-{current_base}...")

        # Convert to the current base
        converted_base = binary_to_base(current_data, current_base)

        # Convert back to binary
        current_data = base_to_binary(converted_base, current_base)

        # Debug: Output the current binary length
        print(f"Base-{current_base}: Binary Length = {len(current_data)}")

        # If the binary data fits within the target length, stop the loop
        if len(current_data) <= target_length:
            break

    # Truncate to the target length if necessary
    if len(current_data) > target_length:
        current_data = current_data[:target_length]

    # Convert final binary to hexadecimal
    final_hex = binary_to_hex(current_data)
    print(f"\nFinal Base: {current_base}")
    print(f"Final Compressed Data in Hexadecimal: {final_hex}")

    return padded_data, current_data, final_hex


def decompress_large_input(compressed_data, initial_base, target_length=256):
    """
    Decompress the compressed data back to its binary form.
    """
    current_data = compressed_data
    current_base = initial_base

    while len(current_data) < target_length:
        current_base -= 1
        print(f"\nDecompressing from Base-{current_base}...")

        # Convert binary to the current base
        converted_base = binary_to_base(current_data, current_base)

        # Convert base back to binary
        current_data = base_to_binary(converted_base, current_base)

        # Debug: Output the current binary length
        print(f"Base-{current_base}: Binary Length = {len(current_data)}")

    return current_data


# Example usage
seed = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua."
padded_binary, compressed_data, compressed_hex = compress_large_input(seed)

# Output the padded binary and the final binary after compression
print("\nOriginal Seed Padded Binary:")
print("".join(map(str, padded_binary)))

print("\nFinal Compressed Data in Binary:")
print("".join(map(str, compressed_data)))

# Decompress the data
decompressed_binary = decompress_large_input(compressed_data, initial_base=128)

print("\nDecompressed Binary:")
print("".join(map(str, decompressed_binary)))



# Example usage
seed = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua."
padded_binary, compressed_data, compressed_hex = compress_large_input(seed)

# Output the padded binary and the final binary after compression
print("\nOriginal Seed Padded Binary:")
print("".join(map(str, padded_binary)))

print("\nFinal Compressed Data in Binary:")
print("".join(map(str, compressed_data)))

# Decompress the data
decompressed_binary = decompress_large_input(compressed_data, initial_base=128)

print("\nDecompressed Binary:")
print("".join(map(str, decompressed_binary)))



# Example usage
seed = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua."
padded_binary, compressed_data, compressed_hex = compress_large_input(seed)

# Output the padded binary and the final binary after compression
print("\nOriginal Seed Padded Binary:")
print("".join(map(str, padded_binary)))

print("\nFinal Compressed Data in Binary:")
print("".join(map(str, compressed_data)))

# Decompress the data
decompressed_binary = decompress_large_input(compressed_data, initial_base=128)

print("\nDecompressed Binary:")
print("".join(map(str, decompressed_binary)))






# Example usage
seed = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."
compressed_data, compressed_hex = compress_large_input(seed)

print("\nFinal Hash in Hexadecimal:", compressed_hex)

# Decompress the data
decompressed_text = decompress_large_input(compressed_data)

