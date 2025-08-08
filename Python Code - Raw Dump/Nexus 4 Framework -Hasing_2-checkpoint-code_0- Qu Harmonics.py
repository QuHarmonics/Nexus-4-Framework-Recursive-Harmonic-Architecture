import numpy as np
from math import ceil, log2

def is_power_of_two(n):
    """
    Checks if 'n' is a power of two.
    """
    return (n & (n - 1) == 0) and (n != 0)

def zero_pad_to_block(seed, block_size):
    """
    Convert 'seed' to binary, then pad with zeros until 'len(seed_binary)'
    is a multiple of 'block_size'.
    """
    seed_binary = "".join(f"{ord(c):08b}" for c in seed)
    while len(seed_binary) % block_size != 0:
        seed_binary += "0"
    return [int(bit) for bit in seed_binary]

def determine_container_size(data_length):
    """
    Determine the smallest power-of-two container size that can hold 'data_length' bits.
    """
    base = int(ceil(log2(data_length)))
    container_size = 2 ** base
    return base, container_size

def binary_to_base(data, base):
    """
    Interpret the entire 'data' list of bits as one integer, then convert
    that integer into an array of digits in 'base'.
    """
    # Join bits into a string, interpret as base-2 integer
    value = int("".join(str(bit) for bit in data), 2)
    if value == 0:
        return [0]
    result = []
    while value > 0:
        result.append(value % base)
        value //= base
    return result[::-1]

def base_to_binary(data, base, total_bits=None):
    """
    Convert an array of digits in 'base' back to a list of bits.
    If 'total_bits' is specified, zero-pad the front so the output
    always has 'total_bits' bits.
    """
    # Reconstruct the integer from the base-d digits
    value = 0
    for digit in data:
        value = value * base + digit

    # Convert to binary (string), removing '0b' prefix
    binary_str = bin(value)[2:]

    # Enforce a fixed total bit-length if provided: THIS is the critical update
    if total_bits is not None:
        binary_str = binary_str.zfill(total_bits)

    return [int(bit) for bit in binary_str]

def binary_to_hex(binary_data):
    """
    Convert a list of bits to a hexadecimal string.
    """
    if not binary_data:
        return ""
    binary_string = "".join(str(bit) for bit in binary_data)
    hex_string = hex(int(binary_string, 2))[2:].upper()
    return hex_string

def binary_to_text(binary_data):
    """
    Convert a list of bits to text, taking 8 bits at a time as ASCII codes.
    Leading zeros are stripped before grouping into bytes.
    """
    if not binary_data:
        return ""

    # Convert bits to a string, then strip leading zeros
    binary_string = "".join(str(bit) for bit in binary_data).lstrip("0")
    if not binary_string:
        return ""

    # Group every 8 bits into one character
    chars = []
    for i in range(0, len(binary_string), 8):
        chunk = binary_string[i:i+8]
        if len(chunk) < 8:
            break
        chars.append(chr(int(chunk, 2)))

    # Strip possible padding nulls
    return "".join(chars).rstrip("\x00")

def compress_large_input(seed, initial_base=16, final_base=32):
    """
    Convert 'seed' to a 512-bit (or nearest block) container, then iteratively
    re-encode only through power-of-two bases from 'initial_base' to 'final_base'.
    """
    # 1. Zero-pad the seed to a multiple of 512 bits.
    padded_data = zero_pad_to_block(seed, block_size=512)
    data_length = len(padded_data)

    # 2. Determine power-of-two container size.
    _, container_size = determine_container_size(data_length)

    # 3. Initialize container and copy padded data into it.
    container = [0] * container_size
    container[:len(padded_data)] = padded_data[:container_size]

    current_data = container

    # 4. Convert only through bases that are powers of two, from initial_base to final_base.
    for b in range(initial_base, final_base + 1):
        if not is_power_of_two(b):
            continue

        # a) Convert current_data (binary) -> base b
        converted_base = binary_to_base(current_data, b)

        # b) Convert base b -> binary, enforcing container_size bits
        current_data = base_to_binary(converted_base, b, total_bits=container_size)

        # c) Copy into a fresh container for consistency
        new_container = [0] * container_size
        length_to_copy = min(len(current_data), container_size)
        new_container[:length_to_copy] = current_data[:length_to_copy]
        current_data = new_container

    # 5. Render final data as hexadecimal
    final_hex = binary_to_hex(current_data)
    return current_data, final_hex

def decompress_large_input(compressed_data, initial_base=32, final_base=16):
    """
    Reverse the re-encoding process by iterating down from 'initial_base'
    to 'final_base' (only for bases that are powers of two).
    """
    current_data = compressed_data
    container_size = len(current_data)

    # Move downward from initial_base to final_base
    for b in range(initial_base, final_base - 1, -1):
        if not is_power_of_two(b):
            continue

        # a) Convert current_data (binary) -> base b
        converted_base = binary_to_base(current_data, b)

        # b) Convert base b -> binary, zero-padding to the container size
        current_data = base_to_binary(converted_base, b, total_bits=container_size)

    # Remove trailing zeros up to the last '1' bit
    if 1 in current_data:
        last_one_index = len(current_data) - 1 - current_data[::-1].index(1)
        current_data = current_data[:last_one_index + 1]

    # Convert the final bitstream to text
    original_text = binary_to_text(current_data)
    return original_text

if __name__ == "__main__":
    seed = "Hello"
    compressed_data, compressed_hex = compress_large_input(seed, 16, 32)
    print("Compressed Data (Binary):", compressed_data)
    print("Compressed Hex:", compressed_hex)

    decompressed_text = decompress_large_input(compressed_data, 32, 16)
    print("Decompressed Text:", decompressed_text)
