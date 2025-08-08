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
    value = 0
    for digit in data:
        value = value * base + digit

    binary_str = bin(value)[2:]  # remove '0b' prefix

    if total_bits is not None:
        # Pad at the front to keep alignment
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
    
    - We NO LONGER strip all leading zeros from the entire binary string.
    - Instead, we only remove leftover bits at the end that don't form a full byte.
    - Then convert each 8-bit chunk to one character.
    - Finally, strip trailing nulls if any.
    """
    if not binary_data:
        return ""

    # Convert the entire bit list to a string (keeping leading zeros intact!)
    binary_string = "".join(str(bit) for bit in binary_data)

    # If there's a partial byte at the END, remove it:
    remainder = len(binary_string) % 8
    if remainder != 0:
        binary_string = binary_string[:-remainder]

    # Convert each 8-bit chunk into a character
    chars = []
    for i in range(0, len(binary_string), 8):
        chunk = binary_string[i : i + 8]
        val = int(chunk, 2)
        chars.append(chr(val))

    return "".join(chars).rstrip("\x00")

def compress_large_input(seed, initial_base=16, final_base=32):
    """
    Convert 'seed' to a 512-bit (or nearest block) container, then iteratively
    re-encode only through power-of-two bases from 'initial_base' to 'final_base'.
    """
    padded_data = zero_pad_to_block(seed, block_size=512)
    data_length = len(padded_data)

    _, container_size = determine_container_size(data_length)
    container = [0] * container_size
    container[:len(padded_data)] = padded_data[:container_size]

    current_data = container
    for b in range(initial_base, final_base + 1):
        if not is_power_of_two(b):
            continue
        converted_base = binary_to_base(current_data, b)
        current_data = base_to_binary(converted_base, b, total_bits=container_size)

        new_container = [0] * container_size
        length_to_copy = min(len(current_data), container_size)
        new_container[:length_to_copy] = current_data[:length_to_copy]
        current_data = new_container

    final_hex = binary_to_hex(current_data)
    return current_data, final_hex

def decompress_large_input(compressed_data, initial_base=32, final_base=16):
    """
    Reverse the re-encoding process by iterating down from 'initial_base'
    to 'final_base' (only for bases that are powers of two).
    """
    current_data = compressed_data
    container_size = len(current_data)

    for b in range(initial_base, final_base - 1, -1):
        if not is_power_of_two(b):
            continue
        converted_base = binary_to_base(current_data, b)
        current_data = base_to_binary(converted_base, b, total_bits=container_size)

    # Remove trailing zeros up to the last '1' bit
    if 1 in current_data:
        last_one_index = len(current_data) - 1 - current_data[::-1].index(1)
        current_data = current_data[:last_one_index + 1]

    original_text = binary_to_text(current_data)
    return original_text

if __name__ == "__main__":
    seed = "Hello my friend, this worksssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssdfdfasdfasdfasdfasdfasdssss3434w34w3r"
    compressed_data, compressed_hex = compress_large_input(seed, 16, 32)
    print("Compressed Data (Binary):", compressed_data)
    print("Compressed Hex:", compressed_hex)

    decompressed_text = decompress_large_input(compressed_data, 32, 16)
    print("Decompressed Text:", decompressed_text)
