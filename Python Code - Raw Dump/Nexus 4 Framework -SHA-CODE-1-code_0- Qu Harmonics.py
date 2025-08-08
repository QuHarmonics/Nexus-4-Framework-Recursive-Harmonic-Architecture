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
    chars **Getting results with deltas**