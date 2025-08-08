def iterative_correction(input_data, target_data, target_length=None):
    """
    Unfolds and harmonizes the transformation of input_data to match target_data.

    Args:
        input_data (str): The hashed input string (e.g., hash-like value).
        target_data (str): The pre-hash value of the seed with padding.
        target_length (int, optional): Length to pad target_data to match input_data.

    Returns:
        str: The unfolded and harmonized string.
    """
    # Convert input data (hash) to bytes
    input_bytes = list(bytearray.fromhex(input_data))

    # Convert target data to padded bytes
    target_bytes = list(target_data.encode('utf-8'))

    # If the target needs to be padded, repeat or zero-pad it
    if target_length is None:
        target_length = len(input_bytes)

    # Extend or trim the target_bytes to match input length
    if len(target_bytes) < target_length:
        target_bytes += [0] * (target_length - len(target_bytes))
    elif len(target_bytes) > target_length:
        target_bytes = target_bytes[:target_length]

    # Initialize registers (guess initial values)
    eax, ebx, edx, esi, esp, edi = 0x30, 0x30, 0x30, 0x30, 0x30, 0x30
    memory = [0] * 256  # Simulated memory

    unfolded_bytes = []
    for i, byte in enumerate(input_bytes):
        # Apply transformations
        eax ^= byte  # XOR with EAX
        memory[i] = eax & 0xFF  # AND result, store in memory

        esi ^= memory[i]  # XOR with ESI
        esp ^= memory[i]  # XOR with ESP
        eax &= 0xFF  # AND with 0xFF to limit to one byte

        # Check against the target byte
        target_byte = target_bytes[i]
        if memory[i] != target_byte:
            # Adjust EAX to reduce mismatch
            correction = target_byte ^ memory[i]
            eax ^= correction
            memory[i] = eax & 0xFF

        unfolded_bytes.append(memory[i])

        # Print intermediate states
        print(f"Step {i + 1}: Byte: {byte}, Target: {target_byte}, "
              f"EAX: {eax}, Memory[{i}]: {memory[i]}")

    # Convert unfolded bytes back to string
    unfolded_data = bytes(unfolded_bytes).decode('utf-8', errors='replace')
    return unfolded_data


# Test the iterative correction function
input_data = "e52d388f44d390cd721c3020c00b1d34ec08cb97a844abc7ca076522d7aaf541"  # Example hash
target_data = "2+2=4"  # Pre-hash value of the seed
refined_output = iterative_correction(input_data, target_data, target_length=len(input_data) // 2)
print(f"Input: {input_data}")
print(f"Target: {target_data}")
print(f"Final Unfolded Data: {refined_output}")
