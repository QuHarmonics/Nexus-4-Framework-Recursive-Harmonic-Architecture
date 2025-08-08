def iterative_correction(input_data, target_data):
    input_bytes = list(input_data.encode('utf-8'))
    target_bytes = list(target_data.encode('utf-8'))

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
        print(f"Step {i + 1}: Byte: {chr(byte)}, Target: {chr(target_byte)}, "
              f"EAX: {eax}, Memory[{i}]: {memory[i]}")

    # Convert unfolded bytes back to string
    unfolded_data = bytes(unfolded_bytes).decode('utf-8', errors='replace')
    return unfolded_data


# Test the iterative correction function
input_data = "2+2=4"
target_data = "2+2=4"
refined_output = iterative_correction(input_data, target_data)
print(f"Input: {input_data}")
print(f"Target: {target_data}")
print(f"Final Unfolded Data: {refined_output}")
