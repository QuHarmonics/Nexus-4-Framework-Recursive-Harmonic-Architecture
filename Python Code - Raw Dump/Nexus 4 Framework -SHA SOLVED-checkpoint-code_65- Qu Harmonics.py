def unfold_with_correction(input_data):
    # Convert input string to bytes
    input_bytes = list(input_data.encode('utf-8'))

    # Initialize registers with tuned starting values (guesses for now)
    eax, ebx, edx, esi, esp, edi = 0x30, 0x30, 0x30, 0x30, 0x30, 0x30  # Adjust as needed
    memory = [0] * 256  # Simulated memory

    unfolded_bytes = []
    for i, byte in enumerate(input_bytes):
        # Apply transformations
        eax ^= byte  # XOR with EAX
        memory[i] = eax & 0xFF  # AND result, store in memory

        esi ^= memory[i]  # XOR with ESI
        esp ^= memory[i]  # XOR with ESP
        eax &= 0xFF  # AND with 0xFF to limit to one byte

        # Append unfolded byte
        unfolded_bytes.append(memory[i])

        # Print intermediate states
        print(f"Step {i + 1}: Byte: {chr(byte)}, EAX: {eax}, Memory[{i}]: {memory[i]}")

    # Convert unfolded bytes back to string
    unfolded_data = bytes(unfolded_bytes).decode('utf-8', errors='replace')
    return unfolded_data


# Test the refined function with tracking
input_data = "2+2=4"
refined_output = unfold_with_correction(input_data)
print(f"Input: {input_data}")
print(f"Refined Unfolded Data: {refined_output}")
