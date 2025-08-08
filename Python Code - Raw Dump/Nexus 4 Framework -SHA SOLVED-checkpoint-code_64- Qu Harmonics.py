def refine_unfolding(input_data):
    # Convert input string to bytes
    input_bytes = list(input_data.encode('utf-8'))

    # Initialize registers with guessed starting values
    eax, ebx, edx, esi, esp, edi = 0x32, 0x20, 0x20, 0x20, 0x20, 0x20  # Example guesses
    memory = [0] * 256  # Simulated memory

    # Simulate unfolding process
    unfolded_bytes = []
    for i, byte in enumerate(input_bytes):
        eax ^= byte  # XOR with EAX
        memory[i] = eax & 0xFF  # AND result, store in memory

        # Simulate other operations
        esi ^= memory[i]
        esp ^= memory[i]
        eax &= 0xFF

        # Append the unfolded byte
        unfolded_bytes.append(memory[i])

    # Convert to string
    unfolded_data = bytes(unfolded_bytes).decode('utf-8', errors='replace')
    return unfolded_data


# Test the refined function
input_data = "2+2=4"
refined_output = refine_unfolding(input_data)
print(f"Input: {input_data}")
print(f"Refined Unfolded Data: {refined_output}")
