def simulate_hash_operations(input_data):
    # Convert the input string to a list of bytes
    input_bytes = list(input_data.encode('utf-8'))

    # Initialize registers and memory
    eax, ebx, edx, esi, esp, edi = 0, 0, 0, 0, 0, 0
    memory = [0] * 256  # Simulated memory

    # Example: Perform cyclic XOR and AND operations
    for i, byte in enumerate(input_bytes):
        eax ^= byte  # XOR with EAX
        memory[i] = eax & 0xFF  # AND result, store in memory

        esi ^= memory[i]  # XOR with ESI
        esp ^= memory[i]  # XOR with ESP
        eax &= 0xFF  # Apply AND to EAX

    # Construct the "unfolded" data
    unfolded_data = bytes(memory[:len(input_bytes)]).decode('utf-8', errors='replace')
    return unfolded_data


# Test the function
input_data = "2+2=4"
unfolded = simulate_hash_operations(input_data)
print(f"Input: {input_data}")
print(f"Unfolded Data: {unfolded}")
