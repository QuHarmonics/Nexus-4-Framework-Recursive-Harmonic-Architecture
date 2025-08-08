def zero_pad_to_block(seed, container_size):
    """
    Pad the binary representation of the seed with zeros to fit the container size.
    """
    # Convert seed to binary representation
    seed_binary = "".join(f"{ord(c):08b}" for c in seed)

    # Ensure the length matches the container size by adding zeros
    if len(seed_binary) > container_size:
        raise ValueError("Seed data exceeds the fixed container size.")
    
    padding_needed = container_size - len(seed_binary)
    padded_binary = seed_binary + "0" * padding_needed

    return [int(bit) for bit in padded_binary], len(seed_binary)  # Return padded data and original length

def compress_large_input(seed, container_size=1024, initial_base=16, final_base=18):
    """
    Compress large input using a fixed container size.
    """
    # Convert the seed to binary and pad to the container size
    padded_data, original_length = zero_pad_to_block(seed, container_size)
    
    print(f"Input Length: {original_length}, Container Size: {container_size}")

    # Initialize the container
    container = [0] * container_size
    container[:len(padded_data)] = padded_data  # Fill the container

    print("\nInitial Container (Binary):", container)

    # Iteratively perform base conversions
    current_data = container
    for b in range(initial_base, final_base + 1):
        # Step 1: Convert to next base
        converted_base = binary_to_base(current_data, b)
        print(f"\nBase-{b} Representation:", converted_base)

        # Step 2: Convert back to binary
        current_data = base_to_binary(converted_base, b)
        print(f"Binary Representation after Base-{b} Conversion:", current_data)

        # Step 3: Store the converted data in the container
        container = [0] * container_size  # Reset container
        container[:len(current_data)] = current_data[:container_size]  # Ensure it fits
        print(f"Container after Base-{b} Conversion:", container)

    # Convert final binary container to hexadecimal
    final_hex = binary_to_hex(container)
    print("\nFinal Compressed Data in Hexadecimal:", final_hex)

    return container, final_hex, seed

# Example usage
seed = "Watson, I need you"
container_size = 8192
compressed_data, compressed_hex, adjusted_seed = compress_large_input(seed, container_size)

print("\nFinal Hash in Hexadecimal:", compressed_hex)
