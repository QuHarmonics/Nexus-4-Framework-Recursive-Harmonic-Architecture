import numpy as np

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
    return [int(bit) for bit in binary.zfill(len(data) * int(np.log2(base)))]  # Ensure proper bit length

def compress_with_fixed_container(data, initial_base=16, final_base=17, container_size=265):
    """
    Perform compression using a fixed container and iterative base conversions.
    """
    # Pad data to fit the initial container size
    container = [0] * container_size
    data_bits = [int(bit) for bit in "".join(f"{ord(c):08b}" for c in data)]  # Convert seed to binary
    container[:len(data_bits)] = data_bits  # Fill the initial container with data

    print("Initial Container (Binary):", container)

    # Iteratively perform base conversions
    current_data = container
    for base in range(initial_base, final_base + 1):
        # Step 1: Convert to next base
        converted_base = binary_to_base(current_data, base)
        print(f"Base-{base} Representation:", converted_base)

        # Step 2: Convert back to binary
        current_data = base_to_binary(converted_base, base)
        print(f"Binary Representation after Base-{base} Conversion:", current_data)

        # Step 3: Store the converted data in the fixed container
        container = [0] * container_size  # Reset container
        container[:len(current_data)] = current_data[:container_size]  # Ensure it fits
        print(f"Container after Base-{base} Conversion:", container)

    return container

# Example usage
seed = "Watson, I need you"
compressed_data = compress_with_fixed_container(seed, initial_base=16, final_base=18)

print("\nFinal Compressed Data in Fixed Container:", compressed_data)
