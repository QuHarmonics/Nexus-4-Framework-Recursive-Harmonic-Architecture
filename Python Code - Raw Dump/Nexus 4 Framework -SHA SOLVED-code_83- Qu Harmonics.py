import hashlib
import numpy as np
import matplotlib.pyplot as plt

def generate_hash(data):
    """Generate SHA-256 hash of the given data."""
    return hashlib.sha256(data.encode()).hexdigest()

def grow_from_byte1(hash_value, max_length):
    """
    Use Byte1 (the first byte) to grow the hash data recursively.
    """
    # Extract the first byte of the hash
    byte1 = hash_value[:2]
    expanded_data = bytearray([int(byte1, 16)])

    # Grow the data iteratively
    for _ in range(max_length - 1):
        new_byte = (expanded_data[-1] * 137) % 256  # Example growth logic
        expanded_data.append(new_byte)

    return expanded_data

def visualize_growth(original_hash, expanded_data):
    """Visualize the growth process."""
    fig, ax = plt.subplots(2, 1, figsize=(12, 8))

    # Convert the original hash's first 32 characters into integers
    original_bytes = [int(original_hash[i:i+2], 16) for i in range(0, 64, 2)]

    # Visualize the original hash
    ax[0].bar(range(len(original_bytes)), original_bytes, color='blue', alpha=0.6)
    ax[0].set_title("Original Hash (First 32 Bytes)")
    ax[0].set_xlabel("Byte Index")
    ax[0].set_ylabel("Value (Hexadecimal)")

    # Visualize the expanded data
    ax[1].plot(range(len(expanded_data)), expanded_data, label="Expanded Data from Byte1", color='red')
    ax[1].set_title("Expanded Data Grown from Byte1")
    ax[1].set_xlabel("Byte Index")
    ax[1].set_ylabel("Value (Decimal)")
    ax[1].legend()

    plt.tight_layout()
    plt.show()

def main():
    # Input data
    data = "Hello, Byte1 and SHA-256!"
    hash_value = generate_hash(data)

    # Expand data using Byte1 principle
    expanded_data = grow_from_byte1(hash_value, max_length=64)

    # Visualize the results
    print(f"Original Hash: {hash_value}")
    print(f"Expanded Data (First 20 Bytes): {list(expanded_data[:20])}")
    visualize_growth(hash_value, expanded_data)

if __name__ == "__main__":
    main()
