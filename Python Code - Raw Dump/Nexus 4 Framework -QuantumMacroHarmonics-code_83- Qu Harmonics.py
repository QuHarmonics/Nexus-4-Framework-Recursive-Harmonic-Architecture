import numpy as np
import matplotlib.pyplot as plt

# Constants for Mark1 framework
EXPANSION_FACTOR = 1.5

# Function to convert a hash value to its binary representation
def hash_to_binary(hash_value):
    """
    Converts a hexadecimal hash into its binary representation.
    """
    return ''.join(format(int(char, 16), '04b') for char in hash_value)

# Step 1: Encode the Binary Data into H (Storage)
def store_in_H(binary_data, expansion_factor=EXPANSION_FACTOR):
    """
    Encodes binary data into H array with cumulative storage.
    """
    harmonics = np.cumsum(binary_data.astype(np.float64) * expansion_factor)
    return harmonics

# Step 2: Reverse the Process (Retrieve Original Data)
def retrieve_from_H(harmonics, expansion_factor=EXPANSION_FACTOR):
    """
    Retrieves binary data from H array.
    """
    first_value = harmonics[0] / expansion_factor
    reversed_data = np.diff(harmonics) / expansion_factor
    reversed_data = np.insert(reversed_data, 0, first_value)
    return np.round(reversed_data).astype(np.uint8)

# Input hash value
hash_value = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"

# Convert hash to binary representation
binary_string = hash_to_binary(hash_value)
binary_data = np.array([int(bit) for bit in binary_string], dtype=np.uint8)

# Store the binary data into H
harmonics = store_in_H(binary_data)

# Retrieve the binary data from H
retrieved_data = retrieve_from_H(harmonics)

# Visualize H
def visualize_harmonics(harmonics):
    plt.figure(figsize=(10, 6))
    plt.plot(harmonics, label="H(n)", color='blue')
    plt.title("Harmonic Storage Representation")
    plt.xlabel("Index")
    plt.ylabel("H(n)")
    plt.legend()
    plt.grid()
    plt.show()

visualize_harmonics(harmonics)

# Validate the process
original_binary_string = ''.join(map(str, binary_data))
retrieved_binary_string = ''.join(map(str, retrieved_data))
print("Original Binary Data (First 100 bits):", original_binary_string[:100])
print("Retrieved Binary Data (First 100 bits):", retrieved_binary_string[:100])

if np.array_equal(binary_data, retrieved_data):
    print("Data matches successfully!")
else:
    print("Data mismatch detected.")
