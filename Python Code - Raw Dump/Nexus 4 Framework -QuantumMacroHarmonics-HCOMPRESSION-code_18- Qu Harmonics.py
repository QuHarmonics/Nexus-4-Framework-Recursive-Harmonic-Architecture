import hashlib
import numpy as np
import matplotlib.pyplot as plt

# Function to convert a string to binary
def string_to_binary(string):
    return ''.join(format(ord(c), '08b') for c in string)

# Function to hash a string and return binary
def hash_to_binary(string):
    hash_hex = hashlib.sha256(string.encode()).hexdigest()
    return ''.join(format(int(c, 16), '04b') for c in hash_hex)

# Simulated H Analyzer function
def H_analyzer(binary_data):
    iterations = np.arange(len(binary_data))
    H_values = np.sin(2 * np.pi * iterations / len(binary_data)) * (binary_data - 0.5)
    return iterations, H_values

# Generate original padded binary
original_binary = string_to_binary("abc") + '1'  # Add single 1 for padding
remaining_zeros = 512 - len(original_binary) - 64  # Calculate remaining padding
original_padded = original_binary + '0' * remaining_zeros + format(len("abc") * 8, '064b')  # Append length in 64 bits
padded_binary = [int(b) for b in original_padded]

# Hash binary for "abc"
hashed_binary = [int(b) for b in hash_to_binary("abc")]

# Convert to numpy arrays for processing
hashed = np.array(hashed_binary)
original_padded_np = np.array(padded_binary)

# Append zeros iteratively and analyze
max_iterations = 260
converged = False

for i in range(max_iterations):
    hashed = np.append(hashed, 0)  # Append one zero
    iterations, H_hashed = H_analyzer(hashed)
    _, H_original = H_analyzer(original_padded_np[:len(hashed)])  # Match lengths of data

    # Plot the waveform with the original padded binary overlay
    if i % 2 == 0 or i == max_iterations - 1:
        plt.figure(figsize=(12, 6))
        plt.plot(iterations, H_hashed, label=f"Hashed Binary (Iteration {i+1})", color='blue')
        plt.plot(iterations, H_original, label="Original Padded Binary", color='red', linestyle='dashed')
        plt.xlabel("Iteration (n)")
        plt.ylabel("H(n)")
        plt.title("H Analyzer Output - Appending Zeros with Original Binary Overlay")
        plt.legend()
        plt.grid()
        plt.show()

    # Check if the waveform matches the length of the original padded binary
    if len(H_hashed) == len(original_padded_np):
        print(f"Converged at iteration {i+1}.")
        converged = True
        
        # Convert H_hashed back into binary
        binary_reconstructed = [1 if h > 0 else 0 for h in H_hashed]
        print("Reconstructed Binary from H:")
        print(''.join(map(str, binary_reconstructed[:512])))  # Display first 512 bits for clarity
        break

if not converged:
    print("Did not converge within the maximum number of iterations.")
