import hashlib
import numpy as np
import matplotlib.pyplot as plt

# Function to convert a string to its binary representation
def string_to_binary(string):
    return ''.join(format(ord(c), '08b') for c in string)

# Function to hash a string using SHA-256 and return its binary representation
def hash_to_binary(string):
    hashed = hashlib.sha256(string.encode()).hexdigest()
    return ''.join(format(int(c, 16), '04b') for c in hashed)

# Padding the binary representation of "abc"
original_binary = string_to_binary("abc") + '1'  # Add the '1' as part of the padding
padded_binary = original_binary + '0' * (512 - len(original_binary))  # Example padding to 512 bits

# Hashing "abc" and converting to binary
hashed_binary = hash_to_binary("abc")

# Simulating the H analyzer
def h_analyzer(binary_string):
    n = len(binary_string)
    x = np.linspace(0, 10, n)
    y = [1 if bit == '1' else -1 for bit in binary_string]
    return x, np.sin(x) * y

# Analyze the padded binary
x_padded, y_padded = h_analyzer(padded_binary)

# Analyze the hashed binary
x_hashed, y_hashed = h_analyzer(hashed_binary)

# Plot the results
plt.figure(figsize=(12, 6))

# Plot for padded binary
plt.subplot(1, 2, 1)
plt.plot(x_padded, y_padded, label="Padded Binary")
plt.title("H Analyzer Output - Padded Binary")
plt.xlabel("Iteration (n)")
plt.ylabel("H(n)")
plt.legend()

# Plot for hashed binary
plt.subplot(1, 2, 2)
plt.plot(x_hashed, y_hashed, label="Hashed Binary", color='orange')
plt.title("H Analyzer Output - Hashed Binary")
plt.xlabel("Iteration (n)")
plt.ylabel("H(n)")
plt.legend()

plt.tight_layout()
plt.show()
