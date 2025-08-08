import numpy as np

# Initialize the hash
input_hash = "0753e9bcdebcbe1bbbf19b6bfe6c06bd9606ce89f3c6d57d6d03a6c257f44b3f"  # Replace with your SHA-256 hash
binary_hash = [int(bin(int(c, 16))[2:].zfill(4)) for c in input_hash]  # Convert hash to binary representation
hash_array = np.array([int(bit) for group in binary_hash for bit in str(group)], dtype=np.uint8)

# Initial pointers
current_position = len(hash_array) // 2 + 1  # Start at mid-point + 1
current_hash_pointer = 1
exp_hash_pointer = 1
center_pointer = 0

# Parameters for iterations
iterations = 512  # Set this to control the extent of the growth
adjustments = []  # Store adjustments for debugging

def adjust_array(hash_array, insert_position, value):
    """Insert a value into the hash array."""
    new_array = np.insert(hash_array, insert_position, value)
    return new_array

# Loop through the steps
for _ in range(iterations):
    insert_position = current_position - (current_hash_pointer ** 2)
    if insert_position < 0 or insert_position >= len(hash_array):
        break  # Stop if we go out of bounds

    # Insert into array based on even or odd condition
    if (current_hash_pointer ** 2) % 2 == 0:
        hash_array = adjust_array(hash_array, insert_position, 1)
    else:
        hash_array = adjust_array(hash_array, insert_position, 0)

    # Store adjustment for debugging
    adjustments.append((current_position, insert_position, current_hash_pointer ** 2))

    # Shift current position and pointers
    current_position += 1
    current_hash_pointer += 1
    exp_hash_pointer *= 2
    center_pointer -= 2

# Visualization
print("Initial Hash (binary):", ''.join(map(str, hash_array[:256])))  # Show the original binary hash
print("Expanded Hash Array:", ''.join(map(str, hash_array[:512])))  # Show the expanded binary array
print("Adjustments Made:", adjustments[:10])  # Print the first 10 adjustments for debugging
