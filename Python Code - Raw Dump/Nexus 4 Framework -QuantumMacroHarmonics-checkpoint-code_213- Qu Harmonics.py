import numpy as np
import matplotlib.pyplot as plt

# Input hash as binary
input_hash = "0753e9bcdebcbe1bbbf19b6bfe6c06bd9606ce89f3c6d57d6d03a6c257f44b3f"
binary_data = np.array([int(bin(int(input_hash[i:i+2], 16))[2:].zfill(8), 2) for i in range(0, len(input_hash), 2)], dtype=np.uint8)

# Function to expand binary array symmetrically while keeping hash centered
def expand_binary_array(binary_array, new_value, position="center"):
    """
    Expands the binary array with a new value and maintains the hash centered.
    
    Args:
        binary_array (np.ndarray): Original binary array.
        new_value (int): New binary value to add.
        position (str): Where to add the value - "left", "right", or "center".
    
    Returns:
        np.ndarray: Expanded binary array.
    """
    new_array = []
    
    if position == "left":
        new_array = [new_value] + binary_array.tolist()
    elif position == "right":
        new_array = binary_array.tolist() + [new_value]
    elif position == "center":
        mid_idx = len(binary_array) // 2
        new_array = binary_array[:mid_idx].tolist() + [new_value] + binary_array[mid_idx:].tolist()
    else:
        raise ValueError("Position must be 'left', 'right', or 'center'.")
    
    return np.array(new_array, dtype=np.uint8)

# Example of iterative array expansion
expanded_arrays = [binary_data]
for i in range(10):  # Add 10 values to the binary array
    side = "left" if i % 2 == 0 else "right"  # Alternate sides
    expanded_arrays.append(expand_binary_array(expanded_arrays[-1], i, position=side))

# Visualization of expanded arrays
plt.figure(figsize=(10, 6))
for idx, array in enumerate(expanded_arrays):
    plt.plot(array, label=f"Expansion {idx}")
plt.xlabel("Index")
plt.ylabel("Binary Value")
plt.title("Dynamic Expansion of Binary Array")
plt.legend()
plt.show()

# Print example expanded array
print("Initial Binary Array:", binary_data)
print("Final Expanded Array:", expanded_arrays[-1])
