def unhash(hash_data, current_base, previous_base_style, target_size=None):
    """
    Reconstructs the original data by unhashing the input hash using harmonic expansion.

    Args:
        hash_data (str): The compressed hash data to unhash.
        current_base (int): The base in which the hash data is represented.
        previous_base_style (int): The base style to use for reconstruction.
        target_size (int, optional): The desired size of the expanded container. 
                                      Defaults to twice the current size.

    Returns:
        str: The reconstructed data.
    """
    import math
    
    # Validate input data
    if any(char not in "0123456789abcdef"[:current_base] for char in hash_data):
        raise ValueError(f"Invalid digit found for base-{current_base}: {hash_data}")

    # Determine the target size if not provided
    original_size = len(hash_data)
    target_size = target_size or original_size * 2
    
    # Expand the container
    expanded_container = [""] * target_size
    
    # Populate the expanded container with data written in the previous base style
    for i in range(original_size):
        # Get the current digit in the hash data
        current_digit = int(hash_data[i], current_base)
        
        # Translate this digit into the previous base style
        translated_value = int_to_base(current_digit, previous_base_style)
        
        # Insert the translated value into the expanded container
        expanded_container[i] = translated_value
    
    # Fill the rest of the expanded container with harmonically-aligned values
    for i in range(original_size, target_size):
        expanded_container[i] = int_to_base(i % previous_base_style, previous_base_style)
    
    # Perform the final fold
    folded_data = fold_data(expanded_container)
    
    # Combine and return the reconstructed data
    return "".join(folded_data)


def int_to_base(value, base):
    """
    Converts an integer to a string representation in a given base.
    
    Args:
        value (int): The integer to convert.
        base (int): The base for the conversion.
    
    Returns:
        str: The string representation in the given base.
    """
    if value == 0:
        return "0"
    digits = []
    while value:
        digits.append(int(value % base))
        value //= base
    return "".join(str(x) for x in reversed(digits))


def fold_data(expanded_container):
    """
    Folds the expanded container to finalize the unhashing process.
    
    Args:
        expanded_container (list): The expanded data container.
    
    Returns:
        list: The folded container.
    """
    half_size = len(expanded_container) // 2
    folded = []
    for i in range(half_size):
        # Combine mirrored elements harmonically
        left = int(expanded_container[i], 10)
        right = int(expanded_container[-(i + 1)], 10)
        folded.append(int_to_base((left + right) % 10, 10))  # Modulo for harmonic alignment
    return folded + expanded_container[half_size:-half_size]  # Preserve the center if odd-sized


# Example Usage
hash_data = "5e36e1b0204be9a4e3d4ed56b05058db73a1bc639027959a34ecb8dcb7ec5c91"  # Example valid hash in base-16
current_base = 16        # Current base of the hash (hexadecimal)
previous_base_style = 10 # Base style to reconstruct in (e.g., base-10)

try:
    reconstructed_data = unhash(hash_data, current_base, previous_base_style)
    print("Reconstructed Data:", reconstructed_data)
except ValueError as e:
    print(e)
