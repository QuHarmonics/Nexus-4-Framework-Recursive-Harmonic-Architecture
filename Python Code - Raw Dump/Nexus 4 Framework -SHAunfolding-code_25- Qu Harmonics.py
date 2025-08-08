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
    if any(int(char) >= current_base for char in hash_data):
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
    
    # Combine and return the reconstructed data
    return "".join(expanded_container)


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


# Example Usage
hash_data = "123456788"  # Example valid hash in base-9
current_base = 9         # Current base of the hash
previous_base_style = 8  # Base style to reconstruct in

try:
    reconstructed_data = unhash(hash_data, current_base, previous_base_style)
    print("Reconstructed Data:", reconstructed_data)
except ValueError as e:
    print(e)
