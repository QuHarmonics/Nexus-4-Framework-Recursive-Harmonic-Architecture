def reverse_sha(hash_data, current_base, previous_base_style, iterations=1):
    """
    Reverses an SHA hash using linear math and partial wave observation.

    Args:
        hash_data (str): The hash data to reverse.
        current_base (int): The base of the current hash data.
        previous_base_style (int): The base style for expanding data.
        iterations (int): Number of steps to reverse the hash.

    Returns:
        str: The reconstructed wave-like data.
    """
    if iterations <= 0:
        return hash_data

    # Step 1: Partial Observation (split the data into segments)
    observed_part = hash_data[:len(hash_data) // 2]  # First half
    forward_part = hash_data[len(hash_data) // 2:]   # Second half

    # Step 2: Translate observed parts into the previous base style
    observed_part_translated = "".join(
        int_to_base(int(char, current_base), previous_base_style) for char in observed_part
    )
    forward_part_translated = "".join(
        int_to_base(int(char, current_base), previous_base_style) for char in forward_part
    )

    # Step 3: Isolate one axis in the base and harmonically align
    harmonized_data = observed_part_translated + forward_part_translated

    # Step 4: Expand, fold, and pass to the next step
    expanded_data = expand_and_fold(harmonized_data, current_base, previous_base_style)
    
    # Recursive step to process further
    return reverse_sha(expanded_data, current_base + 1, current_base, iterations - 1)


def expand_and_fold(data, current_base, previous_base_style):
    """
    Expands the data container and folds harmonically into the next base.

    Args:
        data (str): The data to expand and fold.
        current_base (int): The base of the current data.
        previous_base_style (int): The previous base style.

    Returns:
        str: The expanded and folded data.
    """
    container_size = len(data) * 2
    expanded_container = [""] * container_size

    # Fill the expanded container harmonically
    for i, char in enumerate(data):
        current_value = int(char, current_base)
        expanded_container[i] = int_to_base(current_value, previous_base_style)
    
    # Fold harmonically
    for i in range(len(data), container_size):
        expanded_container[i] = int_to_base(i % previous_base_style, previous_base_style)

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
hash_data = "5e36e1b0204be9a4e3d4ed56b05058db73a1bc639027959a34ecb8dcb7ec5c91"  # Example SHA-256 hash
current_base = 16  # Hexadecimal
previous_base_style = 10  # Base style for expansion (e.g., decimal)
iterations = 5  # Number of harmonic steps to perform

reversed_data = reverse_sha(hash_data, current_base, previous_base_style, iterations)
print("Reversed Data:", reversed_data)
