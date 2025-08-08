def H_transform(base10_data):
    """
    Perform the H transformation: Convert Base 10 data to binary, process into Base 9, and store in a Base 11 container.
    
    Args:
        base10_data (list of int): Input data in Base 10 format.
        
    Returns:
        dict: A dictionary containing the Base 11 array, intermediate Base 9 results, and binary states.
    """
    # Step 1: Base 10 → Fixed 4-bit binary strings
    binary_states = [format(num, '04b') for num in base10_data]
    
    # Step 2: Binary string → Base 9 value (mod 9)
    base9_results = [int(b, 2) % 9 for b in binary_states]
    
    # Step 3: Base 9 values → Fixed 3-bit binary
    base10_reconstruction = [list(map(int, format(num, '03b'))) for num in base9_results]
    
    # Step 4: Flatten bits to create Base 11 array
    flat_bits = [bit for bits in base10_reconstruction for bit in bits]
    
    # Step 5: Pad to match Base 11 container size (2x original)
    padding_size = (2 * len(base10_data)) - len(flat_bits)
    base11_container = flat_bits + [""] * padding_size
    
    return {
        "Base 11 Container": base11_container,
        "Binary States": binary_states,
        "Base 9 Results": base9_results,
        "Base 10 Reconstruction": base10_reconstruction,
    }

# Example Usage
base10_data = [7, 2, 5, 8]
result = H_transform(base10_data)

# Print result nicely
import pprint
pprint.pprint(result)
