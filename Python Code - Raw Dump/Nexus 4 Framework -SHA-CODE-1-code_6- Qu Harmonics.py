import numpy as np

def fold_binary(binary, block_size):
    """
    Fold binary data into blocks of the specified size and represent it as hexadecimal.
    Ensures binary is padded to a multiple of block size.
    """
    # Pad binary to ensure length is a multiple of block size
    padding_length = (block_size - (len(binary) % block_size)) % block_size
    binary = binary + "0" * padding_length  # Add trailing zeros if necessary
    
    blocks = [binary[i:i+block_size] for i in range(0, len(binary), block_size)]
    folded_hex = [hex(int(block, 2))[2:].upper().zfill(block_size // 4) for block in blocks]
    return folded_hex

def unfold_binary(folded_hex, block_size):
    """
    Unfold hexadecimal blocks back into binary.
    """
    unfolded_binary = ''.join([bin(int(h, 16))[2:].zfill(block_size) for h in folded_hex])
    return unfolded_binary.rstrip("0")  # Remove any extra padding zeros

# Original binary string
binary_data = "101010111100110111101111000010101011110011011110111100001010"

# Test with multiple block sizes
block_sizes = [4, 8, 16, 32]

results = []
for block_size in block_sizes:
    # Fold binary data into hexadecimal blocks
    folded_hex = fold_binary(binary_data, block_size)
    
    # Unfold back into binary
    unfolded_binary = unfold_binary(folded_hex, block_size)
    
    # Verify integrity
    matches_original = binary_data == unfolded_binary
    results.append((block_size, folded_hex, matches_original))

# Display results
for block_size, folded_hex, matches_original in results:
    print(f"Block Size: {block_size}")
    print(f"Folded Hex: {' '.join(folded_hex)}")
    print(f"Unfolded Binary Matches Original: {matches_original}")
    print("-")
