import math

def fold_data_with_metadata(binary_data, block_size):
    """
    Folds binary data into a hex representation with embedded metadata.
    The block size determines the folding structure and is included in the output.
    """
    # Pad binary data to match the block size
    padding_length = (block_size - (len(binary_data) % block_size)) % block_size
    padded_binary = binary_data + '0' * padding_length

    # Split binary data into blocks and convert each to hex
    blocks = [padded_binary[i:i + block_size] for i in range(0, len(padded_binary), block_size)]
    folded_hex = ''.join(hex(int(block, 2))[2:].zfill(block_size // 4).upper() for block in blocks)

    # Embed metadata (block size and padding length)
    metadata = f"{block_size:04X}{padding_length:04X}"  # Block size and padding in 4-digit hex
    return folded_hex + metadata

def unfold_data_with_metadata(folded_hex, block_size=None):
    """
    Unfolds hex data back into binary using embedded metadata.
    If block_size is not provided, it will extract it from the metadata.
    """
    # Extract metadata
    metadata_start = -8
    metadata = folded_hex[metadata_start:]
    folded_hex = folded_hex[:metadata_start]

    # Decode block size and padding from metadata
    block_size = block_size or int(metadata[:4], 16)
    padding_length = int(metadata[4:], 16)

    # Split folded hex into blocks and convert each back to binary
    binary_blocks = [
        bin(int(folded_hex[i:i + block_size // 4], 16))[2:].zfill(block_size)
        for i in range(0, len(folded_hex), block_size // 4)
    ]
    unfolded_binary = ''.join(binary_blocks)

    # Remove padding
    unfolded_binary = unfolded_binary[:len(unfolded_binary) - padding_length]

    return unfolded_binary

def test_data_folding():
    original_binary = '110100101011001010101001110100101011001010101001110100101011001010101001110100101011001010101001110100101011001010101001110100101011001010101001110100101011001010101001110100101011001010101001110100101011001010101001110100101011001010101001110100101011001010101001'  # Example binary input
    block_size = 32  # Exponential block size

    print(f"Original Binary: {original_binary}")
    folded_hex = fold_data_with_metadata(original_binary, block_size)
    print(f"Folded Hex: {folded_hex}")

    unfolded_binary = unfold_data_with_metadata(folded_hex)
    print(f"Unfolded Binary: {unfolded_binary}")
    print(f"Unfolded Binary Matches Original: {unfolded_binary == original_binary}")

# Run the test
test_data_folding()
