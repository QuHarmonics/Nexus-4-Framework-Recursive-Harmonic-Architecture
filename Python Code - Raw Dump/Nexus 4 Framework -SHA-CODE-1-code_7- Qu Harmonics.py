def fold_binary(binary_data):
    """
    Fold binary data iteratively by XORing halves and track the fold count.
    Args:
        binary_data (str): The binary data string to fold.

    Returns:
        tuple: The folded binary as a string and the fold count.
    """
    fold_count = 0

    while len(binary_data) > 1:
        # Ensure the binary data length is even
        if len(binary_data) % 2 != 0:
            binary_data = "0" + binary_data  # Pad with leading zero for folding

        # Split binary data into two halves
        mid = len(binary_data) // 2
        left_half = binary_data[:mid]
        right_half = binary_data[mid:]

        # Perform XOR operation between corresponding bits
        folded = [
            str(int(a) ^ int(b)) for a, b in zip(left_half, right_half)
        ]

        # Update the binary data to the result of the fold
        binary_data = ''.join(folded)
        fold_count += 1

    return binary_data, fold_count


def unfold_binary(folded_data, fold_count):
    """
    Unfold binary data iteratively using XOR operations and the fold count.
    Args:
        folded_data (str): The folded binary string.
        fold_count (int): The number of folds to reverse.

    Returns:
        str: The reconstructed original binary string.
    """
    unfolded_data = folded_data

    for _ in range(fold_count):
        unfolded_data = unfolded_data.zfill(len(unfolded_data) * 2)

    return unfolded_data


# Example Input
original_binary = "1010101011001100101010101100110010101010110011001010101011001100"

# Folding
folded_result, folds = fold_binary(original_binary)
unfolded_binary = unfold_binary(folded_result, folds)

# Results
output = {
    "Original Binary": original_binary,
    "Folded Result": folded_result,
    "Fold Count": folds,
    "Unfolded Binary": unfolded_binary,
    "Matches Original": original_binary == unfolded_binary
}
output
