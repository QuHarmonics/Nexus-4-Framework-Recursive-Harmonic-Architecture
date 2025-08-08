import hashlib
import itertools

def hash_string(input_string):
    """Hash the input string using SHA-256 and return the hexadecimal digest."""
    return hashlib.sha256(input_string.encode('utf-8')).hexdigest()

def reverse_sha256(target_hash, charset, max_length):
    """
    Attempt to reverse the SHA-256 hash by iterating through all possible inputs.

    Args:
    - target_hash: The hash to reverse.
    - charset: Characters to use for generating inputs.
    - max_length: Maximum length of inputs to test.

    Returns:
    - The input string that generates the target hash, or None if not found.
    """
    for length in range(1, max_length + 1):
        for candidate in itertools.product(charset, repeat=length):
            candidate_string = ''.join(candidate)
            candidate_hash = hash_string(candidate_string)

            # Log progress for debugging (optional)
            print(f"Testing: {candidate_string} -> {candidate_hash}")

            if candidate_hash == target_hash:
                print(f"\nMatch found! Input: {candidate_string}")
                return candidate_string

    print("\nNo match found within the given constraints.")
    return None

# Example Usage
if __name__ == "__main__":
    # Target hash to reverse
    target_hash = "9e818eb02297f30a183e06b0649ca61d80d559e555399fec7717acb2ea7e170e"
    
    # Character set to test and maximum length
    charset = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ ,."  # Adjust as needed
    max_length = 30  # Limit for input length

    # Run the reverse SHA-256 process
    original_input = reverse_sha256(target_hash, charset, max_length)

    if original_input:
        print(f"Original Input: {original_input}")
    else:
        print("Failed to reverse the hash.")
