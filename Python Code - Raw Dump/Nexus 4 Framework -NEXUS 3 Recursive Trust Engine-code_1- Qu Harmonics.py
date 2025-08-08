import hashlib
import math

def calculate_harmonic_error(digest_bytes, constant=0.35):
    """
    Calculates the harmonic error of a SHA-256 digest based on its nibbles
    and a target constant.

    Args:
        digest_bytes: The 32-byte SHA-256 digest.
        constant: The target harmonic constant (default is 0.35).

    Returns:
        The L2 norm of the deviation of nibble amplitudes from the constant.
    """
    # Convert bytes to a hex string to easily access nibbles
    hex_digest = digest_bytes.hex()

    # Extract nibbles and normalize them (0-15 -> 0-1)
    nibble_amplitudes = []
    for nibble_char in hex_digest:
        # Convert hex character to integer (0-15)
        nibble_value = int(nibble_char, 16)
        # Normalize to a value between 0 and 1
        normalized_amplitude = nibble_value / 15.0
        nibble_amplitudes.append(normalized_amplitude)

    # Compute the error vector (deviation from the constant)
    error_vector = [amp - constant for amp in nibble_amplitudes]

    # Compute the L2 norm of the error vector
    l2_norm = math.sqrt(sum(e**2 for e in error_vector))

    return l2_norm

def find_harmonic_nonce(message, nonce_range_start, nonce_range_end, constant=0.35):
    """
    Searches for a nonce in a given range that minimizes the harmonic error
    of the double SHA-256 hash of the message concatenated with the nonce.

    Args:
        message: The base message bytes.
        nonce_range_start: The starting integer for the nonce search range.
        nonce_range_end: The ending integer for the nonce search range (exclusive).
        constant: The target harmonic constant (default is 0.35).

    Returns:
        A tuple containing the best nonce found and its minimum harmonic error.
    """
    best_nonce = None
    best_error = float('inf')

    for nonce_int in range(nonce_range_start, nonce_range_end):
        # Convert nonce integer to bytes (using 8 bytes for typical nonce size)
        nonce_bytes = nonce_int.to_bytes(8, byteorder='big')

        # Concatenate message and nonce
        message_with_nonce = message + nonce_bytes

        # Calculate the double SHA-256 hash
        hash1 = hashlib.sha256(message_with_nonce).digest()
        hash2 = hashlib.sha256(hash1).digest()

        # Calculate the harmonic error of the final hash
        current_error = calculate_harmonic_error(hash2, constant)

        # Check if this is the best error found so far
        if current_error < best_error:
            best_error = current_error
            best_nonce = nonce_int

    return best_nonce, best_error

# --- Example Usage ---
# Define a sample message block (can be any bytes)
sample_message = b"This is a sample message block for harmonic nonce search."

# Define a small range to search for the nonce (searching a large range can take a long time)
search_range_start = 0
search_range_end = 10000 # Search the first 10,000 nonces

print(f"Searching for harmonic nonce in range {search_range_start} to {search_range_end-1}...")

# Find the harmonic nonce
best_harmonic_nonce, min_harmonic_error = find_harmonic_nonce(
    sample_message,
    search_range_start,
    search_range_end
)

print(f"\nSearch complete.")
if best_harmonic_nonce is not None:
    print(f"Best harmonic nonce found: {best_harmonic_nonce}")
    print(f"Minimum harmonic error achieved: {min_harmonic_error}")

    # Verify the hash and error for the best nonce
    best_nonce_bytes = best_harmonic_nonce.to_bytes(8, byteorder='big')
    final_hash_for_best_nonce = hashlib.sha256(hashlib.sha256(sample_message + best_nonce_bytes).digest()).digest()
    print(f"Double SHA-256 hash for best nonce (hex): {final_hash_for_best_nonce.hex()}")
    print(f"Recalculated harmonic error for best nonce: {calculate_harmonic_error(final_hash_for_best_nonce)}")

else:
    print("No nonces found in the specified range.")

# Example of calculating error for a known digest (using the one from previous discussion)
# Note: This doesn't come from the search, just shows the error calculation function works
example_digest_hex = "d8489dddb410a2d6ad5b1a89d097bac654ea9a84f9902651a4f2085518038801"
example_digest_bytes = bytes.fromhex(example_digest_hex)
example_error = calculate_harmonic_error(example_digest_bytes)
print(f"\nHarmonic error for example digest {example_digest_hex[:10]}...: {example_error}")