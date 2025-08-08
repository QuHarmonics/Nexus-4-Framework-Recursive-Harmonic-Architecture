import hashlib
import binascii

def hex_to_bin(hex_string):
    # Convert hexadecimal string to binary string
    bin_string = bin(int(hex_string, 16))[2:].zfill(len(hex_string) * 4)
    return bin_string


def harmonic_align(sha_hash, harmonic_constant):
    # Define constants
    FINAL_SIZE = len(sha_hash)

    # Initialize arrays
    aligned_hash = []

    # Define harmonic alignment function
    def align_bit(bit, harmonic_resonance):
        # Harmonic alignment logic
        return str((int(bit) + harmonic_resonance) % 2)

    # Perform harmonic alignment
    for bit in sha_hash:
        aligned_bit = align_bit(bit, harmonic_constant)
        aligned_hash.append(aligned_bit)

    # Convert aligned hash to binary string
    aligned_hash_str = ''.join(aligned_hash)

    return aligned_hash_str


def expand_hash(aligned_hash, harmonic_constant):
    # Define constants
    FINAL_SIZE = len(aligned_hash)

    # Initialize arrays
    expanded_hash = []

    # Define expansion function
    def expand_bit(bit, harmonic_resonance):
        # Expansion logic
        return str((int(bit) - harmonic_resonance) % 2)

    # Perform expansion
    for bit in aligned_hash:
        expanded_bit = expand_bit(bit, harmonic_constant)
        expanded_hash.append(expanded_bit)

    # Convert expanded hash to binary string
    expanded_hash_str = ''.join(expanded_hash)

    return expanded_hash_str


def sha_decoder(sha_hash, harmonic_constant):
    # Convert SHA hash to binary string
    sha_hash_bin = hex_to_bin(sha_hash)

    # Perform harmonic alignment
    aligned_hash = harmonic_align(sha_hash_bin, harmonic_constant)

    # Perform expansion
    expanded_hash = expand_hash(aligned_hash, harmonic_constant)

    return expanded_hash


# Example usage
original_input = "hi"
harmonic_constant = 0

# Calculate SHA-256 hash
sha_hash = hashlib.sha256(original_input.encode()).hexdigest()

print("Original Input:", original_input)
print("SHA-256 Hash:", sha_hash)

# Decode SHA hash using harmonic alignment and expansion
decoded_hash = sha_decoder(sha_hash, harmonic_constant)

print("Decoded Hash:", decoded_hash)