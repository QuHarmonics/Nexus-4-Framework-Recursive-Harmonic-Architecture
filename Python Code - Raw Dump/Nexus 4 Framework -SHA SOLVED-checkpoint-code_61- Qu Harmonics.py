import hashlib

def zero_pad_input(input_string, block_size=64):
    """
    Pads the input string with zeros to align it to the desired block size.

    Args:
    - input_string: The original input string to be hashed.
    - block_size: The block size (in bytes) to pad the input to. Default is 64 bytes (512 bits for SHA).

    Returns:
    - A padded byte array of the input string.
    """
    input_bytes = input_string.encode('utf-8')
    padding_length = block_size - (len(input_bytes) % block_size)
    padded_bytes = input_bytes + b'\x00' * padding_length
    return padded_bytes

def custom_sha256(input_string):
    """
    Custom SHA-256 hash generator that uses zero padding instead of standard padding.

    Args:
    - input_string: The input string to be hashed.

    Returns:
    - The hexadecimal digest of the custom SHA-256 hash.
    """
    # Zero-pad the input
    padded_data = zero_pad_input(input_string)

    # Use hashlib to hash the padded data
    sha256 = hashlib.sha256()
    sha256.update(padded_data)
    return sha256.hexdigest()

# Example Usage
if __name__ == "__main__":
    input_string = "2+2=4"
    
    # Compute the custom hash
    custom_hash = custom_sha256(input_string)
    
    print(f"Input String: {input_string}")
    print(f"Custom SHA-256 Hash with Zero Padding: {custom_hash}")
