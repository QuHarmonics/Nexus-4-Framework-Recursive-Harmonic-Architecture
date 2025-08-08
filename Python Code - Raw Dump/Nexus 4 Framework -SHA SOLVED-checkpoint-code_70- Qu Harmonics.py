import numpy as np

def oversample_to_fixed_size(hex_data, target_size=512):
    """
    Oversamples hex data to produce a fixed-size output (e.g., 512 bytes).

    Args:
        hex_data (str): Input hexadecimal string representing audio-like data.
        target_size (int): Desired output size in bytes (default: 512).

    Returns:
        str: Hexadecimal string of the oversampled data.
    """
    # Convert hex string to array of integers
    input_bytes = [int(hex_data[i:i+2], 16) for i in range(0, len(hex_data), 2)]
    audio_data = np.array(input_bytes, dtype=np.uint8)
    
    # Normalize to [-1, 1] for interpolation
    normalized_data = (audio_data / 255.0) * 2 - 1

    # Determine dynamic oversampling factor
    current_size = len(normalized_data)
    factor = target_size / current_size

    # Interpolate to produce the exact target size
    original_time = np.arange(current_size)
    oversampled_time = np.linspace(0, current_size - 1, target_size)
    oversampled_audio = np.interp(oversampled_time, original_time, normalized_data)

    # Convert back to 8-bit range [0, 255]
    oversampled_bytes = ((oversampled_audio + 1) / 2 * 255).astype(np.uint8)

    # Convert to hexadecimal format
    oversampled_hex = ''.join(f'{byte:02x}' for byte in oversampled_bytes)
    return oversampled_hex

# Input data
hex_data = "18c84c4e92bc9c408bdc0738577c3ddc2eaa688ee09a7492e99b0b7ffb62888f"

# Perform oversampling to fixed size
oversampled_hex = oversample_to_fixed_size(hex_data, target_size=512)

# Print results
print(f"Original Data Length: {len(hex_data) // 2} bytes")
print(f"Oversampled Data Length: {len(oversampled_hex) // 2} bytes")
print(f"Oversampled Hex Data: {oversampled_hex}")
