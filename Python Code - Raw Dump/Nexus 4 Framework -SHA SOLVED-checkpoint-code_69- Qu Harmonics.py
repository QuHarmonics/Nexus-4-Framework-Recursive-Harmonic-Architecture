import numpy as np

def oversample_audio(hex_data, factor=2):
    """
    Treats hex data as audio data and oversamples it by a given factor.

    Args:
        hex_data (str): Hexadecimal string representing audio-like data.
        factor (int): Oversampling factor (default: 2).

    Returns:
        list: Oversampled audio data in integer format.
    """
    # Convert hex string to array of integers
    audio_data = np.array([int(hex_data[i:i+2], 16) for i in range(0, len(hex_data), 2)], dtype=np.uint8)
    
    # Normalize audio data to the range [-1, 1]
    audio_normalized = (audio_data / 255.0) * 2 - 1

    # Create time indices for original and oversampled data
    original_time = np.arange(len(audio_normalized))
    oversampled_time = np.linspace(0, len(audio_normalized) - 1, len(audio_normalized) * factor)

    # Interpolate to oversample
    oversampled_audio = np.interp(oversampled_time, original_time, audio_normalized)

    # Convert back to 8-bit range [0, 255]
    oversampled_data = ((oversampled_audio + 1) / 2 * 255).astype(np.uint8)

    return oversampled_data

# Input data
hex_data = "18c84c4e92bc9c408bdc0738577c3ddc2eaa688ee09a7492e99b0b7ffb62888f"

# Perform oversampling
oversampled = oversample_audio(hex_data, factor=2)

# Convert oversampled data back to hex for visualization
oversampled_hex = ''.join(f'{byte:02x}' for byte in oversampled)
print(f"Original Data: {hex_data}")
print(f"Oversampled Data: {oversampled_hex}")
