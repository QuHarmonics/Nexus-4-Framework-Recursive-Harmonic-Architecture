import numpy as np
from math import ceil, log2

###############################################
# Existing wave-based code (truncated for brevity)
###############################################

def compress_wave_input(seed, initial_base=16, final_base=32):
    """
    Similar to your existing compress function, returning a 512-bit container.
    We'll skip re-listing all details. We'll assume it returns:
      final_container_512_bits, metadata_ledger, original_length
    """
    # ... your compression steps ...
    # (For demonstration, we assume the result is exactly 512 bits in final_container)
    final_container_512_bits = [0]*512  # placeholder
    metadata_ledger = []
    original_length = 0
    return final_container_512_bits, metadata_ledger, original_length

def decompress_wave_input(compressed_512_bits, metadata_ledger, original_len):
    """
    Reverses the wave-based compression, assuming we still have 512 bits.
    """
    # ... your decompression steps ...
    # This returns the original text
    return "Hello Turbulent World (original)"

###############################################
# New folding logic
###############################################

def fold_512_to_256(full_512_bits):
    if len(full_512_bits) != 512:
        raise ValueError("Expected exactly 512 bits to fold.")

    folded_256_bits = [0] * 256
    parted_256_bits = [0] * 256

    for i in range(256):
        parted_256_bits[i] = full_512_bits[256 + i]
        folded_256_bits[i] = full_512_bits[i] ^ parted_256_bits[i]

    return folded_256_bits, parted_256_bits

def unfold_256_to_512(folded_256_bits, parted_256_bits):
    if len(folded_256_bits) != 256 or len(parted_256_bits) != 256:
        raise ValueError("folded_256_bits or parted_256_bits must each be 256 bits.")

    full_512_bits = [0] * 512

    for i in range(256):
        full_512_bits[i] = folded_256_bits[i] ^ parted_256_bits[i]
        full_512_bits[256 + i] = parted_256_bits[i]

    return full_512_bits

###############################################
# Revised "final" compress/decompress with folding
###############################################

def final_compress_with_fold(seed):
    """
    1) Wave-based compression -> 512 bits
    2) Fold 512 -> 256 bits
    3) Return final 256-bit array + parted half + ledger
    """
    container_512, ledger, orig_len = compress_wave_input(seed)

    # Step 2: Fold from 512 down to 256
    folded_256, parted_256 = fold_512_to_256(container_512)

    # We store parted_256 so we can unfold later, plus ledger + original length
    # Also store an integer "full_size" for the final container that was folded
    final_flags = {
        'parted_256': parted_256,
        'full_size': 512,
        'ledger': ledger,
        'orig_len': orig_len
    }

    return folded_256, final_flags

def final_decompress_with_unfold(folded_256, final_flags):
    """
    1) Unfold 256 -> 512 bits
    2) Wave-based decompression
    3) Return final text
    """
    parted_256 = final_flags['parted_256']
    full_size = final_flags['full_size']  # e.g., 512
    ledger = final_flags['ledger']
    orig_len = final_flags['orig_len']

    if full_size != 512:
        raise ValueError("This example expects a final container of 512 bits originally.")

    # Step 1: Unfold
    container_512 = unfold_256_to_512(folded_256, parted_256)

    # Step 2: Decompress wave-based
    original_text = decompress_wave_input(container_512, ledger, orig_len)
    return original_text

###############################################
# Example usage
###############################################

if __name__ == "__main__":
    seed = "Hello Turbulent World"

    # Final compression with folding
    folded_hash_256, flags = final_compress_with_fold(seed)

    # In a real scenario, you would store folded_hash_256 (256 bits) + flags somewhere.
    # The parted_256 can be stashed in the flags, or appended, etc.

    # Then to get the original data:
    recovered_text = final_decompress_with_unfold(folded_hash_256, flags)

    print("Original:", seed)
    print("Recovered:", recovered_text)
    print("Match:", seed == recovered_text)
