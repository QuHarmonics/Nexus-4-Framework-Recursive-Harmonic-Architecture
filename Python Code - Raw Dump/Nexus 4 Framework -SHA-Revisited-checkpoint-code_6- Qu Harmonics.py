import math

###############################################################################
# 1. Helper Functions: Bit-Level Operations & Basic Conversions
###############################################################################

def text_to_bits_512_padded(seed, block_size=512):
    """
    Convert 'seed' (string) into a list of bits, then pad to 'block_size' bits.
    If the string requires more than 'block_size' bits, we pad up to the next
    multiple of 'block_size'.

    Returns: (bit_list, original_bit_length, final_block_size)
    """
    # Convert each character to 8-bit ASCII
    bit_str = "".join(f"{ord(c):08b}" for c in seed)
    original_len = len(bit_str)

    if original_len <= block_size:
        # Pad up to exactly block_size
        pad_needed = block_size - original_len
        bit_str += "0" * pad_needed
        final_block_size = block_size
    else:
        # For demonstration, if data exceeds 512 bits, pad up to next multiple
        new_size = math.ceil(original_len / block_size) * block_size
        bit_str = bit_str.ljust(new_size, "0")
        final_block_size = new_size

    bit_list = [int(b) for b in bit_str]
    return bit_list, original_len, final_block_size

def bits_to_text(bit_list, original_len):
    """
    Convert 'bit_list' to text using only the first 'original_len' bits.
    Any partial byte is truncated at the end.
    """
    truncated_bits = bit_list[:original_len]  # use only the original length
    # If there's a partial byte, drop it
    remainder = len(truncated_bits) % 8
    if remainder != 0:
        truncated_bits = truncated_bits[:-remainder]

    text = ""
    for i in range(0, len(truncated_bits), 8):
        chunk = truncated_bits[i:i+8]
        val = int("".join(str(b) for b in chunk), 2)
        text += chr(val)
    return text

def bits_to_int(bit_list):
    """
    Interpret a list of bits as one large integer (base 2).
    """
    return int("".join(str(b) for b in bit_list), 2)

def int_to_bits(value, total_bits=None):
    """
    Convert an integer to a list of bits (base-2).
    If total_bits is given, pad with leading zeros to that bit length.
    """
    bin_str = bin(value)[2:]  # remove '0b'
    if total_bits is not None:
        bin_str = bin_str.zfill(total_bits)
    return [int(x) for x in bin_str]

###############################################################################
# 2. Chunk-Based Base Conversion (One Step: base 2 -> base 17)
###############################################################################

def chunked_binary_to_base17(bit_list, chunk_size=8):
    """
    Interpret 'bit_list' in base 2 in chunks of 'chunk_size' bits,
    convert each chunk to an integer, then represent that integer in base 17.
    We store the result as a list of digits (each in [0..16]).
    
    Returns: (digits_out, chunk_lengths)
      digits_out: concatenated digits in base 17
      chunk_lengths: how many digits each chunk produced (to invert the process)
    """
    # Ensure bit_list length is multiple of chunk_size
    remainder = len(bit_list) % chunk_size
    if remainder != 0:
        bit_list += [0]*(chunk_size - remainder)

    digits_out = []
    chunk_lengths = []

    for i in range(0, len(bit_list), chunk_size):
        chunk = bit_list[i:i+chunk_size]
        val = int("".join(str(b) for b in chunk), 2)  # interpret chunk as base 2 integer

        # Convert that val to base-17 digits
        if val == 0:
            digits_out.append(0)
            chunk_lengths.append(1)
        else:
            tmp = val
            chunk_digs = []
            while tmp > 0:
                chunk_digs.append(tmp % 17)
                tmp //= 17
            chunk_digs.reverse()
            digits_out.extend(chunk_digs)
            chunk_lengths.append(len(chunk_digs))

    return digits_out, chunk_lengths

def chunked_base17_to_binary(digits17, chunk_lengths, chunk_size=8):
    """
    Reverse of chunked_binary_to_base17.
    We read digits17 according to chunk_lengths. For each chunk:
      1) Convert those base-17 digits to an integer val.
      2) Convert val -> chunk_size bits in base-2.
    """
    bits_out = []
    idx = 0
    for length in chunk_lengths:
        sub = digits17[idx: idx+length]
        idx += length

        val = 0
        for d in sub:
            val = val*17 + d

        # Convert val to chunk_size bits
        bin_str = bin(val)[2:].zfill(chunk_size)
        bits_out.extend(int(x) for x in bin_str)
    return bits_out

###############################################################################
# 3. Single-Step Wave Compression to 512 bits (Base 2 -> Base 17 -> Base 2)
###############################################################################

def wave_compress_512(seed):
    """
    1) Convert seed -> bits, pad up to 512 (or multiple).
    2) Chunk-based convert (base 2 -> base 17).
    3) Convert those digits back to bits (base 17 -> base 2).
    4) Pad/truncate to 512 bits.
    5) Return final bits, plus metadata for reversing.
    """
    bit_list, original_len, final_block_size = text_to_bits_512_padded(seed, block_size=512)

    # Step A: chunk-based 2 -> 17
    digits17, chunk_lengths = chunked_binary_to_base17(bit_list, chunk_size=8)

    # Step B: back from base 17 -> base 2
    recovered_bits = chunked_base17_to_binary(digits17, chunk_lengths, chunk_size=8)

    # If recovered_bits is shorter/longer than final_block_size, pad or truncate
    if len(recovered_bits) < final_block_size:
        recovered_bits += [0]*(final_block_size - len(recovered_bits))
    elif len(recovered_bits) > final_block_size:
        recovered_bits = recovered_bits[:final_block_size]

    # Build metadata ledger
    ledger = {
        'chunk_size': 8,
        'final_block_size': final_block_size,
        'chunk_lengths': chunk_lengths,
    }

    return recovered_bits, ledger, original_len

def wave_decompress_512(bits_512, ledger, original_bit_len):
    """
    Reverse wave_compress_512:
      1) We interpret bits_512 as base-2 -> base-17 digits (the same chunk pattern).
      2) Then re-chunk from base 17 -> base 2, hopefully recovering the original bits
      3) Convert to text with 'original_bit_len'.
    """
    final_block_size = ledger['final_block_size']
    chunk_lengths = ledger['chunk_lengths']
    chunk_size = ledger['chunk_size']

    # Step A: chunk-based 2 -> 17
    # Actually, we are going from bits_512 => digits17, using the same logic:
    # We'll do exactly chunked_binary_to_base17 to see the digits
    digits17, lengths_ignored = chunked_binary_to_base17(bits_512, chunk_size=chunk_size)
    # The chunk_lengths might differ from the original if there's any mismatch,
    # but let's trust we didn't lose data.

    # Step B: base 17 -> base 2 using the stored chunk_lengths
    recovered_bits = chunked_base17_to_binary(digits17, chunk_lengths, chunk_size=chunk_size)

    # Truncate/pad to final_block_size
    if len(recovered_bits) < final_block_size:
        recovered_bits += [0]*(final_block_size - len(recovered_bits))
    elif len(recovered_bits) > final_block_size:
        recovered_bits = recovered_bits[:final_block_size]

    # Step C: convert to text using original_bit_len
    return bits_to_text(recovered_bits, original_bit_len)

###############################################################################
# 4. Folding & Unfolding
###############################################################################

def fold_512_to_256(full_512_bits):
    """
    XOR-based fold: 
      - parted_256[i] = full_512_bits[256 + i]
      - folded_256[i] = full_512_bits[i] XOR parted_256[i]
    """
    if len(full_512_bits) != 512:
        raise ValueError("Expected exactly 512 bits to fold.")
    folded_256_bits = [0]*256
    parted_256_bits = [0]*256
    for i in range(256):
        parted_256_bits[i] = full_512_bits[256 + i]
        folded_256_bits[i] = full_512_bits[i] ^ parted_256_bits[i]
    return folded_256_bits, parted_256_bits

def unfold_256_to_512(folded_256_bits, parted_256_bits):
    if len(folded_256_bits) != 256 or len(parted_256_bits) != 256:
        raise ValueError("folded_256_bits or parted_256_bits must each be 256 bits.")
    full_512_bits = [0]*512
    for i in range(256):
        full_512_bits[i] = folded_256_bits[i] ^ parted_256_bits[i]
        full_512_bits[256 + i] = parted_256_bits[i]
    return full_512_bits

###############################################################################
# 5. Final "Compress with Fold" & "Decompress"
###############################################################################

def final_compress_with_fold(seed):
    """
    1) wave_compress_512 -> 512 bits
    2) fold_512_to_256
    3) Store parted_256 and ledger + original length
    """
    final_512_bits, ledger, orig_len = wave_compress_512(seed)

    folded_256, parted_256 = fold_512_to_256(final_512_bits)
    final_flags = {
        'parted_256': parted_256,
        'full_size': 512,
        'ledger': ledger,
        'orig_len': orig_len
    }
    return folded_256, final_flags

def final_decompress_with_unfold(folded_256, final_flags):
    parted_256 = final_flags['parted_256']
    full_size = final_flags['full_size']
    ledger = final_flags['ledger']
    orig_len = final_flags['orig_len']

    if full_size != 512:
        raise ValueError("This example only handles 512-bit containers.")
    
    # 1) Unfold
    container_512 = unfold_256_to_512(folded_256, parted_256)
    # 2) wave_decompress_512
    recovered_text = wave_decompress_512(container_512, ledger, orig_len)
    return recovered_text

###############################################################################
# 6. Demo: Feed Output Into Next Input
###############################################################################

def demo_double_compression(seed):
    """
    1) Compress 'seed' into folded_256 + parted_256.
    2) Convert folded_256 bits to hex, feed that hex as new seed in Round 2.
    3) Observe results.
    """
    print("=== ROUND 1 ===")
    folded_256, flags = final_compress_with_fold(seed)
    
    # Convert folded_256 to an integer for demonstration
    folded_int = 0
    for bit in folded_256:
        folded_int = (folded_int << 1) | bit
    # produce hex
    folded_hex = hex(folded_int)[2:].upper()
    print(f"Folded 256 bits as hex (truncated): {folded_hex[:64]}...")

    # Decompress to ensure we can get the original
    recovered_text = final_decompress_with_unfold(folded_256, flags)
    print("Recovered Round1:", recovered_text)
    print("Match Round1:", recovered_text == seed)

    # Round 2: feed folded_hex as new seed
    print("\n=== ROUND 2 ===")
    folded_256b, flags2 = final_compress_with_fold(folded_hex)
    recovered_text2 = final_decompress_with_unfold(folded_256b, flags2)
    print("Recovered Round2:", recovered_text2)
    print("Note: Round2 won't match the original seed, because we compressed the hex.")
    return folded_256b, flags2

###############################################################################
# Main
###############################################################################

if __name__ == "__main__":
    seed = "Hello Turbulent World"
    print("Original Seed:", seed)

    # 1) Single pass compress + fold
    folded_out, flags_out = final_compress_with_fold(seed)
    recovered = final_decompress_with_unfold(folded_out, flags_out)
    print("\nFinal Recovered:", recovered)
    print("Match:", recovered == seed)

    # 2) Double compression demo
    print("\nNow demonstrating double compression with feed-back loop:\n")
    demo_double_compression(seed)
