#!/usr/bin/env python3
"""
Bitwise Nonce Finder for Bitcoin Proof-of-Work (no ASIC)

This script builds the 32-bit nonce one bit at a time (MSB→LSB),
testing each candidate against the full double-SHA256 target check.
Total hash calls ≈ 2×32 = 64, versus brute-force 2**32 ≈ 4e9.

Requirements:
    - Python 3.6+
    - hashlib (built-in)
"""

import hashlib
import struct

def double_sha256(data: bytes) -> bytes:
    """Compute SHA256(SHA256(data))."""
    return hashlib.sha256(hashlib.sha256(data).digest()).digest()

def bits_to_target(bits: int) -> int:
    """
    Decode the compact 'bits' field to a full target.
    bits: 4-byte compact representation (big-endian uint32)
    """
    exponent = bits >> 24
    coeff    = bits & 0x007fffff
    # If the sign bit (0x00800000) is set, coeff is negative in original spec,
    # but Bitcoin never uses negative coefficients.
    return coeff * (1 << (8 * (exponent - 3)))

def find_nonce_bitwise(header_prefix: bytes, target: int) -> int:
    """
    Find the 32-bit nonce by setting one bit at a time.
    header_prefix: first 76 bytes of the header (everything but the nonce)
    target: integer target threshold
    Returns the integer nonce that yields hash ≤ target.
    """
    nonce = 0
    # Work from most-significant bit (31) down to least (0)
    for bit in reversed(range(32)):
        # Try setting this bit to 1
        trial = nonce | (1 << bit)
        # Pack trial nonce in little-endian and hash
        header = header_prefix + struct.pack('<I', trial)
        h = int.from_bytes(double_sha256(header), 'little')
        if h <= target:
            # Keep the bit set
            nonce = trial
        # Else leave it cleared (bit remains 0)
    # Final verification
    final_header = header_prefix + struct.pack('<I', nonce)
    final_hash   = int.from_bytes(double_sha256(final_header), 'little')
    if final_hash > target:
        raise ValueError(f"No valid nonce found; last hash {final_hash:064x} > target")
    return nonce

def main():
    # --- Configure your block header fields here ---
    version      = struct.pack('<I', 0x20000000)  # version 536870912 (example)
    prev_hash    = bytes.fromhex(
        '0000000000000000000b4d0f5f2e8e2c4b1a6f3d2c1e0f9a8b7c6d5e4f3a2b1c'
    )  # little-endian previous block hash
    merkle_root  = bytes.fromhex(
        '4e3b2a1d0c9f8e7d6c5b4a3928171615041322110f0e1d2c3b4a5d6e7f8a9b0c'
    )  # example merkle root (little-endian)
    timestamp    = struct.pack('<I', 1622505600)  # 2021-06-01 00:00:00 UTC
    bits         = struct.pack('<I', 0x170d3f6f)  # example compact target
    
    # --- Build prefix and compute target ---
    header_prefix = version + prev_hash + merkle_root + timestamp + bits
    bits_int      = struct.unpack('>I', bits)[0]   # interpret bits as big-endian uint32
    target        = bits_to_target(bits_int)
    
    print(f"Target (hex): {target:064x}")
    print("Searching nonce bit by bit...")
    
    nonce = find_nonce_bitwise(header_prefix, target)
    
    # Assemble final header and compute final hash
    final_header = header_prefix + struct.pack('<I', nonce)
    final_hash   = double_sha256(final_header)[::-1].hex()  # display big-endian
    
    print(f"✔ Found nonce: {nonce} (0x{nonce:08x})")
    print(f"✔ Resulting hash: {final_hash}")

if __name__ == "__main__":
    main()
