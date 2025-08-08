#!/usr/bin/env python3
"""
bitwise_then_bruteforce.py —  
1) Attempt bit-by-bit nonce derivation;  
2) If it “saturates” at 0xFFFFFFFF, abandon and do a full 0…2^32-1 scan.
"""

import hashlib
import struct
import sys

def double_sha256(data: bytes) -> bytes:
    return hashlib.sha256(hashlib.sha256(data).digest()).digest()

def bits_to_target(bits: int) -> int:
    exp = bits >> 24
    coeff = bits & 0x007FFFFF
    return coeff << (8 * (exp - 3))

def find_nonce_bitwise(prefix: bytes, target: int) -> int:
    nonce = 0
    for bit in reversed(range(32)):
        trial = nonce | (1 << bit)
        header = prefix + struct.pack('<I', trial)
        h = int.from_bytes(double_sha256(header), 'little')
        if h <= target:
            nonce = trial
    # If we end up with the all-ones value, treat it as failure
    if nonce == 0xFFFFFFFF:
        return None
    return nonce

def find_nonce_bruteforce(prefix: bytes, target: int) -> int:
    for nonce in range(2**32):
        header = prefix + struct.pack('<I', nonce)
        if int.from_bytes(double_sha256(header), 'little') <= target:
            return nonce
    return None

def main():
    # --- Example header fields (replace with real values) ---
    version     = struct.pack('<I', 0x20000000)
    prev_hash   = bytes.fromhex('00'*32)
    merkle_root = bytes.fromhex('4e3b2a1d0c9f8e7d6c5b4a3928171615041322110f0e1d2c3b4a5d6e7f8a9b0c')
    timestamp   = struct.pack('<I', 1622505600)
    bits        = struct.pack('<I', 0x170d3f6f)

    prefix   = version + prev_hash + merkle_root + timestamp + bits
    bits_int = struct.unpack('>I', bits)[0]
    target   = bits_to_target(bits_int)

    print(f"Decoded target: 0x{target:064x}\n")

    # 1) Try greedy bitwise
    print("Attempting bit-wise derivation…")
    nonce = find_nonce_bitwise(prefix, target)
    if nonce is not None:
        print(f"✔ Bit-wise success: nonce = 0x{nonce:08x}")
    else:
        print("⚠ Bit-wise approach failed (returned 0xFFFFFFFF).")
        print("Falling back to full brute-force scan…")
        # 2) Full-space scan
        nonce = find_nonce_bruteforce(prefix, target)
        if nonce is None:
            print("✘ No valid nonce found in 0…2^32−1.", file=sys.stderr)
            sys.exit(1)
        print(f"✔ Brute-force success: nonce = 0x{nonce:08x}")

    # Final verification
    final_header = prefix + struct.pack('<I', nonce)
    final_hash   = double_sha256(final_header)[::-1].hex()
    print(f"Resulting hash (big-endian): {final_hash}")

if __name__ == '__main__':
    main()
