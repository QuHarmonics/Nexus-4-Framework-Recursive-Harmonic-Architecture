#!/usr/bin/env python3
"""
full_sha256.py

A correct, from-scratch implementation of the SHA-256 hash algorithm in pure Python.

Usage:
    python full_sha256.py

Produces the SHA-256 digest of "hello world" as a demonstration.
"""

from typing import Iterator, List

# SHA-256 round constants
K: List[int] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b,
    0x59f111f1, 0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01,
    0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7,
    0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152,
    0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147,
    0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
    0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819,
    0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116, 0x1e376c08,
    0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f,
    0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]

# Initial hash values (first 32 bits of the fractional parts of the square roots of the first 8 primes)
H0: List[int] = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

def right_rotate(x: int, n: int) -> int:
    return ((x >> n) | (x << (32 - n))) & 0xFFFFFFFF

def pad_message(message: bytes) -> bytes:
    """Pad the message to a multiple of 512 bits."""
    ml = len(message) * 8
    message += b'\x80'
    # pad with zeros until length ≡ 448 mod 512 (i.e. 56 bytes mod 64)
    message += b'\x00' * ((56 - (len(message) % 64)) % 64)
    # append 64-bit big-endian length
    message += ml.to_bytes(8, 'big')
    return message

def chunk_blocks(padded: bytes) -> Iterator[bytes]:
    """Yield successive 64-byte blocks."""
    for i in range(0, len(padded), 64):
        yield padded[i:i+64]

def sha256(message: bytes) -> str:
    # 1) Preprocessing
    padded = pad_message(message)

    # 2) Initialize working variables
    H = H0.copy()

    # 3) Process each 512-bit block
    for block in chunk_blocks(padded):
        # prepare message schedule W[0..63]
        W: List[int] = [0] * 64
        # first 16 words from block
        for t in range(16):
            W[t] = int.from_bytes(block[t*4:(t+1)*4], 'big')
        # remaining words
        for t in range(16, 64):
            s0 = (right_rotate(W[t-15], 7) ^
                  right_rotate(W[t-15], 18) ^
                  (W[t-15] >> 3))
            s1 = (right_rotate(W[t-2], 17) ^
                  right_rotate(W[t-2], 19) ^
                  (W[t-2] >> 10))
            W[t] = (W[t-16] + s0 + W[t-7] + s1) & 0xFFFFFFFF

        # initialize a..h
        a, b, c, d, e, f, g, h = H

        # main compression loop
        for t in range(64):
            S1 = (right_rotate(e, 6) ^
                  right_rotate(e, 11) ^
                  right_rotate(e, 25))
            ch = (e & f) ^ ((~e) & g)
            temp1 = (h + S1 + ch + K[t] + W[t]) & 0xFFFFFFFF
            S0 = (right_rotate(a, 2) ^
                  right_rotate(a, 13) ^
                  right_rotate(a, 22))
            maj = (a & b) ^ (a & c) ^ (b & c)
            temp2 = (S0 + maj) & 0xFFFFFFFF

            h = g
            g = f
            f = e
            e = (d + temp1) & 0xFFFFFFFF
            d = c
            c = b
            b = a
            a = (temp1 + temp2) & 0xFFFFFFFF

        # add to intermediate hash value
        H = [
            (H[0] + a) & 0xFFFFFFFF,
            (H[1] + b) & 0xFFFFFFFF,
            (H[2] + c) & 0xFFFFFFFF,
            (H[3] + d) & 0xFFFFFFFF,
            (H[4] + e) & 0xFFFFFFFF,
            (H[5] + f) & 0xFFFFFFFF,
            (H[6] + g) & 0xFFFFFFFF,
            (H[7] + h) & 0xFFFFFFFF,
        ]

    # produce final hash (big-endian)
    return ''.join(f'{value:08x}' for value in H)

if __name__ == "__main__":
    input_data = b"hello world"
    print(f"SHA-256 Hash of '{input_data.decode()}': {sha256(input_data)}")
