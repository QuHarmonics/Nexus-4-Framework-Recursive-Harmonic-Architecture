import hashlib
import struct

# SHA-256 Logical Functions
def rotr(x, n, bits=32):
    """Bitwise rotate right."""
    return ((x >> n) | (x << (bits - n))) & (2**bits - 1)

def sha256_logic_functions(a, b, c, d, e, f, g, h, k, w):
    """Simulates one SHA-256 round."""
    # SHA-256 logical functions
    ch = (e & f) ^ (~e & g)  # Choose
    maj = (a & b) ^ (a & c) ^ (b & c)  # Majority
    s1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)  # Sigma1
    s0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)  # Sigma0

    temp1 = (h + s1 + ch + k + w) & 0xFFFFFFFF
    temp2 = (s0 + maj) & 0xFFFFFFFF
    h = g
    g = f
    f = e
    e = (d + temp1) & 0xFFFFFFFF
    d = c
    c = b
    b = a
    a = (temp1 + temp2) & 0xFFFFFFFF

    return a, b, c, d, e, f, g, h

def reverse_sha256_round(a, b, c, d, e, f, g, h, k, w):
    """Attempts to reverse a SHA-256 round."""
    # Reverse the transformations step-by-step
    temp1 = (a - ((rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)) + ((e & f) ^ (~e & g)) + k + w)) & 0xFFFFFFFF
    temp2 = (temp1 - ((rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)) + ((a & b) ^ (a & c) ^ (b & c)))) & 0xFFFFFFFF
    h = temp1
    g = h
    f = g
    e = d - temp1
    d = c
    c = b
    b = a
    a = temp2

    return a, b, c, d, e, f, g, h

# Test the Round
def test_sha256_round():
    # Test data (dummy example)
    a, b, c, d, e, f, g, h = (0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                              0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19)
    k = 0x428a2f98  # Example constant
    w = 0x61626364  # First word ("abcd" in hex)

    # Forward round
    new_a, new_b, new_c, new_d, new_e, new_f, new_g, new_h = sha256_logic_functions(a, b, c, d, e, f, g, h, k, w)
    print("Forward SHA-256 Round:")
    print(f"Output: {new_a:08x} {new_b:08x} {new_c:08x} {new_d:08x} {new_e:08x} {new_f:08x} {new_g:08x} {new_h:08x}")

    # Reverse round
    rev_a, rev_b, rev_c, rev_d, rev_e, rev_f, rev_g, rev_h = reverse_sha256_round(new_a, new_b, new_c, new_d, new_e, new_f, new_g, new_h, k, w)
    print("Reverse SHA-256 Round:")
    print(f"Recovered: {rev_a:08x} {rev_b:08x} {rev_c:08x} {rev_d:08x} {rev_e:08x} {rev_f:08x} {rev_g:08x} {rev_h:08x}")

test_sha256_round()
