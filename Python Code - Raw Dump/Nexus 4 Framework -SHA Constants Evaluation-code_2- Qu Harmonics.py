# Implementation of the full SHA-256 breakdown in pure Python

def sha256_process_block(state, block):
    import struct

    # Constants defined by the SHA-256 algorithm
    K = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
    ]

    # Initialize working variables
    a, b, c, d, e, f, g, h = state

    # Prepare the message schedule
    w = list(struct.unpack('>16L', block)) + [0] * (64 - 16)

    # Extend the first 16 words into the remaining 48 words of the message schedule array
    for i in range(16, 64):
        s0 = (w[i - 15] >> 7 | w[i - 15] << (32 - 7)) ^ \
             (w[i - 15] >> 18 | w[i - 15] << (32 - 18)) ^ \
             (w[i - 15] >> 3)
        s1 = (w[i - 2] >> 17 | w[i - 2] << (32 - 17)) ^ \
             (w[i - 2] >> 19 | w[i - 2] << (32 - 19)) ^ \
             (w[i - 2] >> 10)
        w[i] = (w[i - 16] + s0 + w[i - 7] + s1) & 0xFFFFFFFF

    # Compression function main loop
    for i in range(64):
        s1 = (e >> 6 | e << (32 - 6)) ^ \
             (e >> 11 | e << (32 - 11)) ^ \
             (e >> 25 | e << (32 - 25))
        ch = (e & f) ^ (~e & g)
        temp1 = (h + s1 + ch + K[i] + w[i]) & 0xFFFFFFFF
        s0 = (a >> 2 | a << (32 - 2)) ^ \
             (a >> 13 | a << (32 - 13)) ^ \
             (a >> 22 | a << (32 - 22))
        maj = (a & b) ^ (a & c) ^ (b & c)
        temp2 = (s0 + maj) & 0xFFFFFFFF

        # Update working variables
        h = g
        g = f
        f = e
        e = (d + temp1) & 0xFFFFFFFF
        d = c
        c = b
        b = a
        a = (temp1 + temp2) & 0xFFFFFFFF

    # Update the hash state
    state[0] = (state[0] + a) & 0xFFFFFFFF
    state[1] = (state[1] + b) & 0xFFFFFFFF
    state[2] = (state[2] + c) & 0xFFFFFFFF
    state[3] = (state[3] + d) & 0xFFFFFFFF
    state[4] = (state[4] + e) & 0xFFFFFFFF
    state[5] = (state[5] + f) & 0xFFFFFFFF
    state[6] = (state[6] + g) & 0xFFFFFFFF
    state[7] = (state[7] + h) & 0xFFFFFFFF

    return state

# Example: Processing a single 512-bit block
# Initial hash values
initial_hash = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
]

# 512-bit block (padded message)
block = b"hello" + b'\x80' + b'\x00' * (64 - len(b"hello") - 9) + b'\x00\x00\x00\x00\x00\x00\x00\x28'

# Compute the new state
new_state = sha256_process_block(initial_hash, block)

# Output the new state
new_state
