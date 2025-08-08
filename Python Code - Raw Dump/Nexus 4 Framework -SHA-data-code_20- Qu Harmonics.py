import hashlib
import struct

# The 13‐cycle you found
cycle = [21, 8, 154, 65, 253, 155, 94, 203, 29, 205, 209, 153, 88]

def build_schedule(byte_val):
    # Step A: single byte message
    m = bytes([byte_val])
    # Step B: padding (0x80, zeros, then length 8-bit)
    ml = len(m) * 8
    pad = b'\x80' + b'\x00' * ((56 - (len(m)+1)) % 64) + struct.pack('>Q', ml)
    block = m + pad  # exactly 64 bytes

    # Parse into 16 big-endian 32-bit words
    W = list(struct.unpack('>16I', block))
    # Expand to 64 words
    def rotr(x,n): return ((x >> n) | (x << (32-n))) & 0xFFFFFFFF
    def sigma0(x): return (rotr(x,7) ^ rotr(x,18) ^ (x >> 3)) & 0xFFFFFFFF
    def sigma1(x): return (rotr(x,17) ^ rotr(x,19) ^ (x >> 10)) & 0xFFFFFFFF

    for t in range(16,64):
        val = (sigma1(W[t-2]) + W[t-7] + sigma0(W[t-15]) + W[t-16]) & 0xFFFFFFFF
        W.append(val)
    return W

# Collect which W_t mod256 hits each cycle-state
matches = {s: [] for s in cycle}

for s in cycle:
    W = build_schedule(s)
    W_mod = [w & 0xFF for w in W]
    for t, byte in enumerate(W_mod):
        if byte == s:
            matches[s].append(t)

print("State → matching schedule‐indices:", matches)
