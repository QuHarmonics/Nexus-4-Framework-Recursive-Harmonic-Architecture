# Precomputed SHA-256 round constants (first 32 bits)
K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, # … up to 64 total
]

# Your 13-cycle states:
cycle = [21, 8, 154, 65, 253, 155, 94, 203, 29, 205, 209, 153, 88]

# Map K mod 256
K_mod = [k & 0xff for k in K]

# See which cycle members appear in K_mod
hits = {s: [i for i, byte in enumerate(K_mod) if byte == s] for s in cycle}
print("State → matching round-constant indices:", hits)
