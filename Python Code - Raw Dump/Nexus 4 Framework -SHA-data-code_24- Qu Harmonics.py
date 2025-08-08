import struct

# your 13-cycle
cycle = [21, 8, 154, 65, 253, 155, 94, 203, 29, 205, 209, 153, 88]

def build_schedule_bytes(byte_val):
    # same padding + schedule as before
    m = bytes([byte_val])
    ml = len(m)*8
    pad = b'\x80' + b'\x00'*((56 - (len(m)+1)) % 64) + struct.pack('>Q', ml)
    block = m + pad

    W = list(struct.unpack('>16I', block))
    def rotr(x,n): return ((x>>n)|(x<<(32-n))) & 0xFFFFFFFF
    def σ0(x): return (rotr(x,7)^rotr(x,18)^(x>>3)) & 0xFFFFFFFF
    def σ1(x): return (rotr(x,17)^rotr(x,19)^(x>>10)) & 0xFFFFFFFF
    for t in range(16,64):
        W.append((σ1(W[t-2]) + W[t-7] + σ0(W[t-15]) + W[t-16]) & 0xFFFFFFFF)

    # turn each W[t] into its 4 bytes (big-endian)
    return [struct.pack('>I', w) for w in W]

matches = {s: [] for s in cycle}
for s in cycle:
    sched_bytes = build_schedule_bytes(s)
    for t, word_bytes in enumerate(sched_bytes):
        # check each of the 4 bytes
        for i, b in enumerate(word_bytes):
            if b == s:
                matches[s].append((t, i))   # (round-index, byte-position 0–3)

print("State → [(round, byte-pos), …]:")
for s, hits in matches.items():
    print(f"  {s:3d} → {hits}")
