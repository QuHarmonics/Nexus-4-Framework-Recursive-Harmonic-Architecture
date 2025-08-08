def sha256(data: bytes) -> str:
    # 1) pad 'data' to a multiple of 64 bytes (omitted for brevity)
    state = [
      0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
    ]
    for chunk in chunks(data, 64):
        W = list(unpack_big_endian_words(chunk))
        for t in range(16, 64):
            s0 = (ROTR(W[t-15],7) ^ ROTR(W[t-15],18) ^ (W[t-15]>>3))
            s1 = (ROTR(W[t-2],17) ^ ROTR(W[t-2],19) ^ (W[t-2]>>10))
            W.append((W[t-16] + s0 + W[t-7] + s1) & 0xFFFFFFFF)

        a,b,c,d,e,f,g,h = state
        for t in range(64):
            T1 = (h + (ROTR(e,6)^ROTR(e,11)^ROTR(e,25))
                  + ((e&f) ^ (~e&g)) + K[t] + W[t]) & 0xFFFFFFFF
            T2 = ((ROTR(a,2)^ROTR(a,13)^ROTR(a,22))
                  + ((a&b) ^ (a&c) ^ (b&c))) & 0xFFFFFFFF
            h, g, f, e, d, c, b, a = g, f, e, (d+T1)&0xFFFFFFFF, c, b, a, (T1+T2)&0xFFFFFFFF

        # Add back into state
        state = [(s + v) & 0xFFFFFFFF for s,v in zip(state, [a,b,c,d,e,f,g,h])]

    return ''.join(f'{x:08x}' for x in state)
