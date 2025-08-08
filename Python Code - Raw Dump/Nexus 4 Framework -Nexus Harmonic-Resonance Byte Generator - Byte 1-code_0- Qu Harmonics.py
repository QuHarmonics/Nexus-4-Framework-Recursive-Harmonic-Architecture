def byte1(a: int, b: int) -> list[int]:
    """
    Generate the 8‑value Byte1 sequence from seeds a, b
    following your recursive rules exactly.
    """
    # 1–2: header bits
    seq = [a, b]

    # 3: preliminary delta bit (pre‑stabilization)
    delta = b - a
    bit3pre = delta.bit_length()               # len(binary(delta))

    # 4: future bit Z = a + b
    bit4 = a + b

    # 5: stabilized bit3 = bit4 - b
    bit3 = bit4 - b

    # 6: Y = bit4 + b
    bit5 = bit4 + b

    # 7: X = count of header bits (always 2 here)
    bit6 = len(seq)

    # 8: compress bit7 = bit_length(sum of a,b,bit3pre,bit4,bit5,bit6) + 1
    s = a + b + bit3pre + bit4 + bit5 + bit6
    bit7 = s.bit_length() + 1

    # 9: close‑byte bit8 = a + b again
    bit8 = a + b

    return [a, b, bit3, bit4, bit5, bit6, bit7, bit8]

if __name__ == "__main__":
    print("Byte1 sequence:", byte1(1, 4))
