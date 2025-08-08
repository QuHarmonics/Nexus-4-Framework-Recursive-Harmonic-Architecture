import math

def pi_hex_digits(start, count):
    """Return `count` hex digits of π beginning at 1-indexed position `start`."""
    def S(j, n):
        # first sum: k=0..n  (modular)
        total = 0.0
        for k in range(n+1):
            total += pow(16, n-k, 8*k + j) / (8*k + j)
        # tail: k=n+1..∞ until term tiny
        term = 0.0
        k = n + 1
        while True:
            t = 16**(n-k) / (8*k + j)
            if t < 1e-17:
                break
            term += t
            k += 1
        return (total + term) % 1

    hex_str = ""
    for i in range(count):
        n = start + i - 1
        x = (4*S(1,n) - 2*S(4,n) -  S(5,n) -  S(6,n)) % 1
        hex_str += "{:X}".format(int(x * 16))
    return hex_str

# Your positions:
positions = [35714677, 40718987, 197596491]
for pos in positions:
    print(f"π hex @ {pos:,} →", pi_hex_digits(pos, 8))
