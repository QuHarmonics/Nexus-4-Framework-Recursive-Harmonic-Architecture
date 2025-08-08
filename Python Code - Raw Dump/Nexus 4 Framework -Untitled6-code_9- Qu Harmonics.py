from mpmath import mp

# Set precision for high-position digits
mp.dps = 5000

def compute_bbp_digit(n):
    """Compute the nth hexadecimal digit of π using the BBP formula."""
    s = mp.mpf(0)
    # Sum up to n-1
    for k in range(n):
        exponent = n - 1 - k
        if exponent >= 0:
            term = (mp.power(16, exponent) % (8*k + 1)) / (8*k + 1) - \
                   (mp.power(16, exponent) % (8*k + 4)) / (8*k + 4) - \
                   (mp.power(16, exponent) % (8*k + 5)) / (8*k + 5) - \
                   (mp.power(16, exponent) % (8*k + 6)) / (8*k + 6)
        else:
            term = mp.power(16, exponent) * (4 / (8*k + 1) - 2 / (8*k + 4) - 
                                             1 / (8*k + 5) - 1 / (8*k + 6))
        s += term
    # Approximate tail sum
    for k in range(n, n + 50):
        term = mp.power(16, n - 1 - k) * (4 / (8*k + 1) - 2 / (8*k + 4) - 
                                          1 / (8*k + 5) - 1 / (8*k + 6))
        s += term
    # Extract digit
    s = s % 1
    digit = int(16 * s)
    return hex(digit)[2:].upper()

# Example: Compute digits for positions 1 to 10
positions = range(1, 11)
digits = [compute_bbp_digit(n) for n in positions]
print("Position to Digit Mapping:")
for n, digit in zip(positions, digits):
    print(f"Position {n}: {digit}")