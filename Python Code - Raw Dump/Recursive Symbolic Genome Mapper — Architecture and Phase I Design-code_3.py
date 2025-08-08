from decimal import Decimal, getcontext

# Set precision
getcontext().prec = 50

# Byte system
bytes_full = [
    [1, 4, 1, 5, 9, 2, 6, 5],
    [3, 5, 8, 9, 7, 9, 3, 2],
    [3, 8, 4, 6, 2, 6, 4, 3],
    [3, 8, 3, 2, 7, 9, 5, 0],
    [2, 8, 8, 4, 1, 9, 7, 1],
    [6, 9, 3, 9, 9, 3, 7, 5],
    [1, 0, 5, 8, 2, 0, 9, 7],
    [4, 5, 9, 2, 3, 0, 7, 8]
]

# Flatten to single control stream
control = [n for row in bytes_full for n in row]

# Generator function
def byte_pi_nilakantha(iterations=100):
    pi = Decimal(3)
    sign = 1
    for i in range(1, iterations):
        a = Decimal(control[i % len(control)])
        b = Decimal(control[(i+1) % len(control)])
        c = Decimal(control[(i+2) % len(control)])

        # Skip zeros to avoid div-by-zero
        if a == 0 or b == 0 or c == 0:
            continue

        term = Decimal(4) / (a * b * c)
        pi += sign * term
        sign *= -1  # alternate sign
    return +pi

# Run it
approx_pi = byte_pi_nilakantha(300)
print("Byte-synthesized π approximation:", approx_pi)
