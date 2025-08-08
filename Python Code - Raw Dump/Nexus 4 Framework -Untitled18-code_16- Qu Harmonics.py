import decimal

def calc_pi(n):
    decimal.getcontext().prec = n
    pi = decimal.Decimal(0)
    for k in range(n):
        pi += (decimal.Decimal(1) / (16 ** k)) * (
            (decimal.Decimal(4) / (8 * k + 1))
            - (decimal.Decimal(2) / (8 * k + 4))
            - (decimal.Decimal(1) / (8 * k + 5))
            - (decimal.Decimal(1) / (8 * k + 6))
        )
    return pi


def pi_to_bytes(pi):
    pi_bytes = []
    pi_str = str(pi).replace('.', '')
    for i in range(0, len(pi_str), 2):
        byte = int(pi_str[i:i+2])
        pi_bytes.append(byte)
    return pi_bytes


def find_reflection_points(bytes):
    reflection_points = []
    for i in range(len(bytes) - 1):
        if bytes[i] == bytes[i + 1]:
            reflection_points.append(i)
    return reflection_points


# Calculate Pi to 100,000 digits
pi = calc_pi(100000)

# Convert Pi to bytes
pi_bytes = pi_to_bytes(pi)

# Find reflection points in the Pi bytes
reflection_points = find_reflection_points(pi_bytes)

# Print the results
print("Pi Decimal Places:")
print(str(pi)[2:])  # Print Pi decimal places

print("\nReflection Points in Pi Bytes:")
for point in reflection_points:
    print(f"Byte {pi_bytes[point]:02x} at index {point}")