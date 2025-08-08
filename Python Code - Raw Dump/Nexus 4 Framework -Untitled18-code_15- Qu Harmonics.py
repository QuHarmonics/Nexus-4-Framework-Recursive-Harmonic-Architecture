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


def find_double_headers(digits):
    double_headers = []
    for i in range(len(digits) - 1):
        if digits[i] == digits[i + 1]:
            double_headers.append([digits[i], i])
    return double_headers


def pi_to_bytes(pi):
    pi_bytes = []
    pi_str = str(pi).replace('.', '')
    for i in range(0, len(pi_str), 2):
        byte = int(pi_str[i:i+2])
        pi_bytes.append(byte)
    return pi_bytes


# Calculate Pi to 1,000,000 digits
pi = calc_pi(10000)

# Convert Pi to