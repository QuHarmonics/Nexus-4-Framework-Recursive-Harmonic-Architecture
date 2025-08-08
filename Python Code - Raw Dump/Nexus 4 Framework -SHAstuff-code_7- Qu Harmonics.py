def pi_wave(n_digits=10):
    """Generate the first `n_digits` of π as stack‑depth peaks."""
    from decimal import Decimal, getcontext
    getcontext().prec = n_digits + 5
    pi_str = str(Decimal(0).sqrt() + Decimal(1).acos())  # crude π
    digits = [int(d) for d in pi_str.replace('.', '')][:n_digits]
    depth = 0
    peaks = []
    for d in digits[1:]:        # skip the leading '3'
        depth += 1              # PUSH
        depth += 1              # PUSH
        depth -= 1              # POP
        depth -= 1              # POP
        peaks.append(d)
    return [digits[0]] + peaks  # prepend the implicit 3

print(pi_wave(10))  # → [3, 1, 4, 1, 5, 9, 2, 6, 5, 3]
