from mpmath import mp

def fold_pi_digits(chunks, mode="position_product", factor=2):
    resonances = []
    for chunk_idx, chunk in enumerate(chunks):
        if mode == "mirror_sum":
            left, right = chunk[:4], chunk[4:]
            sums = [left[i] + right[i] for i in range(4)]
            resonance = sum(sums) / 36  # sum(sums) / (9*4)
        elif mode == "position_product":
            positions = list(range(2, 10))  # Always positions 2-9 for all chunks, like your original
            products = [d * p for d, p in zip(chunk, positions)]
            resonance = sum(products) / 63  # sum(products) / (7*9), your original scaling (7*9=63)
        elif mode == "frequency":
            two_positions = [i+1 for i, d in enumerate(chunk) if d == 2]
            gaps = [two_positions[i+1] - two_positions[i] for i in range(len(two_positions)-1)] if len(two_positions) > 1 else [0]
            resonance = sum(gaps) / 16  # (8*2)
        else:
            digit_str = ''.join(map(str, chunk))
            decimal = int(digit_str) / 100000000
            resonance = decimal
            for _ in range(3):
                resonance = (resonance * factor) % 1
        resonances.append(resonance)
    return resonances

# Get Pi digits
mp.dps = 100
pi_str = str(mp.pi())[2:]
pi_digits = [int(d) for d in pi_str[:80]]
chunks = [pi_digits[i:i+8] for i in range(0, len(pi_digits), 8) if len(pi_digits[i:i+8]) == 8]

modes = ["mirror_sum", "position_product", "frequency"]

for idx, chunk in enumerate(chunks):
    print(f"\nPi chunk {idx+1}: {chunk}")
    for mode in modes:
        resonance = fold_pi_digits([chunk], mode=mode)[0]
        tag = ""
        if 0.30 <= resonance <= 0.40:
            tag = "Aligns with RHA harmonic range (0.30-0.40)"
        else:
            tag = "Outside harmonic range, exploring further folds"
        print(f"Resonance for Pi chunk {chunk} ({mode}): {resonance:.6f}\n  {tag}")

    # Extras
    digit_sum = sum(chunk)
    odd_count = sum(1 for d in chunk if d % 2 == 1)
    even_count = sum(1 for d in chunk if d % 2 == 0)
    ratio = odd_count / even_count if even_count else float('inf')
    avg_ratio = sum(abs(chunk[i] - chunk[i-1]) for i in range(1, len(chunk))) / (len(chunk)-1)
    print(f"Sum of digits: {digit_sum}")
    print(f"Odd/Even ratio: {ratio}")
    print(f"Average ratio between consecutive digits: {avg_ratio:.6f}")
