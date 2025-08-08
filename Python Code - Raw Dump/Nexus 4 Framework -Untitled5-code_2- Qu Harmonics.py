from mpmath import mp

def fold_pi_digits(chunk, mode="position_product", factor=2):
    if mode == "mirror_sum":
        left, right = chunk[:4], chunk[4:]
        sums = [left[i] + right[3-i] for i in range(4)]
        resonance = sum(sums) / (9 * 4)
    elif mode == "position_product":
        positions = list(range(1, 9))
        products = [d * p for d, p in zip(chunk, positions)]
        resonance = sum(products) / (9 * 8 * 4)
    elif mode == "frequency":
        two_positions = [i+1 for i, d in enumerate(chunk) if d == 2]
        gaps = [two_positions[i+1] - two_positions[i] for i in range(len(two_positions)-1)] if len(two_positions) > 1 else [0]
        resonance = sum(gaps) / (8 * 2)
    else:
        digit_str = ''.join(map(str, chunk))
        decimal = int(digit_str) / 10**8
        resonance = decimal
        for _ in range(3):
            resonance = (resonance * factor) % 1
    return resonance

# User controls how many chunks (8-byte blocks) to process
num_bytes = 10  # Change this to however many "bytes" you want to process

mp.dps = num_bytes * 8 + 2  # Just enough digits
pi_digits_str = str(mp.pi)[2:(num_bytes*8)+2]
pi_digits = [int(d) for d in pi_digits_str]
chunks = [pi_digits[i:i+8] for i in range(0, len(pi_digits), 8)]
chunks = [c for c in chunks if len(c) == 8]

modes = ["mirror_sum", "position_product", "frequency"]

for idx, chunk in enumerate(chunks[:num_bytes]):
    print(f"\nPi chunk {idx+1}: {chunk}")
    for mode in modes:
        resonance = fold_pi_digits(chunk, mode=mode)
        print(f"Resonance for Pi chunk {chunk} ({mode}): {resonance:.6f}")
        if 0.30 <= resonance <= 0.40:
            print("  Aligns with RHA harmonic range (0.30-0.40)")
        else:
            print("  Outside harmonic range, exploring further folds")
