from mpmath import mp

def fold_pi_digits(chunks, mode="position_product", factor=2):
    resonances = []
    for chunk_idx, chunk in enumerate(chunks):
        if mode == "mirror_sum":
            left, right = chunk[:4], chunk[4:]
            sums = [left[i] + left[3-i] for i in range(2)]
            resonance = sum(sums) / (9 * 2)
        elif mode == "position_product":
            start_pos = 2 + chunk_idx * 8
            positions = list(range(start_pos, start_pos + 8))
            products = [d * p for d, p in zip(chunk, positions)]
            resonance = sum(products) / (9 * max(positions))
        elif mode == "frequency":
            two_positions = [i+1 for i, d in enumerate(chunk) if d == 2]
            gaps = [two_positions[i+1] - two_positions[i] for i in range(len(two_positions)-1)] if len(two_positions) > 1 else [0]
            resonance = sum(gaps) / (8 * 2)
        else:
            digit_str = ''.join(map(str, chunk))
            decimal = int(digit_str) / 100000000
            resonance = decimal
            for _ in range(3):
                resonance = (resonance * factor) % 1
        resonances.append(resonance)
    return resonances

# Get 10,000 Pi digits after decimal
mp.dps = 10010
pi_str = str(mp.pi)[2:]  # Correct: str(mp.pi) is callable

pi_digits = [int(d) for d in pi_str[:100000]]
chunks = [pi_digits[i:i+8] for i in range(0, len(pi_digits), 8) if len(pi_digits[i:i+8]) == 8]

def odd_even_ratio(chunk):
    odd = sum(1 for d in chunk if d % 2 == 1)
    even = sum(1 for d in chunk if d % 2 == 0)
    return odd / even if even else float('inf')

def avg_ratio_consecutive(chunk):
    return sum(abs(chunk[i] - chunk[i-1]) for i in range(1, len(chunk))) / (len(chunk)-1)

# Print detailed results for N chunks
N = 10000  # or whatever number you want
modes = ["mirror_sum", "position_product", "frequency"]
for idx, chunk in enumerate(chunks[:N]):
    print(f"\nPi chunk {idx+1}: {chunk}")
    for mode in modes:
        r = fold_pi_digits([chunk], mode=mode)[0]
        label = "<- aligns with RHA" if 0.30 <= r <= 0.40 else ""
        print(f"Resonance for Pi chunk {chunk} ({mode}): {r:.6f} {label}")
    print(f"Sum of digits: {sum(chunk)}")
    print(f"Odd/Even ratio: {odd_even_ratio(chunk)}")
    print(f"Average ratio between consecutive digits: {avg_ratio_consecutive(chunk):.6f}")
