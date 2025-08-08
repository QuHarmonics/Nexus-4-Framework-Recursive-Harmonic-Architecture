from mpmath import mp

def fold_pi_digits(chunks, mode="position_product", factor=2):
    """
    Explore harmonic resonance in 8-digit Pi chunks.
    Input: chunks (list of lists, each with 8 integers)
    Modes: 'mirror_sum', 'position_product', 'frequency'
    Output: List of resonance values for each chunk
    """
    resonances = []
    for chunk_idx, chunk in enumerate(chunks):
        if mode == "mirror_sum":
            # Sum mirrored halves (left + right)
            left, right = chunk[:4], chunk[4:]
            sums = [left[i] + right[3-i] for i in range(4)]
            resonance = sum(sums) / (9 * 4)  # Normalize by (max digit 9 * 4 pairs)
        elif mode == "position_product":
            # Digit-position products (positions always 1-8 for chunk length 8)
            positions = list(range(1, 9))
            products = [d * p for d, p in zip(chunk, positions)]
            resonance = sum(products) / (9 * 8 * 4)  # Normalize by (max digit * max pos * half len)
        elif mode == "frequency":
            # Count digit 2 frequency/gaps in chunk
            two_positions = [i+1 for i, d in enumerate(chunk) if d == 2]
            gaps = [two_positions[i+1] - two_positions[i] for i in range(len(two_positions)-1)] if len(two_positions) > 1 else [0]
            resonance = sum(gaps) / (8 * 2)  # Normalize by max gap
        else:
            # Recursive fold on decimal (rarely used)
            digit_str = ''.join(map(str, chunk))
            decimal = int(digit_str) / 10**8
            resonance = decimal
            for _ in range(3):
                resonance = (resonance * factor) % 1
        resonances.append(resonance)
    return resonances

# Get first 10,000 digits of pi after "3."
mp.dps = 10010  # Extra to avoid rounding artifacts
pi_digits_str = str(mp.pi)[2:10002]  # Skip '3.', take 10,000 digits only
pi_digits = [int(d) for d in pi_digits_str]

# Chunk into groups of 8
chunks = [pi_digits[i:i+8] for i in range(0, len(pi_digits), 8)]
chunks = [c for c in chunks if len(c) == 8]  # Only full 8-digit chunks

# Test all modes, print only stats
modes = ["mirror_sum", "position_product", "frequency"]
for mode in modes:
    resonances = fold_pi_digits(chunks, mode=mode)
    harmonic_count = sum(1 for r in resonances if 0.30 <= r <= 0.40)
    avg_resonance = sum(resonances) / len(resonances)
    print(f"Mode: {mode}")
    print(f"  Harmonic range (0.30-0.40) count: {harmonic_count}/{len(chunks)} ({harmonic_count/len(chunks)*100:.2f}%)")
    print(f"  Average resonance: {avg_resonance:.6f}")
