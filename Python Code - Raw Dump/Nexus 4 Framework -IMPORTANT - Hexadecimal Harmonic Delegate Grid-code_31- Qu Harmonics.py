
import hashlib
import math

def echo_expand(anchor=3, curve=0, phi=1.618033988749895):
    x0 = anchor
    x3 = anchor + int(curve * phi)
    midpoint = (x0 + x3) // 2
    delta = abs(x3 - x0) // 2
    x1 = midpoint - delta // 2
    x2 = midpoint + delta // 2
    return [x0, x1, x2, min(x3, 7)]

def triadic_collapse(sequence):
    if len(sequence) % 4 != 0:
        raise ValueError("Sequence length must be multiple of 4.")
    drifts = []
    for i in range(0, len(sequence), 4):
        chunk = sequence[i:i+4]
        drifts.append(chunk[3] - 3 if len(chunk) == 4 else 0)
    return drifts

def compute_metrics(sequence):
    if not sequence:
        return 0, 0
    avg = sum(sequence) / len(sequence)
    variance = sum((x - avg) ** 2 for x in sequence) / len(sequence) if len(sequence) > 1 else 0
    return variance, avg

def sha_byte9_256bit(sequence):
    input_bytes = bytes(sequence)
    sha = hashlib.sha256(input_bytes).hexdigest()
    # Use full 256 bits (64 hex chars), split into 16 chunks
    raw = [int(sha[i:i+2], 16) % 10 for i in range(0, 64, 2)]
    drifts = [(raw[i] % 4) + 1 if i % 4 == 0 else raw[i] % 4 for i in range(4)]  # 4 chunks
    chunks = [echo_expand(3, drifts[i]) for i in range(4)]
    return [x for chunk in chunks for x in chunk]  # Flatten

# PiPhiByte sequence
piphiseq = [3, 1, 4, 1, 5, 9, 2, 6]

# Generate Byte 9 (256-bit SHA, triadic)
byte9 = sha_byte9_256bit(piphiseq)
print(f"Byte 9 (256-bit Triadic SHA): {byte9}")

# Compress to drifts
drifts = triadic_collapse(byte9)
print(f"Compressed Drifts: {drifts}")

# Reconstruct and compress
sequences = [byte9]
metrics = []
for i in range(3):
    if i == 0:
        next_seq = echo_expand(3, drifts[0])
    else:
        next_seq = [3, 3] if i == 1 else [3]
    sequences.append(next_seq)
    metrics.append(compute_metrics(next_seq))

# Print results
print(f"Byte 9.0: {sequences[0]}, Variance: {compute_metrics(sequences[0])[0]:.3f}, Avg Freq: {compute_metrics(sequences[0])[1]:.3f} Hz")
for i in range(1, len(sequences)):
    print(f"Byte 9.{i}: {sequences[i]}, Variance: {metrics[i-1][0]:.3f}, Avg Freq: {metrics[i-1][1]:.3f} Hz")
