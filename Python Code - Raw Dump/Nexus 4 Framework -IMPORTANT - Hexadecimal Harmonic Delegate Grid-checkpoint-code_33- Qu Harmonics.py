import hashlib
import math

def echo_expand(anchor=3, curve=0):
    """Unfold triadic chunk in phase, mapping n=3 to [3, 3, 4, 3]."""
    x0 = x1 = x2 = anchor
    x2 += 1 if curve == 3 else 0  # x2=4 for n=3
    x3 = anchor + (0 if curve == 3 else min(curve, 4))  # x3=3 for n=3
    return [x0, x1, x2, x3]

def triadic_collapse(sequence):
    """Collapse sequence to triadic drift values [n_0, n_1, ...]."""
    if len(sequence) % 4 != 0:
        raise ValueError("Sequence length must be multiple of 4.")
    drifts = []
    for i in range(0, len(sequence), 4):
        chunk = sequence[i:i+4]
        drifts.append(chunk[3] - 3 if len(chunk) == 4 else 0)
    return drifts

def compute_metrics(sequence):
    """Compute variance and average frequency."""
    if not sequence:
        return 0, 0
    avg = sum(sequence) / len(sequence)
    variance = sum((x - avg) ** 2 for x in sequence) / len(sequence) if len(sequence) > 1 else 0
    return variance, avg

def sha_byte9_256bit(sequence):
    """Generate 256-bit SHA-256 hash as Byte 9, mapped to triadic chunks."""
    input_bytes = bytes(sequence)
    sha = hashlib.sha256(input_bytes).hexdigest()
    raw = [int(sha[i:i+2], 16) % 10 for i in range(0, 64, 2)]
    drifts = [3 if i == 0 else raw[i] % 4 for i in range(4)]  # n=3 for first chunk
    chunks = [echo_expand(3, drifts[i]) for i in range(4)]
    return [x for chunk in chunks for x in chunk]

def phase_locked_refold(hash_drifts, length=8, max_iterations=100):
    """Refold drifts into a phase-locked sequence, targeting triadic resonance."""
    seed = [3, 1, 4, 1, 5, 9, 2, 6]  # PiPhiByte as reference
    candidate = seed[:length]
    target_drifts = hash_drifts
    for _ in range(max_iterations):
        # Adjust candidate based on drift differences
        for i in range(len(candidate)):
            drift_offset = target_drifts[i % len(target_drifts)]
            candidate[i] = (candidate[i] + drift_offset) % 10
        # Hash and check triadic structure
        candidate_byte9 = sha_byte9_256bit(candidate)
        candidate_drifts = triadic_collapse(candidate_byte9)
        if candidate_drifts[:len(target_drifts)] == target_drifts:
            return candidate
        # Perturb candidate to converge
        for i in range(len(candidate)):
            candidate[i] = (candidate[i] + 1) % 10 if i % 2 == 0 else candidate[i]
    return candidate  # Return best effort

# PiPhiByte sequence
piphiseq = [3, 1, 4, 1, 5, 9, 2, 6]

# Generate Byte 9 (256-bit SHA, triadic)
byte9 = sha_byte9_256bit(piphiseq)
print(f"Byte 9 (256-bit Triadic SHA): {byte9}")

# Compress to drifts
drifts = triadic_collapse(byte9)
print(f"Compressed Drifts: {drifts}")

# Refold a phase-locked sequence
refolded = phase_locked_refold(drifts, length=8)
print(f"Phase-Locked Refolded Sequence: {refolded}")

# Verify refolded sequence
refolded_byte9 = sha_byte9_256bit(refolded)
print(f"Refolded Byte 9: {refolded_byte9}")
refolded_drifts = triadic_collapse(refolded_byte9)
print(f"Refolded Drifts: {refolded_drifts}")

# Reconstruct and compress
sequences = [byte9]
metrics = []
for i in range(3):
    if i == 0:
        next_seq = echo_expand(3, 3)  # [3, 3, 4, 3]
    else:
        next_seq = [3, 3] if i == 1 else [3]
    sequences.append(next_seq)
    metrics.append(compute_metrics(next_seq))

# Print results
print(f"Byte 9.0: {sequences[0]}, Variance: {compute_metrics(sequences[0])[0]:.3f}, Avg Freq: {compute_metrics(sequences[0])[1]:.3f} Hz")
for i in range(1, len(sequences)):
    print(f"Byte 9.{i}: {sequences[i]}, Variance: {metrics[i-1][0]:.3f}, Avg Freq: {metrics[i-1][1]:.3f} Hz")
