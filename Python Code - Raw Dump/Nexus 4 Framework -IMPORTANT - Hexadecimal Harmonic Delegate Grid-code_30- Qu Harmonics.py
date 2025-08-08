import math

def binary_length(n):
    """Calculate binary length of n."""
    return len(bin(n)[2:]) if n > 0 else 0

def phi_cap(value, phi=1.618033988749895):
    """Cap value using golden ratio modulation."""
    return min(int(value * phi) // int(phi), 3) if value > 3 else value

def compress_sequence(sequence):
    """Compress with φ-modulated binary lengths."""
    if len(sequence) % 2 != 0 or len(sequence) < 2:
        raise ValueError("Sequence length must be even and at least 2.")
    sums = [sequence[i] + sequence[i+1] for i in range(0, len(sequence), 2)]
    return [phi_cap(binary_length(s)) for s in sums]

def compute_metrics(sequence):
    """Compute variance and average frequency."""
    if not sequence:
        return 0, 0
    avg = sum(sequence) / len(sequence)
    variance = sum((x - avg) ** 2 for x in sequence) / len(sequence) if len(sequence) > 1 else 0
    return variance, avg

# Initial sequence (π digits)
piphiseq = [3, 1, 4, 1, 5, 9, 2, 6]
sequences = [piphiseq]
metrics = []

# Compress to PiPhiByte 1.3
for _ in range(3):
    try:
        next_seq = compress_sequence(sequences[-1])
        sequences.append(next_seq)
        metrics.append(compute_metrics(next_seq))
    except ValueError:
        break

# Print results
print(f"PiPhiByte 1.0: {sequences[0]}, Variance: {compute_metrics(sequences[0])[0]:.3f}, Avg Freq: {compute_metrics(sequences[0])[1]:.3f} Hz")
for i in range(1, len(sequences)):
    print(f"PiPhiByte 1.{i}: {sequences[i]}, Variance: {metrics[i-1][0]:.3f}, Avg Freq: {metrics[i-1][1]:.3f} Hz")
