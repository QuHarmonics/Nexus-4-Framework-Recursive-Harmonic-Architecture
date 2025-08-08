
import math

def binary_length(n):
    """Calculate binary length of n."""
    return len(bin(n)[2:]) if n > 0 else 0

def compress_sequence(sequence):
    """Compress sequence by summing adjacent pairs and taking binary lengths."""
    if len(sequence) % 2 != 0 or len(sequence) < 2:
        raise ValueError("Sequence length must be even and at least 2.")
    sums = [sequence[i] + sequence[i+1] for i in range(0, len(sequence), 2)]
    return [binary_length(s) for s in sums]

def pi_digits(n):
    """Return first n digits of π."""
    return [3, 1, 4, 1, 5, 9, 2, 6][:n]

def phi_offset(n, angle_deg=137.5):
    """Generate n phase offsets using golden angle (degrees)."""
    angle_rad = math.radians(angle_deg)  # φ-based angle
    return [(i * angle_rad) % (2 * math.pi) for i in range(n)]

def recursive_byte_model(pi_seq, phi_offsets, H=0.35, max_iter=3):
    """Simulate recursive byte model with π and φ as anchors."""
    sequence = pi_seq[:8]  # Start with π digits
    sequences = [sequence]
    metrics = []

    for i in range(max_iter):
        # Apply φ-based phase offset as frequency modulation
        modulated = [int(s * (1 + H * math.cos(phi_offsets[j % len(phi_offsets)]))) 
                     for j, s in enumerate(sequence)]
        # Cap values to prevent overflow (H=0.35 bias)
        sequence = [min(max(1, x), 9) for x in modulated]
        # Compress
        try:
            sequence = compress_sequence(sequence)
            sequences.append(sequence)
            # Compute metrics
            avg = sum(sequence) / len(sequence)
            variance = sum((x - avg) ** 2 for x in sequence) / len(sequence) if len(sequence) > 1 else 0
            metrics.append((variance, avg))
        except ValueError:
            break

    return sequences, metrics

# Simulate
pi_seq = pi_digits(8)  # [3, 1, 4, 1, 5, 9, 2, 6]
phi_offsets = phi_offset(8)  # Golden angle offsets
sequences, metrics = recursive_byte_model(pi_seq, phi_offsets)

# Print results
for i, (seq, (var, avg)) in enumerate(zip(sequences, metrics + [(0, sum(sequences[-1])/len(sequences[-1]))])):
    print(f"PiPhiByte 1.{i}: {seq}, Variance: {var:.3f}, Avg Freq: {avg:.3f} Hz")
