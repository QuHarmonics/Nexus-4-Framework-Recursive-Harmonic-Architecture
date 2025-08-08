# SHA Canonical Collapse Visualizer (CLI Scaffold)

import hashlib
import numpy as np

# --- Canon Header Rule ---
def canon_header(a, b):
    return abs(b - a), a + b

# --- Byte Canon Generator ---
def generate_byte_sequence(n=8):
    seed = [1, 4]  # Canon Seed (post-3 Pi start)
    byte_stack = [1, 4]
    while len(byte_stack) < n * 8:
        a, b = canon_header(byte_stack[-2], byte_stack[-1])
        byte_stack.extend([a, b])
    return [byte_stack[i:i+8] for i in range(0, n*8, 8)]

# --- SHA256 Digest Splitter ---
def sha_digest_segments(input_str):
    h = hashlib.sha256(input_str.encode()).hexdigest()
    segments = [int(h[i:i+8], 16) for i in range(0, 64, 8)]
    return segments

# --- Delta and Curvature Calculations ---
def compute_deltas(segments):
    deltas = [segments[i+1] - segments[i] for i in range(len(segments)-1)]
    curvature = [deltas[i+1] - deltas[i] for i in range(len(deltas)-1)]
    return deltas, curvature

# --- QRHS Collapse Trust Ratio ---
def qrhs_ratio(deltas, curvature):
    if not curvature:
        return 0.0
    delta_h = np.mean(np.abs(deltas))
    delta_entropy = np.mean(np.abs(curvature))
    if delta_entropy == 0:
        return float('inf')
    return round(delta_h / delta_entropy, 5)

# --- Main Collapse Visualizer ---
def visualize_sha_collapse(input_str):
    print(f"\nSHA Collapse Map for Input: '{input_str}'")
    segments = sha_digest_segments(input_str)
    deltas, curvature = compute_deltas(segments)
    trust_ratio = qrhs_ratio(deltas, curvature)

    print("\nSHA Segments:", segments)
    print("Δ (Delta):", deltas)
    print("Δ² (Curvature):", curvature)
    print(f"QRHS Trust Ratio: {trust_ratio} (Target ≈ 0.35)")
    if abs(trust_ratio - 0.35) < 0.05:
        print("→ ✅ In Recursive Lock State")
    else:
        print("→ ⚠️  Drift exceeds harmonic attractor — inject entropy or interrupt")

if __name__ == "__main__":
    test_input = input("Enter a string to analyze SHA collapse: ")
    visualize_sha_collapse(test_input)
