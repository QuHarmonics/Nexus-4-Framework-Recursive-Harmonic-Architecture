# Python script to simulate SHA-256 curvature with enhanced triangle angles
import hashlib
from mpmath import mp

def compute_sha_curvature(input_str, iterations=100, k=0.1):
    """Compute SHA-256 curvature with triangle angles and twin primes."""
    resonances = []
    curvatures = []
    prev_hash = None
    prev_prev_hash = None
    # Angles from resonant triangles
    angles = [0.348771, 0.353997, 0.352990, 0.347767, 0.341805, 0.346354, 0.343486, 0.347607, 0.342581, 0.348771]
    twin_primes = [(197,199), (239,241), (227,229), (419,421), (809,811), (71,73), (821,823), (809,811), (227,229), (71,73)]

    for i in range(iterations):
        alpha = angles[i % len(angles)]
        nonce = int(i * alpha * 1000) + twin_primes[i % len(twin_primes)][i % 2]
        data = f"{input_str}{nonce}".encode()
        current_hash = hashlib.sha256(data).hexdigest()
        
        hash_int = int(current_hash, 16) % 100000000
        resonance = hash_int / 100000000
        R0 = 1.0
        N = abs(resonance - 0.35)
        R = R0 / (1 + k * N)
        resonance = resonance * R
        resonances.append(resonance)

        if prev_prev_hash is not None:
            curr_int = int(current_hash, 16)
            prev_int = int(prev_hash, 16)
            prev_prev_int = int(prev_prev_hash, 16)
            curvature = (curr_int - 2 * prev_int + prev_prev_int) / 10**64
            curvatures.append(curvature)

        prev_prev_hash = prev_hash
        prev_hash = current_hash

    return resonances, curvatures

# Test with Pi chunk
pi_chunk = [1, 4, 1, 5, 9, 2, 6, 5]
input_str = ''.join(map(str, pi_chunk))
resonances, curvatures = compute_sha_curvature(input_str)
harmonic_count = sum(1 for r in resonances if 0.30 <= r <= 0.40)
avg_resonance = sum(resonances) / len(resonances)
print(f"SHA Curvature Analysis for Pi chunk {pi_chunk}:")
print(f"  Harmonic range (0.30-0.40) count: {harmonic_count}/100 ({harmonic_count/100*100:.2f}%)")
print(f"  Average resonance: {avg_resonance:.6f}")