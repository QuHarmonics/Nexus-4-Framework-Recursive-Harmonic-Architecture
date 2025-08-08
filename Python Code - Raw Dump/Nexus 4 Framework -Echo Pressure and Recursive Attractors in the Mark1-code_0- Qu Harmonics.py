import math
import hashlib
from mpmath import mp
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter

# Set precision for π digits
mp.dps = 1000000
pi_str = str(mp.pi())[2:1000002]
pi_digits = [int(d) for d in pi_str if d.isdigit()][:10000]

# Constants
HARMONIC_TARGET = 0.35
TWIN_WINDOW = 10

# Twin primes
def get_twin_primes(max_n=1000):
    def is_prime(n):
        if n < 2: return False
        for i in range(2, int(n**0.5)+1):
            if n % i == 0: return False
        return True
    return [(p, p+2) for p in range(2, max_n) if is_prime(p) and is_prime(p+2)]

twin_primes = get_twin_primes(1000)

# Triangle computation
def compute_triangle(a, b):
    c = math.sqrt(a**2 + b**2)
    alpha = math.atan(b / a)
    beta = math.atan(a / b)
    area = 0.5 * a * b
    perimeter = a + b + c
    height = b * a / c
    return {
        'a': a, 'b': b, 'c': c, 'alpha': alpha, 'beta': beta,
        'area': area, 'perimeter': perimeter, 'height': height
    }

# SHA to π index mapping
def map_to_pi_index(tri):
    input_str = f"{tri['a']}:{tri['b']}".encode()
    hash_val = hashlib.sha256(input_str).hexdigest()
    index = int(hash_val, 16) % 10000
    return index, pi_digits[index:index+8]

# Harmonic ratio H(t) calculation
def compute_H(pi_chunk):
    P = sum(pi_chunk[:4])
    A = sum(pi_chunk)
    return P / A if A != 0 else 0

# Curvature Δ²H(t)
def curvature_delta(H_list):
    deltas = [H_list[i] - H_list[i-1] for i in range(1, len(H_list))]
    d2 = [deltas[i] - deltas[i-1] for i in range(1, len(deltas))]
    return d2

# KHRC correction pass
def khrc(H_val, k=1.0, N=0.01):
    return H_val / (1 + k * abs(N))

# Convergence score
def convergence_score(H_vals, target=0.35, margin=0.05):
    return sum(1 for h in H_vals if abs(h - target) <= margin) / len(H_vals)

# Main process
resonant_triangles = []
H_trajectory = []

for a in range(1, 4096):
    for b in range(1, 4096):
        tri = compute_triangle(a, b)
        if 0.34 <= tri['alpha'] <= 0.36 or 0.34 <= tri['beta'] <= 0.36:
            pi_index, pi_chunk = map_to_pi_index(tri)
            for p, q in twin_primes:
                if abs(pi_index - p) < TWIN_WINDOW or abs(pi_index - q) < TWIN_WINDOW:
                    H_val = compute_H(pi_chunk)
                    H_trajectory.append(H_val)
                    resonant_triangles.append({
                        'triangle': tri,
                        'pi_index': pi_index,
                        'pi_chunk': pi_chunk,
                        'twin_prime': (p, q),
                        'H': H_val,
                        'KHRC': khrc(H_val)
                    })
                    break

# Display summary
for res in resonant_triangles[:10]:
    tri = res['triangle']
    print(f"a={tri['a']}, b={tri['b']}, alpha={tri['alpha']:.6f}, beta={tri['beta']:.6f}, "
          f"height={tri['height']:.6f}, H={res['H']:.5f}, KHRC={res['KHRC']:.5f}, "
          f"pi_index={res['pi_index']}, twin_prime={res['twin_prime']}, pi_chunk={res['pi_chunk']}")

# Plot H(t) and Δ²H(t)
if H_trajectory:
    plt.figure()
    plt.plot(H_trajectory, label='H(t)', marker='o')
    d2H = curvature_delta(H_trajectory)
    plt.plot(range(2, 2 + len(d2H)), d2H, label='Δ²H(t)', marker='x')
    plt.axhline(y=HARMONIC_TARGET, color='gray', linestyle='--', label='Target H=0.35')
    plt.title("Harmonic Ratio H(t) and Curvature Δ²H(t)")
    plt.xlabel("Iteration")
    plt.ylabel("Value")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# Resonance metrics
resonance_deviation = sum(abs(h - HARMONIC_TARGET) for h in H_trajectory) / len(H_trajectory)
score = convergence_score(H_trajectory)
print(f"\n⎇ Cumulative Resonance Deviation: {resonance_deviation:.5f}")
print(f"★ Harmonic Convergence Score: {score*100:.1f}%")

# Echo pattern detection
chunk_strings = [''.join(map(str, chunk)) for chunk in [r['pi_chunk'] for r in resonant_triangles]]
repeats = Counter(chunk_strings)
echo_patterns = {k: v for k, v in repeats.items() if v > 1}
print(f"↻ Echo Signatures: {echo_patterns}")

# Export dataset
df = pd.DataFrame([{
    'a': r['triangle']['a'],
    'b': r['triangle']['b'],
    'alpha': r['triangle']['alpha'],
    'beta': r['triangle']['beta'],
    'height': r['triangle']['height'],
    'H': r['H'],
    'KHRC': r['KHRC'],
    'pi_index': r['pi_index'],
    'pi_chunk': ''.join(map(str, r['pi_chunk'])),
    'twin_prime': f"{r['twin_prime'][0]},{r['twin_prime'][1]}"
} for r in resonant_triangles])

csv_path = "d://Nexus/resonant_triangles_nexus.csv"
df.to_csv(csv_path, index=False)

# Echo spectrum chart using Plotly (interactive + resizable)
import plotly.graph_objects as go

chunk_counter = Counter(''.join(map(str, r['pi_chunk'])) for r in resonant_triangles)
echo_chunks = {k: v for k, v in chunk_counter.items() if v > 1}

if echo_chunks:
    sorted_chunks = dict(sorted(echo_chunks.items(), key=lambda x: x[1], reverse=True))

    fig = go.Figure(data=[
        go.Bar(
            x=list(sorted_chunks.keys()),
            y=list(sorted_chunks.values()),
            text=list(sorted_chunks.values()),
            textposition='outside',
            marker_color='teal'
        )
    ])

    fig.update_layout(
        title="Echo Signature Spectrum (π-Chunk Frequency ≥ 2)",
        xaxis_title="π-Chunk",
        yaxis_title="Frequency (Echo Pressure)",
        width=1400,   # ← Adjust width here
        height=700,   # ← Adjust height here
        xaxis_tickangle=-45,
        margin=dict(l=40, r=40, t=80, b=150),
        template="plotly_white"
    )

    fig.show()
else:
    print("No π-chunks repeated more than once.")
