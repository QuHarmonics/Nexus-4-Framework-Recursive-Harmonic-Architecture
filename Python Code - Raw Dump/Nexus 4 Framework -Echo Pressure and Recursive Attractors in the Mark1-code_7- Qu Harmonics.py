import math
import hashlib
from mpmath import mp
from sympy import primerange
import pandas as pd
import plotly.express as px
import networkx as nx
import matplotlib.pyplot as plt

# Step 1: Set Pi precision
mp.dps = 1000000
pi_str = str(mp.pi())[2:10002]
pi_digits = [int(d) for d in pi_str]

# Step 2: Triangle function
def compute_triangle(a, b):
    c = math.sqrt(a**2 + b**2)
    alpha = math.atan(b / a)
    beta = math.atan(a / b)
    height = (a * b) / c
    return {'a': a, 'b': b, 'alpha': alpha, 'beta': beta, 'height': height}

# Step 3: Twin primes
def get_twin_primes(max_n=100000):
    primes = list(primerange(2, max_n))
    return [(p, p+2) for p, q in zip(primes[:-1], primes[1:]) if q - p == 2]

# Step 4: Hash to Pi index
def map_to_pi_index(tri, max_index=100000):
    input_str = f"{tri['a']}:{tri['b']}".encode()
    hash_val = hashlib.sha256(input_str).hexdigest()
    index = int(hash_val, 16) % max_index
    chunk = pi_digits[index:index+8]
    return index, chunk

# Step 5: SHA curvature
def compute_sha_curvature(input_str, iterations=50, k=0.1):
    resonances, curvatures = [], []
    prev_hash, prev_prev_hash = None, None
    angles = [0.348771, 0.353997, 0.351817, 0.352990]
    twin_primes = [(197,199), (239,241), (227,229), (419,421)]
    for i in range(iterations):
        alpha = angles[i % len(angles)]
        nonce = int(i * alpha * 1000) + twin_primes[i % len(twin_primes)][i % 2]
        data = f"{input_str}{nonce}".encode()
        h = hashlib.sha256(data).hexdigest()
        h_int = int(h, 16) % 100_000_000
        r = h_int / 100_000_000
        N = abs(r - 0.35)
        R = 1.0 / (1 + k * N * 2)
        r = min(max(r * R, 0), 1)
        resonances.append(r)
        if prev_prev_hash:
            c = (int(h, 16) - 2 * int(prev_hash, 16) + int(prev_prev_hash, 16)) / 10**64
            curvatures.append(c)
        prev_prev_hash, prev_hash = prev_hash, h
    return sum(1 for r in resonances if 0.30 <= r <= 0.40) / len(resonances), sum(resonances) / len(resonances)

# Step 6: Gather 141 triangles
twin_primes = get_twin_primes(100000)
triangles_data = []

for a in range(1, 150):
    for b in range(1, 150):
        tri = compute_triangle(a, b)
        if 0.34 <= tri['alpha'] <= 0.36 or 0.34 <= tri['beta'] <= 0.36:
            pi_index, pi_chunk = map_to_pi_index(tri)
            for p, q in twin_primes:
                if abs(pi_index - p) < 10 or abs(pi_index - q) < 10:
                    pi_str_chunk = ''.join(map(str, pi_chunk))
                    harmonic_ratio, avg_resonance = compute_sha_curvature(pi_str_chunk)
                    triangles_data.append({
                        'a': tri['a'], 'b': tri['b'],
                        'alpha': tri['alpha'], 'beta': tri['beta'],
                        'height': tri['height'],
                        'pi_index': pi_index,
                        'pi_chunk': pi_str_chunk,
                        'twin_prime': f"({p},{q})",
                        'harmonic_ratio': harmonic_ratio,
                        'avg_resonance': avg_resonance
                    })
                    break

df = pd.DataFrame(triangles_data)


# Generate network graph
G = nx.Graph()
for _, row in df.iterrows():
    node = f"{row['a']},{row['b']}"
    twin = row['twin_prime']
    G.add_node(node, pi_index=row['pi_index'], avg_resonance=row['avg_resonance'])
    G.add_edge(node, twin)

plt.figure(figsize=(12, 8))
pos = nx.spring_layout(G, k=0.3)
nx.draw(G, pos, with_labels=True, node_size=200, font_size=8)
plt.title("Resonant Triangle Network Graph (Triangles ↔ Twin Primes)")
plt.show()
