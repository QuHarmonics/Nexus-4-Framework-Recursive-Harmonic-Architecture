import math
import hashlib
from mpmath import mp
from sympy import primerange
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter
import pandas as pd
import os

# --- Parameters ---
depth = 1000000           # Number of Pi digits to extract
max_n = 256               # Maximum leg length for triangle sides
angle_range = (0.34, 0.365) # Resonant angle window in radians
twin_prime_limit = 1000000 # Upper bound for twin prime search
match_window = 2          # How close a Pi index must be to a twin prime

# Set precision to extract the specified number of Pi digits
mp.dps = depth
pi_str = str(mp.pi())[2:2 + depth]  # Skip '3.' and take 'depth' digits
pi_digits = [int(d) for d in pi_str]  # Convert string digits to integers

# --- Triangle Computation ---
def compute_triangle(a, b):
    c = math.sqrt(a**2 + b**2)
    alpha = math.atan(b / a)
    beta = math.atan(a / b)
    height = (a * b) / c
    return {'a': a, 'b': b, 'alpha': alpha, 'beta': beta, 'height': height}

# --- Twin Prime Generation ---
def get_twin_primes(limit):
    primes = list(primerange(2, limit))
    twin_primes = [(p, primes[i + 1]) for i, p in enumerate(primes[:-1]) if primes[i + 1] - p == 2]
    return twin_primes

# --- SHA Curvature Feedback ---
def sha_resonance_adjustment(chunk):
    digits = [int(d) for d in chunk]
    harmonic = sum(digits[:4]) / (sum(digits[4:]) + 1e-6)
    return round(harmonic, 4)

# --- Hash to Pi Index ---
def hash_triangle_to_pi_index(a, b, max_index):
    h = hashlib.sha256(f"{a}:{b}".encode()).hexdigest()
    index = int(h, 16) % max_index
    pi_chunk = pi_digits[index:index + 8]
    return index, ''.join(map(str, pi_chunk))

# --- Core Discovery ---
resonant_triangles = []
twin_primes = get_twin_primes(twin_prime_limit)

for a in range(1, max_n + 1):
    for b in range(1, max_n + 1):
        tri = compute_triangle(a, b)
        if angle_range[0] <= tri['alpha'] <= angle_range[1] or angle_range[0] <= tri['beta'] <= angle_range[1]:
            pi_index, pi_chunk = hash_triangle_to_pi_index(a, b, depth)
            for p, q in twin_primes:
                if abs(pi_index - p) < match_window or abs(pi_index - q) < match_window:
                    resonant_triangles.append({
                        'triangle': f"{a},{b}",
                        'alpha': tri['alpha'],
                        'beta': tri['beta'],
                        'height': tri['height'],
                        'pi_index': pi_index,
                        'pi_chunk': pi_chunk,
                        'harmonic_factor': sha_resonance_adjustment(pi_chunk),
                        'twin_prime': (p, q)
                    })
                    break

# Frequency Analysis for recursive attractors
chunk_counter = Counter([res['pi_chunk'] for res in resonant_triangles])
for res in resonant_triangles:
    res['shared_hits'] = chunk_counter[res['pi_chunk']]

# --- Graph Construction ---
G = nx.Graph()
for res in resonant_triangles:
    tri_node = res['triangle']
    tp_node = f"TP:{res['twin_prime'][0]},{res['twin_prime'][1]}"
    G.add_node(tri_node, type='triangle', curvature=res['harmonic_factor'])
    G.add_node(tp_node, type='twin_prime')
    G.add_edge(tri_node, tp_node)

# --- Visualization ---
pos = nx.spring_layout(G, seed=42)
colors = ['skyblue' if G.nodes[n]['type'] == 'triangle' else 'salmon' for n in G.nodes]

plt.figure(figsize=(18, 12))
nx.draw_networkx(G, pos, with_labels=True, node_color=colors, font_size=6, node_size=300)
plt.title("Resonant Triangles and Twin Prime Network")
plt.axis('off')
plt.show()

# Export to CSV for external tools
df = pd.DataFrame(resonant_triangles)
csv_path = "/mnt/data/resonant_triangles_output.csv"
df.to_csv(csv_path, index=False)
csv_path
