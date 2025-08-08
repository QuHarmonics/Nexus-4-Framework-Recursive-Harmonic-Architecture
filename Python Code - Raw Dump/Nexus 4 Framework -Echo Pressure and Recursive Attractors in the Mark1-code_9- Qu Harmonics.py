import math
import hashlib
from mpmath import mp
from sympy import primerange
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# Parameters
depth = 1000000000  # Number of Pi digits to generate
max_n = 32       # Max value for triangle sides a and b
angle_range = (0.345, 0.365)  # Resonance window

# Generate Pi digits
mp.dps = depth
pi_str = str(mp.pi())[2:depth + 2]
pi_digits = [int(d) for d in pi_str if d.isdigit()]

# Twin primes
def get_twin_primes(max_n):
    primes = list(primerange(2, max_n))
    return [(primes[i], primes[i + 1]) for i in range(len(primes) - 1) if primes[i + 1] - primes[i] == 2]

twin_primes = get_twin_primes(depth)

# Compute triangle properties
def compute_triangle(a, b):
    c = math.sqrt(a**2 + b**2)
    alpha = math.atan(b / a)
    beta = math.atan(a / b)
    height = (a * b) / c
    return {'a': a, 'b': b, 'alpha': alpha, 'beta': beta, 'height': height}

# Map to Pi index
def map_to_pi_index(tri, pi_digits, max_index):
    input_str = f"{tri['a']}:{tri['b']}".encode()
    hash_val = hashlib.sha256(input_str).hexdigest()
    index = int(hash_val, 16) % max_index
    pi_chunk = pi_digits[index:index+8] if index + 8 <= len(pi_digits) else pi_digits[index:]
    return index, pi_chunk

# Find resonant triangles
resonant_triangles = []
for a in range(1, max_n + 1):
    for b in range(1, max_n + 1):
        tri = compute_triangle(a, b)
        if angle_range[0] <= tri['alpha'] <= angle_range[1] or angle_range[0] <= tri['beta'] <= angle_range[1]:
            pi_index, pi_chunk = map_to_pi_index(tri, pi_digits, depth)
            for p, q in twin_primes:
                if abs(pi_index - p) < 10 or abs(pi_index - q) < 10:
                    resonant_triangles.append({
                        'a': tri['a'], 'b': tri['b'],
                        'alpha': tri['alpha'], 'beta': tri['beta'],
                        'height': tri['height'], 'pi_index': pi_index,
                        'pi_chunk': ''.join(map(str, pi_chunk)),
                        'twin_prime': f"{p},{q}"
                    })
                    break

# Create network graph
G = nx.Graph()
for res in resonant_triangles:
    tri_node = f"{res['a']},{res['b']}"
    twin_node = f"TP:{res['twin_prime']}"
    G.add_node(tri_node, type='triangle')
    G.add_node(twin_node, type='twin_prime')
    G.add_edge(tri_node, twin_node)

# Draw graph
plt.figure(figsize=(14, 12))
pos = nx.spring_layout(G, seed=42, k=0.15)
triangle_nodes = [n for n, attr in G.nodes(data=True) if attr['type'] == 'triangle']
twin_nodes = [n for n, attr in G.nodes(data=True) if attr['type'] == 'twin_prime']

nx.draw_networkx_nodes(G, pos, nodelist=triangle_nodes, node_color='skyblue', label='Triangles', node_size=100)
nx.draw_networkx_nodes(G, pos, nodelist=twin_nodes, node_color='salmon', label='Twin Primes', node_size=120)
nx.draw_networkx_edges(G, pos, alpha=0.5)
nx.draw_networkx_labels(G, pos, font_size=12)

plt.title("Resonant Triangles and Twin Prime Network")
plt.axis('off')
plt.legend()
plt.tight_layout()
plt.show()
