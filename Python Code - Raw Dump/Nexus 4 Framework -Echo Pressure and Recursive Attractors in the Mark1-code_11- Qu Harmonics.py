import math
import hashlib
from mpmath import mp
from sympy import primerange
import networkx as nx
import matplotlib.pyplot as plt

# Set high precision for Pi digits
depth = 10000000
mp.dps = depth
pi_str = str(mp.pi())[2:depth + 2]  # Skip "3." and get digits
pi_digits = [int(d) for d in pi_str]

# Generate twin prime pairs up to max_n
def get_twin_primes(max_n):
    primes = list(primerange(2, max_n))
    return [(p, q) for p, q in zip(primes, primes[1:]) if q - p == 2]

# Compute triangle properties
def compute_triangle(a, b):
    c = math.sqrt(a**2 + b**2)
    alpha = math.atan(b / a)
    beta = math.atan(a / b)
    height = (a * b) / c
    return {'a': a, 'b': b, 'alpha': alpha, 'beta': beta, 'height': height}

# Map triangle to a Pi index using SHA-256
def map_to_pi_index(tri, pi_digits, max_index):
    input_str = f"{tri['a']}:{tri['b']}".encode()
    hash_val = hashlib.sha256(input_str).hexdigest()
    index = int(hash_val, 16) % max_index
    start = index
    end = min(len(pi_digits), index + 8)
    pi_chunk = ''.join(map(str, pi_digits[start:end]))
    return index, pi_chunk

# Generate and connect the network
def generate_network_graph(depth, max_n):
    twin_primes = get_twin_primes(max_n)
    G = nx.Graph()

    for a in range(1, max_n + 1):
        for b in range(1, max_n + 1):
            tri = compute_triangle(a, b)
            if 0.346 <= tri['alpha'] <= 0.365 or 0.346 <= tri['beta'] <= 0.365:
                pi_index, pi_chunk = map_to_pi_index(tri, pi_digits, depth)
                triangle_label = f"{a}|{b}"
                G.add_node(triangle_label, type='triangle')

                for p, q in twin_primes:
                    if abs(pi_index - p) < 10 or abs(pi_index - q) < 10:
                        tp_label = f"TP:{p},{q}"
                        G.add_node(tp_label, type='prime')
                        G.add_edge(triangle_label, tp_label)
                        break

    return G

# Generate the graph
G = generate_network_graph(depth=1000000, max_n=64)

# Draw the graph
plt.figure(figsize=(18, 12))
pos = nx.spring_layout(G, seed=42)
colors = ['skyblue' if G.nodes[n]['type'] == 'triangle' else 'salmon' for n in G.nodes]
nx.draw_networkx(G, pos, with_labels=True, node_color=colors, font_size=6, node_size=300)
plt.title("Resonant Triangles and Twin Prime Network (Real Data Only)")
plt.axis('off')
plt.tight_layout()
plt.show()
