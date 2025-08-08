import math
import hashlib
from mpmath import mp
from sympy import primerange
import networkx as nx
import matplotlib.pyplot as plt

# --- Parameters ---
depth = 1000000           # Number of Pi digits to extract
max_n = 256               # Maximum leg length for triangle sides
angle_range = (0.345, 0.365) # Resonant angle window in radians
twin_prime_limit = 10000000 # Upper bound for twin prime search
match_window = 2          # How close a Pi index must be to a twin prime

# --- Pi Digit Extraction ---
# Set precision to extract the specified number of Pi digits
mp.dps = depth
pi_str = str(mp.pi())[2:2 + depth]  # Skip '3.' and take 'depth' digits
pi_digits = [int(d) for d in pi_str]  # Convert string digits to integers

# --- Triangle Computation ---
def compute_triangle(a, b):
    """
    Compute properties of a right triangle with legs a and b.
    Returns a dictionary with side lengths, angles, and height.
    """
    c = math.sqrt(a**2 + b**2)  # Hypotenuse
    alpha = math.atan(b / a)    # Angle opposite side b
    beta = math.atan(a / b)     # Angle opposite side a
    height = (a * b) / c        # Height from right angle to hypotenuse
    return {'a': a, 'b': b, 'alpha': alpha, 'beta': beta, 'height': height}

# --- Twin Prime Generation ---
def get_twin_primes(limit):
    """
    Generate a list of twin prime pairs up to the specified limit.
    Twin primes are pairs of primes differing by 2 (e.g., (3, 5)).
    """
    primes = list(primerange(2, limit))  # Generate primes up to limit
    twin_primes = [(p, primes[i + 1]) for i, p in enumerate(primes[:-1]) 
                   if primes[i + 1] - p == 2]
    return twin_primes

# --- Hash to Pi Index ---
def hash_triangle_to_pi_index(a, b, max_index):
    """
    Map triangle sides to an index in Pi digits using SHA-256 hash.
    Returns the index and an 8-digit chunk of Pi starting at that index.
    """
    h = hashlib.sha256(f"{a}:{b}".encode()).hexdigest()  # Hash the sides
    index = int(h, 16) % max_index  # Convert hash to integer and bound it
    pi_chunk = pi_digits[index:index + 8]  # Extract 8 digits from Pi
    return index, ''.join(map(str, pi_chunk))

# --- Core Discovery ---
resonant_triangles = []
twin_primes = get_twin_primes(twin_prime_limit)

# Iterate over all possible triangle side pairs
for a in range(1, max_n + 1):
    for b in range(1, max_n + 1):
        tri = compute_triangle(a, b)
        # Check if either angle falls within the resonant range
        if (angle_range[0] <= tri['alpha'] <= angle_range[1] or 
            angle_range[0] <= tri['beta'] <= angle_range[1]):
            pi_index, pi_chunk = hash_triangle_to_pi_index(a, b, depth)
            # Check proximity to twin primes
            for p, q in twin_primes:
                if (abs(pi_index - p) < match_window or 
                    abs(pi_index - q) < match_window):
                    resonant_triangles.append({
                        'triangle': f"{a},{b}",
                        'alpha': tri['alpha'],
                        'beta': tri['beta'],
                        'height': tri['height'],
                        'pi_index': pi_index,
                        'pi_chunk': pi_chunk,
                        'twin_prime': (p, q)
                    })
                    break  # Stop after finding the first match

# --- Graph Construction ---
G = nx.Graph()

# Add nodes and edges to the graph
for res in resonant_triangles:
    tri_node = res['triangle']  # e.g., "1,2"
    tp_node = f"TP:{res['twin_prime'][0]},{res['twin_prime'][1]}"  # e.g., "TP:3,5"
    G.add_node(tri_node, type='triangle')
    G.add_node(tp_node, type='twin_prime')
    G.add_edge(tri_node, tp_node)

# --- Visualization ---
pos = nx.spring_layout(G, seed=42)  # Compute node positions with fixed seed
colors = ['skyblue' if G.nodes[n]['type'] == 'triangle' else 'salmon' 
          for n in G.nodes]  # Assign colors based on node type

plt.figure(figsize=(18, 12))  # Set figure size
nx.draw_networkx(G, pos, with_labels=True, node_color=colors, 
                 font_size=6, node_size=300)  # Draw the graph
plt.title("Resonant Triangles and Twin Prime Network")  # Add title
plt.axis('off')  # Hide axes
plt.show()  # Display the graph