import math
import hashlib
from mpmath import mp
from sympy import primerange

# Compute triangle properties: angles and height
def compute_triangle(a, b):
    c = math.sqrt(a**2 + b**2)  # Hypotenuse
    alpha = math.atan(b / a)    # Angle opposite b
    beta = math.atan(a / b)     # Angle opposite a
    height = (a * b) / c        # Height (area / hypotenuse)
    return {'a': a, 'b': b, 'alpha': alpha, 'beta': beta, 'height': height}

# Generate twin prime pairs up to max_n
def get_twin_primes(max_n):
    primes = list(primerange(2, max_n))
    twin_primes = []
    for i in range(len(primes) - 1):
        if primes[i + 1] - primes[i] == 2:
            twin_primes.append((primes[i], primes[i + 1]))
    return twin_primes

# Map triangle to a Pi index using SHA-256
def map_to_pi_index(tri, pi_digits, max_index):
    input_str = f"{tri['a']}:{tri['b']}".encode()  # Unique string for (a, b)
    hash_val = hashlib.sha256(input_str).hexdigest()
    index = int(hash_val, 16) % max_index  # Modulo to fit within Pi digits
    start = index
    end = min(len(pi_digits), index + 8)   # Grab a chunk of 8 digits
    pi_chunk = ''.join(map(str, pi_digits[start:end]))
    return index, pi_chunk

# Main function to find resonant triangles
def generate_resonant_triangles(depth, max_n):
    # Get Pi digits to specified depth
    mp.dps = depth
    pi_str = str(mp.pi())[2:depth + 2]  # Skip "3." and get digits
    pi_digits = [int(d) for d in pi_str]

    # Get twin primes
    twin_primes = get_twin_primes(max_n)
    resonant_triangles = []

    # Search for resonant triangles
    for a in range(1, max_n + 1):
        for b in range(1, max_n + 1):
            tri = compute_triangle(a, b)
            # Check if alpha or beta is in [0.34, 0.36]
            if 0.34 <= tri['alpha'] <= 0.36 or 0.34 <= tri['beta'] <= 0.36:
                pi_index, pi_chunk = map_to_pi_index(tri, pi_digits, depth)
                # Check proximity to twin primes
                for p, q in twin_primes:
                    if abs(pi_index - p) < 10 or abs(pi_index - q) < 10:
                        resonant_triangles.append({
                            'a': tri['a'],
                            'b': tri['b'],
                            'alpha': tri['alpha'],
                            'beta': tri['beta'],
                            'height': tri['height'],
                            'pi_index': pi_index,
                            'pi_chunk': pi_chunk,
                            'twin_prime': (p, q)
                        })
                        break  # Move to next triangle after finding a match
    return resonant_triangles

# Print results
def print_resonant_triangles(resonant_triangles):
    print(f"Found {len(resonant_triangles)} resonant triangles:")
    print("Resonant Triangles with Alpha/Beta in [0.34, 0.36] rad:")
    for res in resonant_triangles[:10]:  # Limit to first 10 for brevity
        print(f"a={res['a']}, b={res['b']}, alpha={res['alpha']:.6f} rad, "
              f"beta={res['beta']:.6f} rad, height={res['height']:.6f}, "
              f"pi_index={res['pi_index']}, pi_chunk={res['pi_chunk']}, "
              f"twin_prime={res['twin_prime']}")

# Run the experiment
depth = 10000  # Number of Pi digits (adjust as needed)
max_n = 200    # Max value for a and b (adjust as needed)
resonant_triangles = generate_resonant_triangles(depth, max_n)
print_resonant_triangles(resonant_triangles)
