import math
import hashlib
from mpmath import mp
from sympy import primerange

def compute_triangle(a, b):
    """Compute triangle properties for a, b pair."""
    c = math.sqrt(a**2 + b**2)
    alpha = math.atan(b/a)
    beta = math.atan(a/b)
    area = 0.5 * a * b
    perimeter = a + b + c
    height = b * a / c
    return {'a': a, 'b': b, 'c': c, 'alpha': alpha, 'beta': beta, 'area': area, 'perimeter': perimeter, 'height': height}

def get_twin_primes(max_n):
    """Generate twin primes up to max_n."""
    primes = list(primerange(2, max_n))
    twin_primes = []
    for i in range(len(primes) - 1):
        if primes[i+1] - primes[i] == 2:
            twin_primes.append((primes[i], primes[i+1]))
    return twin_primes

def map_to_pi_index(tri, pi_digits, max_index):
    """Map triangle to pi digit index via SHA-256."""
    input_str = f"{tri['a']}:{tri['b']}".encode()
    hash_val = hashlib.sha256(input_str).hexdigest()
    index = int(hash_val, 16) % max_index
    start = max(0, index)
    end = min(len(pi_digits), index + 8)
    return index, pi_digits[start:end]

def generate_resonant_triangles(depth, max_n):
    """Generate resonant triangles up to depth and max_n."""
    # Set precision for pi
    mp.dps = depth
    pi_str = str(mp.pi())[2:depth+2]  # Skip '3.'
    pi_digits = [int(d) for d in pi_str if d.isdigit()]

    twin_primes = get_twin_primes(max_n)
    resonant_triangles = []

    for a in range(1, max_n + 1):
        for b in range(1, max_n + 1):
            tri = compute_triangle(a, b)
            if 0.34 <= tri['alpha'] <= 0.36 or 0.34 <= tri['beta'] <= 0.36:
                pi_index, pi_chunk = map_to_pi_index(tri, pi_digits, depth)
                for p, q in twin_primes:
                    if abs(pi_index - p) < 10 or abs(pi_index - q) < 10:
                        resonant_triangles.append({
                            'triangle': tri,
                            'pi_index': pi_index,
                            'pi_chunk': pi_chunk,
                            'twin_prime': (p, q)
                        })
                        break
    return resonant_triangles

def print_resonant_triangles(resonant_triangles):
    """Print the resonant triangles."""
    print("Resonant Triangles with Alpha/Beta in [0.34, 0.36] rad:")
    for res in resonant_triangles[:10]:  # Limit to first 10 for brevity
        tri = res['triangle']
        print(f"a={tri['a']}, b={tri['b']}, alpha={tri['alpha']:.6f} rad, beta={tri['beta']:.6f} rad, "
              f"height={tri['height']:.6f}, pi_index={res['pi_index']}, "
              f"pi_chunk={res['pi_chunk']}, twin_prime={res['twin_prime']}")

# Example usage
depth = 10000  # Set the depth for pi digits
max_n = 200    # Set the maximum value for a and b

resonant_triangles = generate_resonant_triangles(depth, max_n)
print_resonant_triangles(resonant_triangles)