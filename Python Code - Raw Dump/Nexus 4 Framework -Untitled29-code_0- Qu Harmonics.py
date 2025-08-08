import math
import hashlib
from mpmath import mp
import matplotlib.pyplot as plt

# Set precision for π digits
mp.dps = 10000
pi_str = str(mp.pi())[2:10002]
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

# KHRC correction pass (not applied to system, just logged)
def khrc(H_val, k=1.0, N=0.01):
    return H_val / (1 + k * abs(N))

# Main process
resonant_triangles = []
H_trajectory = []

for a in range(1, 200):
    for b in range(1, 200):
        tri = compute_triangle(a, b)
        if 0.34 <= tri['alpha'] <= 0.36 or 0.34 <= tri['beta'] <= 0.36:
            pi_index, pi_chunk = map_to_pi_index(tri)
            for p, q in twin_primes:
                if abs(pi_index - p) < TWIN_WINDOW or abs(pi_index - q) < TWIN_WINDOW:
                    H_val = compute_H(pi_chunk)
                    H_trajectory.append(H_val)
                    res = {
                        'triangle': tri,
                        'pi_index': pi_index,
                        'pi_chunk': pi_chunk,
                        'twin_prime': (p, q),
                        'H': H_val,
                        'KHRC': khrc(H_val)
                    }
                    resonant_triangles.append(res)
                    break

# Print summary
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
