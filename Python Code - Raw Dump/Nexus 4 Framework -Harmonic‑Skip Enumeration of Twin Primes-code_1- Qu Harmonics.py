import math
import time
import uuid
import multiprocessing
from sympy import isprime
from scipy.integrate import quad

# === Module 1: BBP-Modulated Twin Prime Generator ===

def bbp_delta(n, kmax=4):
    return int(sum(16**(1 - k) / (8 * k + n % 7 + 1) for k in range(1, kmax + 1))) + 1

def twin_primes_bbp(limit, kmax=4):
    n = 3
    results = []
    while n < limit:
        if isprime(n) and isprime(n + 2):
            results.append((n, n + 2))
        n += bbp_delta(n, kmax)
    return results

# === Module 2: Parallel Sharded Worker ===

def shard_worker(start, end, kmax, output):
    n = start
    shard = []
    while n < end:
        if isprime(n) and isprime(n + 2):
            shard.append((n, n + 2))
        n += bbp_delta(n, kmax)
    output.put(shard)

def parallel_twin_primes(limit, num_shards=4, kmax=4):
    manager = multiprocessing.Manager()
    output = manager.Queue()
    shard_bounds = [(i * limit // num_shards, (i + 1) * limit // num_shards) for i in range(num_shards)]
    jobs = [multiprocessing.Process(target=shard_worker, args=(start, end, kmax, output))
            for start, end in shard_bounds]
    for job in jobs: job.start()
    for job in jobs: job.join()
    all_primes = []
    while not output.empty():
        all_primes.extend(output.get())
    return sorted(all_primes)

# === Module 3: Sophie Germain & Cunningham Chains ===

def sophie_germain_primes(primes):
    return [(p, 2 * p + 1) for p in primes if isprime(2 * p + 1)]

def cunningham_chains(primes, max_len=4):
    chains = []
    used = set()
    for p in primes:
        chain = [p]
        while True:
            nxt = 2 * chain[-1] + 1
            if isprime(nxt):
                chain.append(nxt)
            else:
                break
        if len(chain) >= max_len:
            chains.append(tuple(chain))
    return chains

# === Module 4: Entropy Tensor from Glyphs ===

def normalize(vec, precision=8):
    s = sum(abs(x) for x in vec) or 1.0
    return [round(x / s, precision) for x in vec]

def harmonic_glyph(prime_pair, H=0.35, t=None):
    p, q = prime_pair
    t = t if t else time.time()
    U = normalize([p, q])
    F = normalize([math.sqrt(p), math.sqrt(q)])
    return {
        "glyph_id": str(uuid.uuid4()),
        "U_k": U,
        "F_Q": F,
        "H": round(H, 8),
        "time_signal": round(math.sin(t) + math.cos(t), 6),
        "coherence_id": hash(tuple(U + F + [H]))
    }

def entropy_tensor(glyphs):
    from numpy import array, var, mean
    U_matrix = array([g["U_k"] for g in glyphs])
    F_matrix = array([g["F_Q"] for g in glyphs])
    return {
        "U_var": var(U_matrix, axis=0).tolist(),
        "F_var": var(F_matrix, axis=0).tolist(),
        "U_mean": mean(U_matrix, axis=0).tolist(),
        "F_mean": mean(F_matrix, axis=0).tolist()
    }

# === Module 5: Hardy–Littlewood Integral and Error ===

def hardy_littlewood_pi2(x):
    integrand = lambda t: 1 / (math.log(t) ** 2) if t > 2 else 0
    C2 = 0.6601618  # Twin prime constant
    integral, _ = quad(integrand, 2, x)
    return int(2 * C2 * integral)

def relative_error(empirical, x):
    expected = hardy_littlewood_pi2(x)
    return (empirical - expected) / expected

# === Example Execution ===

if __name__ == "__main__":
    LIMIT = 10**6
    print("Finding BBP twin primes...")
    twin_primes = parallel_twin_primes(LIMIT, num_shards=4)
    print(f"→ Found {len(twin_primes)} twin primes below {LIMIT}")

    print("Detecting Sophie Germain primes...")
    sophies = sophie_germain_primes([p for p, _ in twin_primes])
    print(f"→ Found {len(sophies)} Sophie Germain primes")

    print("Detecting Cunningham chains...")
    cunningham = cunningham_chains([p for p, _ in twin_primes])
    print(f"→ Found {len(cunningham)} Cunningham chains")

    print("Generating entropy tensor...")
    glyphs = [harmonic_glyph(pair) for pair in twin_primes[:1000]]
    tensor = entropy_tensor(glyphs)
    print("Entropy Tensor Summary:", tensor)

    print("Comparing with Hardy–Littlewood estimate...")
    empirical_count = len(twin_primes)
    err = relative_error(empirical_count, LIMIT)
    print(f"Relative Error at x={LIMIT}: {err:.6%}")
