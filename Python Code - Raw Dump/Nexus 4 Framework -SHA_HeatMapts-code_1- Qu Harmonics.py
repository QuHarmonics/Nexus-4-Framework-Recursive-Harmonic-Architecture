import hashlib

# 1) Generate first 64 primes (simple sieve or trial division)
def first_n_primes(n):
    primes = []
    candidate = 2
    while len(primes) < n:
        is_p = all(candidate % p for p in primes if p*p <= candidate)
        if is_p:
            primes.append(candidate)
        candidate += 1
    return primes

primes64 = first_n_primes(64)

# 2) Find twin-prime indices
twin_indices = [i for i in range(63) if primes64[i+1] - primes64[i] == 2]
# we'll also include both members of each twin pair:
twin_indices = sorted(set(twin_indices + [i+1 for i in twin_indices]))

print("Twin-prime positions (0-based j):", twin_indices)
print("Twin-primes themselves: ", [primes64[j] for j in twin_indices])

# 3) Build the 64 round constants Kj (from FIPS 180-4)
#    Here we use hashlib to re-derive them; in practice you'd hard-code them.
def round_constants():
    # fractional cube roots of primes: fract(cuberoot(prime_j)) * 2**32
    Ks = []
    for p in primes64:
        frac = (p ** (1/3)) % 1
        Kj = int(frac * (1<<32)) & 0xFFFFFFFF
        Ks.append(Kj)
    return Ks

Ks = round_constants()
K_mod256 = [k & 0xFF for k in Ks]

# 4) Your 13-cycle states:
cycle = [21, 8, 154, 65, 253, 155, 94, 203, 29, 205, 209, 153, 88]

# 5) Which states land on a twin-prime-derived constant?
hits = {}
for s in cycle:
    hits[s] = [j for j in twin_indices if (s % 256) == K_mod256[j]]

print("\nState → matching twin-prime constant indices:")
for s in cycle:
    print(f" {s:3} → {hits[s]}")

# 6) Summarize
matched = [s for s in cycle if hits[s]]
print("\nStates that hit a twin-prime constant (mod256):", matched)
