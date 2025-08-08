from sympy import isprime
from mpmath import mp, nstr

mp.dps = 20000  # precision for π digits

def extract_pi_bytes(n_bytes):
    pi_str = nstr(mp.pi, n_bytes * 8 + 2)[2:]  # Remove "3."
    return [pi_str[i:i+8] for i in range(0, len(pi_str), 8)]

def find_twin_primes_in_bytes(pi_bytes):
    twin_primes = []
    for idx, byte in enumerate(pi_bytes):
        for i in range(len(byte)-1):
            try:
                a = int(byte[i])
                b = int(byte[i+1])
                if isprime(a) and isprime(b) and abs(b - a) == 2:
                    twin_primes.append((a, b, idx+1))
            except ValueError:
                continue
    return twin_primes
