#!/usr/bin/env python3
# pi_delta2_twin_scan.py
#
# ------------------------------------------------------------
# (1)  Pure-Python Chudnovsky π generator (decimal module)
# ------------------------------------------------------------
import decimal
from pathlib import Path
from random import randrange, shuffle
from collections import Counter
import secrets
import time

def pi_digits(n: int) -> str:
    """
    Return *exactly* n decimal digits of π (after the decimal point),
    using the Chudnovsky series and Python's `decimal` module.
    Runtime ~45 s for 1 M digits on a 3 GHz laptop (pure Python).
    """
    decimal.getcontext().prec = n + 20          # guard digits
    C = 426880 * decimal.Decimal(10005).sqrt()
    K, M, L, X, S = 6, 1, 13591409, 1, 13591409
    for k in range(1, int(n / 14) + 2):         # each term adds ~14 digits
        M = (M * (decimal.Decimal(6*k-5) *
                  (2*k-1) * (6*k-1))) / (k**3 * 26680 * 640320**3)
        L += 545140134
        X *= -262537412640768000
        S += decimal.Decimal(M * L) / X
    pi_val = C / S
    digits = str(pi_val).replace('.', '')[1:n+1]  # strip "3."
    return digits

# ------------------------------------------------------------
# (2)  Deterministic Miller–Rabin (safe < 3.4e12 with 8 bases)
# ------------------------------------------------------------
def is_probable_prime(n: int, rounds: int = 8) -> bool:
    if n < 2:
        return False
    small_primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29)
    if n in small_primes:
        return True
    if any(n % p == 0 for p in small_primes):
        return False

    # n − 1 = d·2^s
    d, s = n - 1, 0
    while d & 1 == 0:
        d >>= 1
        s += 1
    for _ in range(rounds):
        a = randrange(2, n - 1)
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

# ------------------------------------------------------------
# (3)  Δ = 2 scanner
# ------------------------------------------------------------
def delta2_scan(digit_str: str):
    """
    Scan adjacent digit pairs, logging every |Δ| == 2 event.
    Returns list of tuples:
        (index, left_digit, right_digit, concat_val, flags_tuple)
    where flags_tuple = (prime_left, prime_right, prime_concat, prime_index)
    """
    hits = []
    N = len(digit_str)
    for k in range(N - 1):
        l = int(digit_str[k])
        r = int(digit_str[k + 1])
        if abs(l - r) == 2:
            concat_val = 10 * l + r
            flags = (
                is_probable_prime(l),
                is_probable_prime(r),
                is_probable_prime(concat_val),
                is_probable_prime(k),
            )
            hits.append((k, l, r, concat_val, flags))
    return hits

def summarise(hits):
    cnt = Counter()
    for _, l, r, c, (pL, pR, pC, pK) in hits:
        cnt["left"]   += pL
        cnt["right"]  += pR
        cnt["concat"] += pC
        cnt["index"]  += pK
    total = len(hits)
    return {k: f"{v}/{total}  ({v/total:.3%})" for k, v in cnt.items()}

# ------------------------------------------------------------
# (4)  Control helpers
# ------------------------------------------------------------
def control_stats(digit_str: str):
    return summarise(delta2_scan(digit_str))

def shuffled_sequence(seq: str):
    s = list(seq)
    shuffle(s)
    return ''.join(s)

def random_digits(length: int):
    return ''.join(str(secrets.randbelow(10)) for _ in range(length))

# ------------------------------------------------------------
# (5)  Main experiment driver
# ------------------------------------------------------------
def main(num_digits: int = 100_000):
    print(f"Generating {num_digits:,} π digits …")
    t0 = time.time()
    pi_str = pi_digits(num_digits)
    print(f" π ready in {time.time()-t0:.2f} s")

    print("Scanning Δ = 2 events in π …")
    hits_pi = delta2_scan(pi_str)
    print(f"Total events: {len(hits_pi)}")
    print("π stats      :", summarise(hits_pi))

    # ----- controls -----
    print("\nGenerating e digits for control …")
    t0 = time.time()
    e_str = e_digits(num_digits)
    print(f" e ready in {time.time()-t0:.2f} s")
    print("Scanning e …")
    print("e stats       :", control_stats(e_str))

    print("Scanning shuffled π …")
    shuf_stats = control_stats(shuffled_sequence(pi_str))
    print("π shuffled    :", shuf_stats)

    print("Scanning random digits …")
    rand_stats = control_stats(random_digits(num_digits))
    print("random        :", rand_stats)

# ------------------------------------------------------------
# (6)  Quick e-digit generator (Bailey–Borwein–Plouffe-like)
# ------------------------------------------------------------
def e_digits(n: int) -> str:
    """
    Generate n digits of e (naïve series, adequate for ≤1e6 digits; slow but OK).
    """
    decimal.getcontext().prec = n + 10
    e_val = decimal.Decimal(0)
    factorial = 1
    i = 0
    while len(str(e_val).replace('.', '')) < n + 1:
        e_val += decimal.Decimal(1) / factorial
        i += 1
        factorial *= i
    digits = str(e_val).replace('.', '')[1:n+1]   # strip "2."
    return digits

# ------------------------------------------------------------
# Run module as script
# ------------------------------------------------------------
if __name__ == "__main__":
    main(num_digits=100_000)   # adjust to e.g. 1_000_000 for deeper test
