# ------------------------------------------------------
# 0.  Utilities
# ------------------------------------------------------
from pathlib import Path
from random import randrange, shuffle
from collections import Counter
import math

def is_probable_prime(n, rounds=8):
    """Deterministic Miller–Rabin for n < 3,474,749,660,383 (rounds=8)."""
    if n < 2: return False
    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    if n in small: return True
    if any(n % p == 0 for p in small): return False
    # write n−1 = d·2^s
    d, s = n - 1, 0
    while d & 1 == 0:
        d >>= 1; s += 1
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

# ------------------------------------------------------
# 1.  Load π digits (plain text file with no whitespace)
# ------------------------------------------------------
PI_PATH = Path("pi_1M.txt")   # put 1 000 002 chars: "3." + 1 000 000 digits
digits = PI_PATH.read_text().strip().lstrip("3.").rstrip()
N = len(digits)
print(f"Loaded {N} π digits.")

# ------------------------------------------------------
# 2.  Scan Δ = 2 events
# ------------------------------------------------------
hits = []  # (index k, d_k, left_digit, right_digit, concat_val, flags)
for k in range(N - 1):
    left  = int(digits[k])
    right = int(digits[k + 1])
    d_k   = abs(right - left)
    if d_k == 2:
        concat_val = 10 * left + right
        flags = (
            is_probable_prime(left),
            is_probable_prime(right),
            is_probable_prime(concat_val),
            is_probable_prime(k),
        )
        hits.append((k, left, right, concat_val, flags))

print(f"Total Δ=2 events: {len(hits)}")

# ------------------------------------------------------
# 3.  Statistics helper
# ------------------------------------------------------
def summarise(hit_records):
    cnt = Counter()
    for _, l, r, c, (pL, pR, pC, pK) in hit_records:
        cnt["left"]  += pL
        cnt["right"] += pR
        cnt["concat"]+= pC
        cnt["index"] += pK
    total = len(hit_records)
    return {k: f"{v}/{total}  ({v/total:.3%})" for k, v in cnt.items()}

print("π statistics:", summarise(hits))

# ------------------------------------------------------
# 4.  Control sequences
# ------------------------------------------------------
def shuffled_hits(seq):
    tmp = list(seq)
    shuffle(tmp)
    return tmp

def control_run(seq_digits):
    H = []
    for k in range(len(seq_digits) - 1):
        l, r = int(seq_digits[k]), int(seq_digits[k+1])
        if abs(l - r) == 2:
            c = 10*l + r
            H.append((k, l, r, c, (
                is_probable_prime(l),
                is_probable_prime(r),
                is_probable_prime(c),
                is_probable_prime(k),
            )))
    return summarise(H)

# digits of e
E_PATH = Path("e_1M.txt")
digits_e = E_PATH.read_text().strip().lstrip("2.").rstrip()
print(" e statistics:", control_run(digits_e))

# shuffled π
pi_shuffled = shuffled_hits(digits)
print("π shuffled   :", control_run(pi_shuffled))

# random digits
import secrets
rand_digits = ''.join(str(secrets.randbelow(10)) for _ in range(N))
print(" random       :", control_run(rand_digits))
