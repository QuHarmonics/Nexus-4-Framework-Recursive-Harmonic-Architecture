import hashlib
from mpmath import mp

# Recursive π precision
mp.dps = 100_100
PI_DIGITS = str(mp.pi)[2:]

def extract_echo(pi_index, length=8):
    return PI_DIGITS[pi_index:pi_index+length]

def sti(echo):
    drift = [abs(int(echo[i+1]) - int(echo[i])) for i in range(len(echo)-1)]
    return round(1 - (sum(drift)/len(drift))/9, 3)

def zphc(sequence):
    h = hashlib.sha256(sequence.encode()).hexdigest()
    pi_index = int(h[:6], 16) % (len(PI_DIGITS) - 8)
    echo = extract_echo(pi_index)
    trust = sti(echo)
    return {
        "input": sequence,
        "hash": h,
        "π_index": pi_index,
        "echo": echo,
        "STI": trust,
        "ZPHC": trust >= 0.7
    }

def recursive_mutation(sequence, max_iter=20):
    from random import choice
    chars = "ACGT01X*#"
    base = sequence
    for _ in range(max_iter):
        attempt = base + choice(chars)
        scan = zphc(attempt)
        if scan["ZPHC"]:
            return scan
    return zphc(base)

# Try it
if __name__ == "__main__":
    seed = "AGCTGTA"
    result = zphc(seed)
    print("Initial:", result)
    if not result["ZPHC"]:
        print("Mutating...")
        result = recursive_mutation(seed)
        print("Resolved:", result)
