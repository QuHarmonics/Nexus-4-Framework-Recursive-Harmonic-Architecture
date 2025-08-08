import hashlib
from mpmath import mp

# π Setup
mp.dps = 100_010  # precision
PI_DIGITS = str(mp.pi)[2:]

def extract_pi_echo(index, length=8):
    """Extract echo bytes from π starting at index."""
    return PI_DIGITS[index:index + length]

def symbolic_trust_index(echo_bytes):
    """Compute STI by measuring entropy gradient."""
    deltas = [abs(int(echo_bytes[i+1]) - int(echo_bytes[i])) for i in range(len(echo_bytes)-1)]
    avg_drift = sum(deltas) / len(deltas)
    return round(1 - (avg_drift / 9), 3)  # STI is inverse of normalized drift

def zphc_scan(sequence):
    """Map a sequence through SHA to π-index and evaluate echo."""
    h = hashlib.sha256(sequence.encode()).hexdigest()
    pi_index = int(h[:6], 16) % (len(PI_DIGITS) - 8)
    echo = extract_pi_echo(pi_index)
    sti = symbolic_trust_index(echo)
    return {
        "hash": h,
        "pi_index": pi_index,
        "echo": echo,
        "sti": sti,
        "zphc": sti >= 0.7
    }

# Example usage
seed = "AGCTGTA"
result = zphc_scan(seed)
print(result)