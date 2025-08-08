import pandas as pd

def normalize_hash(h):
    return int.from_bytes(h, 'big') / 2**256

fold_threshold = 0.35
modulus = 64  # You can try others: 32, 128, etc.

records = []
for x in range(256):
    h = get_hash(x)
    h_norm = normalize_hash(h)
    residue = int.from_bytes(h, 'big') % modulus
    bits = [int(b) for b in format(x, '08b')]  # Bit pattern
    folded = h_norm >= fold_threshold
    records.append({
        'state': x,
        'residue': residue,
        'folded': folded,
        'h_norm': h_norm,
        'hamming_weight': sum(bits),
        'bit_pattern': bits
    })

df = pd.DataFrame(records)
display(df.head(10))  # show first few rows
