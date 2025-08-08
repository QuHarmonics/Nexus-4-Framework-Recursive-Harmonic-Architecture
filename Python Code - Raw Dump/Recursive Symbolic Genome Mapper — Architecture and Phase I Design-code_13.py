import itertools
import pandas as pd

def xor_collapse(seq):
    return [seq[i] ^ seq[i+1] for i in range(len(seq)-1)]

def delta_collapse(seq):
    return [abs(seq[i+1] - seq[i]) for i in range(len(seq)-1)]

def is_repeating(seq):
    return len(set(seq)) == 1

def is_valid_fixpoint(seq, delta, xor):
    if seq[0] == 0 and not is_repeating(seq):
        return False  # Reject if starts with 0 unless all zeros
    return xor == delta

# Generate all 4-digit binary sequences
sequences = list(itertools.product([0, 1], repeat=4))

# Analyze each sequence
results = []
for seq in sequences:
    seq_list = list(seq)
    delta = delta_collapse(seq_list)
    xor = xor_collapse(seq_list)
    triplet = delta  # Same as xor if they match
    valid = is_valid_fixpoint(seq_list, delta, xor)
    results.append({
        "Sequence": seq_list,
        "Δ Collapse": delta,
        "XOR Collapse": xor,
        "Triplet": triplet if valid else None,
        "XOR == Δ": xor == delta,
        "Starts With 0": seq[0] == 0,
        "Repeating Phase": is_repeating(seq),
        "All Zeros": seq == (0, 0, 0, 0),
        "Valid Fixpoint": valid
    })

# Convert to DataFrame
df = pd.DataFrame(results)

# Show all fixpoints
print("🔍 All Valid XOR–Δ Echo Fixpoints:\n")
print(df[df["Valid Fixpoint"] == True].to_string(index=False))

# Optionally export
# df.to_csv("nexus_echo_fixpoints.csv", index=False)
