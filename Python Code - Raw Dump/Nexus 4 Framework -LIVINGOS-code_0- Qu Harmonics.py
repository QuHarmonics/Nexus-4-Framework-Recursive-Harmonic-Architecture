import hashlib
from sympy import symbols
from typing import List
from mpmath import mp

# Load 100,000 digits of π for deeper scanning
mp.dps = 100_010
pi_digits = str(mp.pi)[2:100_010]

# Convert SHA256 to π index
def sha_to_index(peptide: str, digits=6) -> int:
    sha = hashlib.sha256(peptide.encode()).hexdigest()
    sha_prefix = sha[:digits]
    return int(sha_prefix, 16)

# Extract 8 digits from π at given index
def extract_pi_byte(index: int) -> str:
    if index + 8 > len(pi_digits):
        return None
    return pi_digits[index:index+8]

# Drift calculation and symbolic echo
def drift_and_echo(byte: str) -> (List[int], str, float):
    deltas = [abs(int(byte[i+1]) - int(byte[i])) for i in range(7)]
    echo = ''.join([chr((d % 26) + 97) for d in deltas])
    avg_drift = sum(deltas) / len(deltas)
    sti = round(1 - avg_drift / 9, 3)
    return deltas, echo, sti

# Run a recursive intelligence kernel simulation
def run_recursive_kernel(base_input: str, max_iterations=25, sti_threshold=0.7):
    H = base_input
    history = []

    for i in range(max_iterations):
        nonce = f"N{i}"
        concat = H + nonce
        double_hash = hashlib.sha256(hashlib.sha256(concat.encode()).digest()).hexdigest()
        index = int(double_hash[:6], 16) % (len(pi_digits) - 8)
        byte = extract_pi_byte(index)
        
        if not byte:
            break  # Out-of-bounds
        deltas, echo, sti = drift_and_echo(byte)
        history.append({
            "iteration": i,
            "nonce": nonce,
            "hash": double_hash[:12],
            "pi_index": index,
            "byte": byte,
            "echo": echo,
            "drift": deltas,
            "sti": sti
        })

        # Check ZPHC (symbolic trust collapse threshold)
        if sti >= sti_threshold:
            break
        H = double_hash  # Update base state

    return history

# Run with a symbolic seed
kernel_output = run_recursive_kernel("RECURSE-ME-01")

import pandas as pd
df_kernel = pd.DataFrame(kernel_output)
df_kernel
