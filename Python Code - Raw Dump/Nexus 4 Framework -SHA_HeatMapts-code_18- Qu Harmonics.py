import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# === (Assuming df is already loaded and contains at least 'state', 'residue', and 'h_norm') ===

# Compute next residue and delta
df['residue_next'] = df['residue'].shift(-1)
df['delta_residue'] = (df['residue_next'] - df['residue']).abs()

# Compute 8-bit binary pattern from 'state'
df['bit_pattern'] = df['state'].apply(lambda x: [int(b) for b in f"{x:08b}"])

# Extract bit 7 (MSB)
df['bit7'] = df['bit_pattern'].apply(lambda bits: bits[0])
df['bit7_next'] = df['bit7'].shift(-1)
df['msb_flip'] = (df['bit7'] != df['bit7_next']).astype(int)

# Fold condition: normalized SHA hash above threshold
df['fold'] = df['h_norm'] >= 0.35

# Define residue jump threshold for prediction
DELTA_THRESHOLD = 30
df['flip_predicted'] = (df['delta_residue'] >= DELTA_THRESHOLD).astype(int)

# === Metrics: Precision and Recall
tp = ((df['flip_predicted'] == 1) & (df['msb_flip'] == 1)).sum()
fp = ((df['flip_predicted'] == 1) & (df['msb_flip'] == 0)).sum()
fn = ((df['flip_predicted'] == 0) & (df['msb_flip'] == 1)).sum()

precision = tp / (tp + fp) if (tp + fp) else 0
recall = tp / (tp + fn) if (tp + fn) else 0

print(f"\n🧪 Predicting MSB Flip via ΔResidue ≥ {DELTA_THRESHOLD}")
print(f"Precision: {precision:.3f}")
print(f"Recall:    {recall:.3f}")
print(f"True Positives: {tp}, False Positives: {fp}, False Negatives: {fn}\n")

# === Plotting
plt.figure(figsize=(12, 5))
plt.plot(df.index, df['delta_residue'], label='ΔResidue', color='blue')
plt.scatter(df[df['msb_flip'] == 1].index, df[df['msb_flip'] == 1]['delta_residue'],
            color='red', label='MSB Flip', marker='x')
plt.axhline(DELTA_THRESHOLD, color='gray', linestyle='--', label=f'Threshold {DELTA_THRESHOLD}')
plt.title("Residue Delta vs. MSB Flip Events")
plt.xlabel("Timestep")
plt.ylabel("ΔResidue")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
