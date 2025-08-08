import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# === Generate or load your input data ===
# Ensure you have a DataFrame `df` with a 'state' column of integer (0–255)

# Example mock setup:
np.random.seed(42)
df = pd.DataFrame({'state': np.random.randint(0, 256, size=300)})

# === Compute Residue (mod 64) and Delta ===
df['residue'] = df['state'] % 64
df['residue_next'] = df['residue'].shift(-1)
df['delta_residue'] = (df['residue_next'] - df['residue']).abs()

# === Compute Bit Pattern and MSB ===
df['bit_pattern'] = df['state'].apply(lambda x: list(map(int, f"{x:08b}")))
df['bit7'] = df['bit_pattern'].apply(lambda bits: bits[0])  # MSB (bit 7)
df['bit7_next'] = df['bit7'].shift(-1)
df['msb_flip'] = (df['bit7'] != df['bit7_next']).astype(int)

# === Compute Hamming distance over previous 3 steps ===
def hamming(a, b):
    return bin(a ^ b).count('1')

df['state_prev1'] = df['state'].shift(1)
df['state_prev2'] = df['state'].shift(2)
df['state_prev3'] = df['state'].shift(3)

df['hamm_1'] = df.apply(lambda r: hamming(int(r['state']), int(r['state_prev1'])) if pd.notnull(r['state_prev1']) else 0, axis=1)
df['hamm_2'] = df.apply(lambda r: hamming(int(r['state_prev1']), int(r['state_prev2'])) if pd.notnull(r['state_prev2']) else 0, axis=1)
df['hamm_3'] = df.apply(lambda r: hamming(int(r['state_prev2']), int(r['state_prev3'])) if pd.notnull(r['state_prev3']) else 0, axis=1)

# === Heuristic to Predict MSB Flip ===
df['predict_msb'] = (
    (df['hamm_1'] >= 4) &
    (df['hamm_2'] >= 3) &
    (df['bit7'] == 0)
).astype(int)

# === Evaluation ===
tp = ((df['predict_msb'] == 1) & (df['msb_flip'] == 1)).sum()
fp = ((df['predict_msb'] == 1) & (df['msb_flip'] == 0)).sum()
fn = ((df['predict_msb'] == 0) & (df['msb_flip'] == 1)).sum()
precision = tp / (tp + fp) if (tp + fp) > 0 else 0
recall = tp / (tp + fn) if (tp + fn) > 0 else 0

print(f"=== Predicting MSB Flip with Hamming Cascade ===")
print(f"Precision: {precision:.3f}")
print(f"Recall:    {recall:.3f}")
print(f"True Positives: {tp}, False Positives: {fp}, False Negatives: {fn}\n")

# === Plot ===
plt.figure(figsize=(12, 6))
plt.plot(df['delta_residue'], label='ΔResidue', color='blue')
plt.scatter(df[df['msb_flip'] == 1].index, df[df['msb_flip'] == 1]['delta_residue'], color='red', label='MSB Flip', marker='x', zorder=5)
plt.scatter(df[df['predict_msb'] == 1].index, df[df['predict_msb'] == 1]['delta_residue'], color='orange', label='Prediction', marker='o', facecolors='none', edgecolors='orange', zorder=4)
plt.axhline(30, color='gray', linestyle='--', alpha=0.5, label='ΔResidue Threshold 30')
plt.xlabel('Timestep')
plt.ylabel('ΔResidue')
plt.title('ΔResidue vs. MSB Flip and Predicted Events')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
