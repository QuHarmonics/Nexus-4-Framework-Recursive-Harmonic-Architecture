import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# === Generate or load input data ===
np.random.seed(42)
df = pd.DataFrame({'state': np.random.randint(0, 256, size=300)})

# === Compute Residue and Delta ===
df['residue'] = df['state'] % 64
df['residue_next'] = df['residue'].shift(-1)
df['delta_residue'] = (df['residue_next'] - df['residue']).abs()

# === Bit Pattern and MSB extraction ===
df['bit_pattern'] = df['state'].apply(lambda x: list(map(int, f"{x:08b}")))
df['bit7'] = df['bit_pattern'].apply(lambda bits: bits[0])  # MSB
df['bit6'] = df['bit_pattern'].apply(lambda bits: bits[1])  # Second MSB
df['bit7_next'] = df['bit7'].shift(-1)
df['msb_flip'] = (df['bit7'] != df['bit7_next']).astype(int)

# === Compute Hamming distances ===
def hamming(a, b):
    return bin(a ^ b).count('1')

df['state_prev1'] = df['state'].shift(1)
df['state_prev2'] = df['state'].shift(2)
df['state_prev3'] = df['state'].shift(3)

df['hamm_1'] = df.apply(lambda r: hamming(int(r['state']), int(r['state_prev1'])) if pd.notnull(r['state_prev1']) else 0, axis=1)
df['hamm_2'] = df.apply(lambda r: hamming(int(r['state_prev1']), int(r['state_prev2'])) if pd.notnull(r['state_prev2']) else 0, axis=1)
df['hamm_3'] = df.apply(lambda r: hamming(int(r['state_prev2']), int(r['state_prev3'])) if pd.notnull(r['state_prev3']) else 0, axis=1)

# === Fold Potential Function φ(t) ===
# Adjustable weights
w1, w2, w3, w4 = 0.6, 0.4, 0.15, 1.2
df['phi'] = (
    w1 * df['hamm_1'] +
    w2 * df['hamm_2'] +
    w3 * df['delta_residue'].fillna(0) +
    w4 * df['bit6']
)

# === Predict MSB Flip from φ(t)
# Fold threshold is heuristic; tune as needed
fold_threshold = df['phi'].quantile(0.85)
df['predict_msb'] = (df['phi'] > fold_threshold).astype(int)

# === Evaluation ===
tp = ((df['predict_msb'] == 1) & (df['msb_flip'] == 1)).sum()
fp = ((df['predict_msb'] == 1) & (df['msb_flip'] == 0)).sum()
fn = ((df['predict_msb'] == 0) & (df['msb_flip'] == 1)).sum()
precision = tp / (tp + fp) if (tp + fp) > 0 else 0
recall = tp / (tp + fn) if (tp + fn) > 0 else 0

print("=== MSB Flip Prediction via φ(t) Fold Potential ===")
print(f"True Positives : {tp}")
print(f"False Positives: {fp}")
print(f"False Negatives: {fn}")
print(f"Precision      : {precision:.2f}")
print(f"Recall         : {recall:.2f}\n")

# === Plotting ΔResidue and Prediction Overlays ===
plt.figure(figsize=(12, 6))
plt.plot(df['delta_residue'], label='ΔResidue', color='steelblue')
plt.scatter(df[df['msb_flip'] == 1].index, df[df['msb_flip'] == 1]['delta_residue'], 
            color='red', label='MSB Flip (Ground Truth)', marker='x', zorder=5)
plt.scatter(df[df['predict_msb'] == 1].index, df[df['predict_msb'] == 1]['delta_residue'], 
            facecolors='none', edgecolors='orange', label='Predicted Flip', marker='o', zorder=4)
plt.axhline(30, color='gray', linestyle='--', alpha=0.5, label='Reference Threshold')
plt.xlabel('Timestep')
plt.ylabel('ΔResidue')
plt.title('Fold-State Tracking via ΔResidue and φ(t) Prediction')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
