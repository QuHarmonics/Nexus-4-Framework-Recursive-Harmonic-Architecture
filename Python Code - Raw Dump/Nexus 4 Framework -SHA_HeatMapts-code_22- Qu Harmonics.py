import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# === STEP 1: Generate or Load Input ===
np.random.seed(42)
df = pd.DataFrame({'state': np.random.randint(0, 256, size=300)})

# === STEP 2: Core Field Properties ===
df['hex'] = df['state'].apply(lambda x: f"{x:02X}")
df['residue'] = df['state'] % 64
df['residue_next'] = df['residue'].shift(-1)
df['delta'] = df['residue_next'] - df['residue']
df['delta_abs'] = df['delta'].abs()
df['acceleration'] = df['delta'].diff()  # ΔΔ — acceleration of residue

# Bit breakdown
df['bits'] = df['state'].apply(lambda x: list(map(int, f"{x:08b}")))
df['bit7'] = df['bits'].apply(lambda b: b[0])
df['bit6'] = df['bits'].apply(lambda b: b[1])
df['bit5'] = df['bits'].apply(lambda b: b[2])
df['bit7_next'] = df['bit7'].shift(-1)
df['msb_flip'] = (df['bit7'] != df['bit7_next']).astype(int)

# === STEP 3: Fold Potential φ(t) Calculation ===
# Fold potential is a synthetic energy metric combining velocity, acceleration, and high bits
weights = {'delta': 0.6, 'accel': 0.3, 'bit6': 0.6, 'bit5': 0.4}
df['φ'] = (
    weights['delta'] * df['delta_abs'].fillna(0) +
    weights['accel'] * df['acceleration'].abs().fillna(0) +
    weights['bit6'] * df['bit6'] +
    weights['bit5'] * df['bit5']
)

# Normalize φ
φ_mean = df['φ'].mean()
φ_std = df['φ'].std()
df['φ_norm'] = (df['φ'] - φ_mean) / φ_std

# === STEP 4: Predict MSB Flip Based on Fold Energy Threshold ===
# Any φ_norm > 1.2 is considered a fold activation
df['predict_msb'] = (df['φ_norm'] > 1.2).astype(int)

# === STEP 5: Evaluation Metrics ===
tp = ((df['predict_msb'] == 1) & (df['msb_flip'] == 1)).sum()
fp = ((df['predict_msb'] == 1) & (df['msb_flip'] == 0)).sum()
fn = ((df['predict_msb'] == 0) & (df['msb_flip'] == 1)).sum()
precision = tp / (tp + fp) if (tp + fp) > 0 else 0
recall = tp / (tp + fn) if (tp + fn) > 0 else 0

# === STEP 6: Plot ===
plt.figure(figsize=(14, 6))
plt.plot(df['φ_norm'], label='φ(t) Normalized', color='blue')
plt.scatter(df[df['msb_flip'] == 1].index, df[df['msb_flip'] == 1]['φ_norm'], color='red', label='MSB Flip (True)', marker='x', zorder=5)
plt.scatter(df[df['predict_msb'] == 1].index, df[df['predict_msb'] == 1]['φ_norm'], color='orange', label='Prediction', marker='o', facecolors='none', edgecolors='orange', zorder=4)
plt.axhline(1.2, color='gray', linestyle='--', alpha=0.5, label='φ(t) Threshold = 1.2')
plt.xlabel('Timestep')
plt.ylabel('Fold Potential φ(t)')
plt.title('Fold Potential φ(t) vs MSB Flip Events')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

(tp, fp, fn, precision, recall)
