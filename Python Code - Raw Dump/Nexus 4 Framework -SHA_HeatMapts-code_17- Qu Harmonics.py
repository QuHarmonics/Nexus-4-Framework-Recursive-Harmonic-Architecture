import pandas as pd

# Assume df contains columns: ['state', 'residue', 'delta_residue', 'bit_pattern']

# Compute MSB flip
df['bit7'] = df['bit_pattern'].apply(lambda bits: bits[0])
df['bit7_next'] = df['bit7'].shift(-1)
df['msb_flip'] = (df['bit7'] != df['bit7_next']).astype(int)

# Compute delta residue
df['residue_next'] = df['residue'].shift(-1)
df['delta_residue'] = (df['residue_next'] - df['residue']).abs()

# Thresholding
THRESH = 30  # Empirically tweak
df['flip_predicted'] = (df['delta_residue'] >= THRESH).astype(int)

# Evaluate predictive power
tp = ((df['flip_predicted'] == 1) & (df['msb_flip'] == 1)).sum()
fp = ((df['flip_predicted'] == 1) & (df['msb_flip'] == 0)).sum()
fn = ((df['flip_predicted'] == 0) & (df['msb_flip'] == 1)).sum()
precision = tp / (tp + fp)
recall = tp / (tp + fn)

print(f"Precision: {precision:.2f}, Recall: {recall:.2f}")

