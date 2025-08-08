import hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import entropy

# === PARAMETERS ===
MODULUS = 64
DELTA_THRESHOLD = 40
FOLD_THRESHOLD = 0.35
NUM_STATES = 256

# === UTILITIES ===
def sha256_bytes(x):
    return hashlib.sha256(bytes([x % 256])).digest()

def normalize_hash(h_bytes):
    return int.from_bytes(h_bytes, 'big') / (2 ** 256)

def bit_array(x):
    return np.array([int(b) for b in format(x, '08b')])

# === MAIN LOOP ===
data = []
for x in range(NUM_STATES):
    h = sha256_bytes(x)
    h_norm = normalize_hash(h)
    residue = int.from_bytes(h, 'big') % MODULUS
    folded = h_norm >= FOLD_THRESHOLD
    bits = bit_array(x)
    hw = np.sum(bits)
    data.append(dict(state=x, residue=residue, h_norm=h_norm, folded=folded,
                     hamming_weight=hw, bits=bits))

df = pd.DataFrame(data)
bit_matrix = np.stack(df['bits'].values)

# === TRANSITION METRICS ===
df_trans = df.iloc[:-1].copy()
df_trans['residue_next'] = df['residue'].iloc[1:].values
df_trans['x+1'] = df['state'].iloc[1:].values
df_trans['delta_residue'] = (df_trans['residue_next'] - df_trans['residue']) % MODULUS
df_trans['bit_flip_count'] = np.sum(np.abs(np.diff(bit_matrix.astype(int), axis=0)), axis=1)

# === FOLDED VS UNFOLDED BIT DRIVER ANALYSIS ===
folded_bits = bit_matrix[df['folded']]
unfolded_bits = bit_matrix[~df['folded']]
bit_delta = folded_bits.mean(axis=0) - unfolded_bits.mean(axis=0)

plt.figure(figsize=(8, 4))
plt.bar(range(8), bit_delta, color='crimson', alpha=0.8)
plt.xlabel("Bit Position (7 = MSB)")
plt.ylabel("Δ Activation (Folded - Unfolded)")
plt.title("Bit Activation Difference: Folded vs Unfolded States")
plt.grid(True)
plt.tight_layout()
plt.show()

# === ENTROPY OF HASHES ACROSS MODULI ===
entropies = []
moduli = range(2, 257)
for m in moduli:
    residues = [(int.from_bytes(sha256_bytes(x), 'big') % m) for x in range(NUM_STATES)]
    counts = np.bincount(residues, minlength=m)
    p = counts / counts.sum()
    entropies.append(entropy(p, base=2))

# === ENTROPY DERIVATIVE
entropy_deriv = np.diff(entropies)
plt.figure(figsize=(8, 4))
plt.plot(moduli[1:], entropy_deriv, label=r'$\Delta H(M)$', color='blue')
plt.axhline(0, linestyle='--', color='black', alpha=0.5)
plt.title("First Derivative of Entropy vs Modulus")
plt.xlabel("Modulus")
plt.ylabel(r'$\Delta H$')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# === HIGH DELTA_RESIDUE EVENTS — BIT DRIVER ANALYSIS
high_jump_rows = df_trans[df_trans['delta_residue'] > DELTA_THRESHOLD]['x'].tolist()
bit_flips = []
for x in high_jump_rows:
    b1 = format(x, '08b')
    b2 = format(x+1, '08b')
    flips = [int(c1 != c2) for c1, c2 in zip(b1, b2)]
    bit_flips.append(flips)

bit_flip_array = np.array(bit_flips)
bit_contrib = bit_flip_array.mean(axis=0)

plt.figure(figsize=(8, 4))
plt.bar(range(8), bit_contrib, color='navy', alpha=0.8)
plt.xlabel("Bit Position")
plt.ylabel("Flip Frequency in High ΔResidue Events")
plt.title("Driver Bits for Large Residue Transitions")
plt.grid(True)
plt.tight_layout()
plt.show()

# === OPTIONAL: FLIP COUNT vs ΔResidue
plt.figure(figsize=(8, 4))
sns.scatterplot(data=df_trans, x='bit_flip_count', y='delta_residue', alpha=0.6)
plt.title("Residue Jump vs Bit Transitions")
plt.xlabel("Bit Flip Count (x → x+1)")
plt.ylabel(r'$\Delta$Residue (mod {})'.format(MODULUS))
plt.tight_layout()
plt.grid(True)
plt.show()
