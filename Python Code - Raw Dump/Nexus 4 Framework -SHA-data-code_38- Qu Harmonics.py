import hashlib
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from scipy.stats import entropy

def get_hash(x):
    return hashlib.sha256(bytes([x])).digest()

def hash_mod_residue(x, mod):
    h = int.from_bytes(get_hash(x), 'big')
    return h % mod

# Entropy for modulus M = 2 to 256
mods = np.arange(2, 257)
entropies = []
for m in mods:
    residues = [hash_mod_residue(x, m) for x in range(256)]
    counts = np.array(list(Counter(residues).values()))
    probs = counts / counts.sum()
    entropies.append(entropy(probs, base=2))

# ΔEntropy: first-order difference
delta_entropy = np.diff(entropies)

plt.figure(figsize=(10, 5))
plt.plot(mods[1:], delta_entropy, label=r'$\Delta H(M)$', color='blue')
plt.axhline(0, color='gray', linestyle='--')
plt.title('First Derivative of Entropy vs Modulus')
plt.xlabel('Modulus')
plt.ylabel(r'$\Delta H$')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
