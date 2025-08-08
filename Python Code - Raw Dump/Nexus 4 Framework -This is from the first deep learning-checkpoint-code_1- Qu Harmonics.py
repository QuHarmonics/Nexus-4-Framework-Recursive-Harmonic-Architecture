import random, math, pandas as pd

# --- 1.  Generate a toy DNA sequence ---------------------------------------
NUC = ['A', 'T', 'G', 'C']
random.seed(42)
seq = ''.join(random.choices(NUC, k=64))  # 64‑nt demo sequence

# --- 2.  Map nucleotides to “potential” and “actualised” energies ----------
#  (toy values just for demonstration purposes)
POTENTIAL = {'A': 2.0, 'T': 2.0, 'G': 3.0, 'C': 3.0}   # H‑bond count proxy
ACTUAL    = {'A': 1.0, 'T': 1.0, 'G': 1.5, 'C': 1.5}   # arbitrary baseline

P = [POTENTIAL[n] for n in seq]
A = [ACTUAL[n]    for n in seq]

# harmonic constant, folding factor, recursion depth per iteration
H = 0.35
F_factor = 0.5
t = 1

# compute the initial “fold contribution” vector  (Pi / Ai * e^{H F t})
import numpy as np
base_contrib = [ (p/a) * math.exp(H * F_factor * t) for p, a in zip(P, A) ]


# --- 3.  Recursive folding --------------------------------------------------
fold_iterations = []
current = base_contrib.copy()
level = 0
while len(current) > 1:
    fold_iterations.append((level, current.copy()))
    # pairwise average (simple symmetric folding)
    current = [ (current[i] + current[i+1]) / 2 for i in range(0, len(current), 2) ]
    level += 1
fold_iterations.append((level, current.copy()))  # final scalar


# --- 4.  Store in dataframe for nice viewing -------------------------------
rows = []
for lvl, vec in fold_iterations:
    rows.append({
        "Iteration": lvl,
        "Vector length": len(vec),
        "First 5 values": ', '.join(f"{v:.3f}" for v in vec[:5]) + ("…" if len(vec) > 5 else "")
    })

df_fold = pd.DataFrame(rows)

import ace_tools as tools; tools.display_dataframe_to_user("Recursive folding overview", df_fold)