import numpy as np
import pandas as pd
import math
from ace_tools import display_dataframe_to_user

# Sample DNA sequence
seq = "ATGCGTACGTTAGC"

# Mapping nucleotides to a "potential energy" proxy
# We'll use hydrogen bond counts scaled up a bit so values aren't all ~1
bond_energy = {"A": 2.0, "T": 2.0, "G": 3.0, "C": 3.0}

# We'll treat "actualized energy" as a simple constant 1 for every base
actual_energy = 1.0

# Build a DataFrame that stores base, P_i, A_i
data = []
for i, base in enumerate(seq, 1):
    P_i = bond_energy[base]
    A_i = actual_energy
    data.append({"Pos": i, "Base": base, "P_i": P_i, "A_i": A_i})

df = pd.DataFrame(data)

# Harmonic constants
H = 0.35
F = 0.5
t = 1

# Compute contribution for each base according to folding formula
df["Fold_contrib"] = (df["P_i"] / df["A_i"]) * np.exp(H * F * t)

# Overall folded quantum state
FQ = df["Fold_contrib"].sum()

# Now demonstrate a *very* simple unfolding: spread the folded value out using phase angles
# We'll assign an equally spaced phase for each position
m = len(df)
theta = np.linspace(0, 2 * math.pi, m, endpoint=False)
df["theta"] = theta
# residual term zeta assumed 0
df["Unfold_contrib"] = (FQ / m) * np.cos(theta)

# Reconstructed (unfolded) sequence energetic profile
UQ = df["Unfold_contrib"].sum()

# Add totals row for clarity
totals = pd.DataFrame([{
    "Pos": "Σ",
    "Base": "",
    "P_i": df["P_i"].sum(),
    "A_i": df["A_i"].sum(),
    "Fold_contrib": FQ,
    "theta": "",
    "Unfold_contrib": UQ
}])

result = pd.concat([df, totals], ignore_index=True)

display_dataframe_to_user("Quantum Folding / Unfolding Example", result)