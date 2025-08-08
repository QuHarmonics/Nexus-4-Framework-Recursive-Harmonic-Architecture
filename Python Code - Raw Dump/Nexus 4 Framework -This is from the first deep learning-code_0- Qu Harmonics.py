import pandas as pd
import numpy as np
import math

# Sample DNA sequence
sequence = "ATGCGTAC"

# Mapping nucleotides to potential and actual energies (toy example)
potential_map = {'A': 1.0, 'T': 2.0, 'G': 3.0, 'C': 4.0}
actual_map    = {'A': 2.0, 'T': 1.0, 'G': 4.0, 'C': 3.0}

# Constants
H = 0.35          # harmonic constant
F_factor = 0.5    # folding factor
t = 1             # single recursion depth for demo

# Build a table of energies
data = []
for idx, base in enumerate(sequence, 1):
    P = potential_map[base]
    A = actual_map[base]
    ratio = P / A
    expo  = math.exp(H * F_factor * t)
    fold_component = ratio * expo
    data.append({
        "Position": idx,
        "Base": base,
        "P_i": P,
        "A_i": A,
        "P_i / A_i": ratio,
        "exp(H*F*t)": expo,
        "Fold Contribution": fold_component
    })

df = pd.DataFrame(data)
display_dataframe_to_user("Quantum Folding Components", df)

# Compute folded quantum state F(Q)
FQ = df["Fold Contribution"].sum()

# Simple unfolding demonstration using cosine phases
# Assign deterministic phases based on position for reproducibility
theta = [math.pi/4 * i for i in range(1, len(sequence)+1)]
unfold_components = [fold * math.cos(th) for fold, th in zip(df["Fold Contribution"], theta)]
UQ = sum(unfold_components)

result_df = pd.DataFrame({
    "Position": df["Position"],
    "Fold Contribution": df["Fold Contribution"],
    "Theta (rad)": theta,
    "cos(theta)": [math.cos(th) for th in theta],
    "Unfold Contribution": unfold_components
})
display_dataframe_to_user("Quantum Unfolding Components", result_df)

# Summary values
summary = pd.DataFrame({
    "Metric": ["Folded Quantum State F(Q)", "Unfolded Quantum State U(Q)"],
    "Value": [FQ, UQ]
})
display_dataframe_to_user("Quantum Folding/Unfolding Summary", summary)
