import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ─────────────── Physical / model constants ───────────────
m       = 1.0                       # mass (arbitrary units)
v       = 10.0                      # linear velocity (same units)
c       = 299_792_458.0             # speed of light (m/s)
γ       = 1.0 / np.sqrt(1 - (v**2)/(c**2))  # Lorentz factor
n       = 1                         # quantum level
ℏ       = 1.054571817e-34           # reduced Planck’s constant (J·s)
L       = 1.0                       # length scale for quantum well
r_chaos = 3.9                       # logistic‑map parameter
x0      = 0.5                       # initial x for chaos map

# ─────────────── Parameter ranges ───────────────
omega_values = np.linspace(0, 10, 20)    # angular velocity sweep
r_values     = np.linspace(0.1, 5, 20)   # radius sweep

# ─────────────── Storage for results ───────────────
ke_classical_results      = []
ke_relativistic_results   = []
e_quantum_results         = []
chaos_results             = []

# ─────────────── Perform the sweeps ───────────────
for ω in omega_values:
    for r in r_values:
        # Classical kinetic energy
        ke_cl = 0.5 * m * (v**2 + (ω**2)*(r**2))
        ke_classical_results.append([ω, r, ke_cl])
        
        # Relativistic kinetic energy (approximate form here)
        ke_rel = (γ - 1) * m * (c**2 + (v**2)*(ω**2)*(r**2))
        ke_relativistic_results.append([ω, r, ke_rel])
        
        # Quantum particle‑in‑a‑box energy, modified by harmonic term
        e_q = (n**2 * np.pi**2 * ℏ**2) / (2 * m * L**2) * (1 + (ω**2)*(r**2)/(c**2))
        e_quantum_results.append([ω, r, e_q])
        
        # One‑step “chaos” logistic map with harmonic modulation
        x_new = r_chaos * x0 * (1 - x0) * (1 + (ω**2)*(r**2))
        chaos_results.append([ω, r, x_new])

# ─────────────── Build DataFrames ───────────────
df_ke_classical    = pd.DataFrame(ke_classical_results,    columns=["ω", "r", "KE_Classical"])
df_ke_relativistic = pd.DataFrame(ke_relativistic_results, columns=["ω", "r", "KE_Relativistic"])
df_e_quantum       = pd.DataFrame(e_quantum_results,       columns=["ω", "r", "E_Quantum"])
df_chaos           = pd.DataFrame(chaos_results,           columns=["ω", "r", "Chaos_Output"])

# ─────────────── Plotting ───────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Classical KE
sc0 = axes[0,0].scatter(df_ke_classical["ω"], df_ke_classical["r"],
                        c=df_ke_classical["KE_Classical"], cmap="inferno", edgecolor="k")
axes[0,0].set(title="Classical KE vs ω, r", xlabel="ω", ylabel="r")
fig.colorbar(sc0, ax=axes[0,0])

# Relativistic KE
sc1 = axes[0,1].scatter(df_ke_relativistic["ω"], df_ke_relativistic["r"],
                        c=df_ke_relativistic["KE_Relativistic"], cmap="plasma", edgecolor="k")
axes[0,1].set(title="Relativistic KE vs ω, r", xlabel="ω", ylabel="r")
fig.colorbar(sc1, ax=axes[0,1])

# Quantum Energy
sc2 = axes[1,0].scatter(df_e_quantum["ω"], df_e_quantum["r"],
                        c=df_e_quantum["E_Quantum"], cmap="viridis", edgecolor="k")
axes[1,0].set(title="Quantum Energy vs ω, r", xlabel="ω", ylabel="r")
fig.colorbar(sc2, ax=axes[1,0])

# Chaos Map Output
sc3 = axes[1,1].scatter(df_chaos["ω"], df_chaos["r"],
                        c=df_chaos["Chaos_Output"], cmap="magma", edgecolor="k")
axes[1,1].set(title="Chaos Map Output vs ω, r", xlabel="ω", ylabel="r")
fig.colorbar(sc3, ax=axes[1,1])

plt.tight_layout()
plt.show()

# ─────────────── Display DataFrames ───────────────
from IPython.display import display

print("First few rows of each sweep:")
display(df_ke_classical.head())
display(df_ke_relativistic.head())
display(df_e_quantum.head())
display(df_chaos.head())
