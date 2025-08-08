# Parameter sweeps for further exploration

# Define parameter ranges
omega_values = np.linspace(0, 10, 20)  # Angular velocity sweep from 0 to 10
r_values = np.linspace(0.1, 5, 20)  # Radius sweep from 0.1 to 5

# Initialize result storage
ke_classical_results = []
ke_relativistic_results = []
e_quantum_results = []
chaos_results = []

# Sweep through values
for omega in omega_values:
    for r in r_values:
        # Classical KE
        ke_classical = (1/2) * m * (v**2 + omega**2 * r**2)
        ke_classical_results.append([omega, r, ke_classical])
        
        # Relativistic KE
        ke_relativistic = (γ - 1) * m * (c**2 + v**2 * omega**2 * r**2)
        ke_relativistic_results.append([omega, r, ke_relativistic])
        
        # Quantum Energy
        e_quantum = (n**2 * π**2 * ℏ**2) / (2 * m * L**2) * (1 + (omega**2 * r**2) / c**2)
        e_quantum_results.append([omega, r, e_quantum])
        
        # Chaos Logistic Map
        x_new = r_chaos * x * (1 - x) * (1 + omega**2 * r**2)
        chaos_results.append([omega, r, x_new])

# Convert to DataFrames
df_ke_classical = pd.DataFrame(ke_classical_results, columns=["ω", "r", "KE_Classical"])
df_ke_relativistic = pd.DataFrame(ke_relativistic_results, columns=["ω", "r", "KE_Relativistic"])
df_e_quantum = pd.DataFrame(e_quantum_results, columns=["ω", "r", "E_Quantum"])
df_chaos = pd.DataFrame(chaos_results, columns=["ω", "r", "Chaos_Map_Output"])

# Visualize results
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Classical KE plot
im1 = axes[0, 0].scatter(df_ke_classical["ω"], df_ke_classical["r"], c=df_ke_classical["KE_Classical"], cmap='inferno', edgecolor='k')
axes[0, 0].set_title("Classical KE vs. ω, r")
axes[0, 0].set_xlabel("ω (Angular Velocity)")
axes[0, 0].set_ylabel("r (Radius)")
fig.colorbar(im1, ax=axes[0, 0])

# Relativistic KE plot
im2 = axes[0, 1].scatter(df_ke_relativistic["ω"], df_ke_relativistic["r"], c=df_ke_relativistic["KE_Relativistic"], cmap='plasma', edgecolor='k')
axes[0, 1].set_title("Relativistic KE vs. ω, r")
axes[0, 1].set_xlabel("ω (Angular Velocity)")
axes[0, 1].set_ylabel("r (Radius)")
fig.colorbar(im2, ax=axes[0, 1])

# Quantum Energy plot
im3 = axes[1, 0].scatter(df_e_quantum["ω"], df_e_quantum["r"], c=df_e_quantum["E_Quantum"], cmap='viridis', edgecolor='k')
axes[1, 0].set_title("Quantum Energy vs. ω, r")
axes[1, 0].set_xlabel("ω (Angular Velocity)")
axes[1, 0].set_ylabel("r (Radius)")
fig.colorbar(im3, ax=axes[1, 0])

# Chaos Map Output plot
im4 = axes[1, 1].scatter(df_chaos["ω"], df_chaos["r"], c=df_chaos["Chaos_Map_Output"], cmap='magma', edgecolor='k')
axes[1, 1].set_title("Chaos Map Output vs. ω, r")
axes[1, 1].set_xlabel("ω (Angular Velocity)")
axes[1, 1].set_ylabel("r (Radius)")
fig.colorbar(im4, ax=axes[1, 1])

# Adjust layout and show
plt.tight_layout()
plt.show()

# Display results in tabular format
tools.display_dataframe_to_user(name="Classical KE Sweep", dataframe=df_ke_classical)
tools.display_dataframe_to_user(name="Relativistic KE Sweep", dataframe=df_ke_relativistic)
tools.display_dataframe_to_user(name="Quantum Energy Sweep", dataframe=df_e_quantum)
tools.display_dataframe_to_user(name="Chaos Map Sweep", dataframe=df_chaos)
