# Scatterplot: Bit flips vs ΔResidue magnitude
plt.figure(figsize=(8, 5))
plt.scatter(df_trans['bit_flip'], df_trans['delta_residue'], alpha=0.6, color='teal')
plt.xlabel('Bit Flip Count (x → x+1)')
plt.ylabel('ΔResidue (mod 64)')
plt.title('Residue Jump vs Bit Transitions')
plt.grid(True)
plt.tight_layout()
plt.show()
