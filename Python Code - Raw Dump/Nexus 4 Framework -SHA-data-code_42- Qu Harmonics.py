import seaborn as sns

# Mean bit activation across folded and non-folded states
bit_matrix = np.array(df['bit_pattern'].tolist())
folded_mask = df['folded'].values

mean_folded = bit_matrix[folded_mask].mean(axis=0)
mean_unfolded = bit_matrix[~folded_mask].mean(axis=0)
bit_delta = mean_folded - mean_unfolded

# Plot bit-level differences
plt.figure(figsize=(8, 5))
plt.bar(range(8), bit_delta, color='crimson', alpha=0.8)
plt.xlabel('Bit Position (7 = MSB)')
plt.ylabel('Δ Activation (Folded − Unfolded)')
plt.title('Bit Activation Difference: Folded vs Unfolded States')
plt.grid(True)
plt.tight_layout()
plt.show()
