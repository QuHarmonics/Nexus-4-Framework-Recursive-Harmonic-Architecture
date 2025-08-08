# Example round constants (first few K-values from SHA-256)
round_constants = [0x428a2f98, 0x71374491, 0xb5c0fbcf, ...]  # Full list omitted for brevity

# Check residue values from data (e.g., a DataFrame 'df')
matched_K = [k for k in round_constants if k % 64 in df['residue'].values]
print("Aligned K-values:", matched_K)