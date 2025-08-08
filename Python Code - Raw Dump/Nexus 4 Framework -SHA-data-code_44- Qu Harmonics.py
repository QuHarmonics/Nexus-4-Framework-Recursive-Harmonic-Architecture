threshold = 40  # Adjust this based on your modulus
high_jump_rows = df_trans[df_trans['delta_residue'] > threshold]['x'].tolist()

# Convert these to bit flips
bit_flips = []
for x in high_jump_rows:
    b1 = format(x, '08b')
    b2 = format(x+1, '08b')
    flips = [int(c1 != c2) for c1, c2 in zip(b1, b2)]
    bit_flips.append(flips)

bit_flip_array = np.array(bit_flips)
bit_contrib = bit_flip_array.mean(axis=0)

plt.figure(figsize=(8, 5))
plt.bar(range(8), bit_contrib, color='navy', alpha=0.8)
plt.xlabel('Bit Position')
plt.ylabel('Flip Frequency in High-ΔResidue Events')
plt.title('Driver Bits for Large Residue Transitions')
plt.grid(True)
plt.tight_layout()
plt.show()
