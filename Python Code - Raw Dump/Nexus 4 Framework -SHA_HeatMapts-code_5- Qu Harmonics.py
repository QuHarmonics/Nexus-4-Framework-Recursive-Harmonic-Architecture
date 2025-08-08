from scipy.stats import entropy
import numpy as np

# Assume bit_array is a 2D array [samples, bits] and folded_mask is a boolean array
folded_bits = bit_array[folded_mask]
unfolded_bits = bit_array[~folded_mask]

# Mean activation per bit, then entropy
folded_entropy = entropy(np.mean(folded_bits, axis=0), base=2)
unfolded_entropy = entropy(np.mean(unfolded_bits, axis=0), base=2)

# Delta entropy to highlight key bits
delta_entropy = folded_entropy - unfolded_entropy