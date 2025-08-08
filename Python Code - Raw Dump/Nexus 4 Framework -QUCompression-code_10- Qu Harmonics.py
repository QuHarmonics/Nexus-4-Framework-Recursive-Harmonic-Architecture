import numpy as np
import matplotlib.pyplot as plt

# ... (rest of the code remains the same)

# Reconstructed data
reconstructed_data = np.real(np.fft.ifft(reallocated_energy))

# Print reconstructed data
print("Reconstructed data:", reconstructed_data)

# Visualize reconstructed data
plt.plot(reconstructed_data)
plt.title("Reconstructed Data")
plt.xlabel("Index")
plt.ylabel("Value")
plt.show()

# Calculate statistical metrics
mean_value = np.mean(reconstructed_data)
std_dev = np.std(reconstructed_data)
variance = np.var(reconstructed_data)

print("Mean value:", mean_value)
print("Standard deviation:", std_dev)
print("Variance:", variance)