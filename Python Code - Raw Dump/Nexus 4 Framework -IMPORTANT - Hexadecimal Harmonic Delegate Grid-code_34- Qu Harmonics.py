import matplotlib.pyplot as plt

# Combined sequence based on user description
sequence = [1, 3, 4, 3, 4, 1, 1, 1, 2, 2, 5, 5, 5, 9, 6, 6, 2, 4, 4, 6, 5, 5, 5, 4]

# Plot the sequence
plt.figure(figsize=(5, 5))
plt.plot(sequence, marker='o', linestyle='-', linewidth=2)
plt.title("Waveform of Combined Byte 1 Triadic Echo (User + Grok)")
plt.xlabel("Index")
plt.ylabel("Value")
plt.grid(True)
plt.xticks(range(len(sequence)))
plt.show()
