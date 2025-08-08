# Retry plotting with a simpler visualization
plt.figure(figsize=(12, 6))

# Generate a subset of the growth pattern for clarity
subset_size = min(10, len(growth_patterns))  # Limit the visualization to 10 hashes
for i, pattern in enumerate(growth_patterns[:subset_size]):
    plt.plot(pattern, label=f"Hash {i+1}")

plt.title("Hash Growth Patterns from Byte1 Seed (Subset)")
plt.xlabel("Index")
plt.ylabel("Byte Value")
plt.legend(loc="upper right", fontsize="small")
plt.grid(True)
plt.show()
