import matplotlib.pyplot as plt

# Time points and corresponding R values
phi = 0.61803398875
H = 0.35
F = 1.0

time_values = [round(phi * i, 6) for i in range(200)]
R_values = [np.exp(H * F * t) for t in time_values]

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(time_values, R_values, marker='o', linestyle='-', linewidth=2)
plt.title("Harmonic Recursive Tone Growth", fontsize=14)
plt.xlabel("t (Golden-Stepped Time)", fontsize=12)
plt.ylabel("R (Recursive Reflection Magnitude)", fontsize=12)
plt.grid(True)
plt.xticks(time_values, rotation=45)
plt.tight_layout()
plt.show()
