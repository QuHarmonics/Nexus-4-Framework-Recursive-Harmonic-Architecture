import matplotlib.pyplot as plt

# Define all rows
rows = [
    [1, 6, 7, 7, 7, 2, 1, 6],   # Row 1
    [1, 4, 5, 5, 6, 6, 6, 7],   # Row 2
    [1, 3, 3, 4, 4, 5, 5, 6],   # Row 3
    [1, 2, 3, 3, 3, 4, 4, 5],   # Row 4
    [1, 2, 2, 2, 3, 3, 3, 4],   # Row 5
    [1, 1, 2, 2, 2, 2, 3, 3],   # Row 6
    [1, 2, 4, 4, 5, 6, 6, 7],   # Derived 1
    [1, 3, 6, 7, 8, 9, 9, 10]   # Derived 2
]

# Process each row to compute recursive Δ-path
delta_rows = []
for row in rows:
    cumulative = [row[0]]
    for i in range(1, len(row)):
        delta = row[i] - row[i - 1]
        cumulative.append(cumulative[-1] + delta)
    delta_rows.append(cumulative)

# Plotting
plt.figure(figsize=(12, 8))
x_vals = [i + 1 for i in range(8)]

for idx, row in enumerate(delta_rows):
    plt.plot(x_vals, row, marker='o', label=f"Row {idx + 1}")

plt.title("Recursive Δ Plot for All Rows (Cumulative State from Local Differences)")
plt.xlabel("Time Step")
plt.ylabel("Cumulative Δ State")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Optional: print all cumulative delta rows
for idx, row in enumerate(delta_rows):
    print(f"Row {idx + 1} cumulative: {row}")
