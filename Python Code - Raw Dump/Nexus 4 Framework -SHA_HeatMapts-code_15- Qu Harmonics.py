import matplotlib.pyplot as plt

# Suppose you have counts for digits '0'–'9':
digit_counts = {'0': 5, '1': 3, '2': 4, '3': 6, '4': 2, '5': 7, '6': 1, '7': 4, '8': 2, '9': 3}
widths = measure_digit_widths(FONT_PATH, FONT_SIZE)

# Prepare data
digits = list(digit_counts.keys())
counts = [digit_counts[d] for d in digits]
bar_widths = [widths[d] for d in digits]

# Normalize widths so the total width fits nicely
total = sum(bar_widths)
bar_widths_norm = [w / total * len(digits) for w in bar_widths]

# Plot
plt.figure(figsize=(8,4))
plt.bar(digits, counts, width=bar_widths_norm, edgecolor='black')
plt.xlabel("Digit")
plt.ylabel("Event Count")
plt.title("Digit Counts with Font‐Width Encoding")
plt.tight_layout()
plt.show()
