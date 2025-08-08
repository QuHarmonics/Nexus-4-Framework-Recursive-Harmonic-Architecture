import numpy as np
import matplotlib.pyplot as plt

# Extract data into a clean numerical array, dropping non-numeric entries
numeric_df = df.apply(pd.to_numeric, errors='coerce').dropna(axis=1, how='all')
collapse_array = numeric_df.values

# Determine when collapse stabilizes at [1] or [0] across steps
final_states = []
for row in collapse_array:
    row_clean = row[~np.isnan(row)]
    if len(row_clean) == 1:
        final_states.append(int(row_clean[0]))
    else:
        final_states.append(None)

# Count total steps until collapse completion
collapse_depth = len(collapse_array)

# Plot collapse length vs. step
fig, ax = plt.subplots(figsize=(10, 6))
step_lengths = [len(row[~np.isnan(row)]) for row in collapse_array]
ax.plot(step_lengths, marker='o', label='Collapse Width')
ax.set_title("Collapse Progression Over Steps")
ax.set_xlabel("Step")
ax.set_ylabel("Width of Row")
ax.grid(True)
ax.legend()

collapse_depth, final_states[-5:], fig
