import pandas as pd

# Make sure all rows and columns are visible
pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.width", None)
pd.set_option("display.max_colwidth", None)

# Collapse until only 2 values remain
def delta_collapse_to_pair(data):
    history = [data[:]]
    current = data[:]
    while len(current) > 2:
        current = [abs(current[i + 1] - current[i]) for i in range(len(current) - 1)]
        history.append(current)
    return history

# Chunk the sequence
def chunk_sequence(data, chunk_size):
    padded_len = ((len(data) + chunk_size - 1) // chunk_size) * chunk_size
    padded = data + [0] * (padded_len - len(data))
    return [padded[i:i+chunk_size] for i in range(0, padded_len, chunk_size)]

# 🔥 User-controlled input (raw digits only)
raw_input = input("Paste your raw decimal digit string (e.g. 335733... no spaces):\n").strip()
chunk_size = int(input("Enter chunk size (e.g. 4, 8, etc): "))

# Parse into real digits
data = [int(d) for d in raw_input if d.isdigit()]

# Process chunks
chunks = chunk_sequence(data, chunk_size)
collapse_rows = []

for chunk in chunks:
    history = delta_collapse_to_pair(chunk)
    row = {
        "Original Chunk": chunk,
        "x0": chunk[0],
        "Δ0": abs(chunk[1] - chunk[0]) if len(chunk) > 1 else 0,
        "Σx": sum(chunk)
    }
    for i, step in enumerate(history[1:], start=1):
        row[f"Step {i}"] = step
    collapse_rows.append(row)

# Display result
df = pd.DataFrame(collapse_rows)
display(df)
