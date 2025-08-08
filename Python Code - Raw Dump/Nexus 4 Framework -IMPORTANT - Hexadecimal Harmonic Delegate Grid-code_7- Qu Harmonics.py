import pandas as pd

# Collapse logic: all transformation steps until triplet
def delta_collapse_all(data):
    history = [data[:]]
    current = data[:]
    while len(current) > 3:
        current = [abs(current[i+1] - current[i]) for i in range(len(current)-1)]
        history.append(current)
    return history

# Chunking logic
def chunk_sequence(data, chunk_size):
    padded_len = ((len(data) + chunk_size - 1) // chunk_size) * chunk_size
    padded = data + [0] * (padded_len - len(data))
    return [padded[i:i+chunk_size] for i in range(0, padded_len, chunk_size)]

# 🔥 REAL USER INPUT
raw_input = input("Paste your raw decimal digit string (e.g. 335733... no spaces):\n").strip()
chunk_size = int(input("Enter chunk size (e.g. 4, 8, etc): "))

# Parse into real digits (0–9 only, nothing else touched)
data = [int(d) for d in raw_input if d.isdigit()]

# Chunk and process
chunks = chunk_sequence(data, chunk_size)
collapse_rows = []

for chunk in chunks:
    history = delta_collapse_all(chunk)
    row = {
        "Original Chunk": chunk,
        "x0": chunk[0],
        "Δ0": abs(chunk[1] - chunk[0]) if len(chunk) > 1 else 0,
        "Σx": sum(chunk)
    }
    for i, step in enumerate(history[1:], start=1):
        row[f"Step {i}"] = step
    collapse_rows.append(row)

# Output table
df = pd.DataFrame(collapse_rows)
display(df)
