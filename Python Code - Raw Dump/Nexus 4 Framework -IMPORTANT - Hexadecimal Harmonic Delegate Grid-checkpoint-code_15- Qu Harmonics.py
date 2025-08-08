import pandas as pd

# Display settings for full visibility
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

# Parse into digits
data = [int(d) for d in raw_input if d.isdigit()]
chunks = chunk_sequence(data, chunk_size)

# Collapse each chunk
collapsed = [delta_collapse_to_pair(chunk) for chunk in chunks]

# Normalize max collapse depth
max_depth = max(len(c) for c in collapsed)

# Reformat as vertical-flow DataFrame
df_dict = {}

for i, col in enumerate(collapsed):
    label = f"Chunk {i} (x0={chunks[i][0]}, Σx={sum(chunks[i])})"
    col_extended = col + [[None]*len(col[-1])] * (max_depth - len(col))  # pad
    df_dict[label] = col_extended

# Final frame: each chunk is a column, collapse flows down
df = pd.DataFrame(df_dict)
display(df)
