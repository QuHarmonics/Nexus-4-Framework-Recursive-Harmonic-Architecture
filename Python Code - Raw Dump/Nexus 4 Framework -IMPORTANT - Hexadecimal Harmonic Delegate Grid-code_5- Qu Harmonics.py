import pandas as pd

# Collapse logic — unchanged
def delta_collapse(data):
    history = [data[:]]
    current = data[:]
    while len(current) > 3:
        next_stage = [abs(current[i + 1] - current[i]) for i in range(len(current) - 1)]
        history.append(next_stage)
        current = next_stage
    return {
        "triplet": current,
        "depth": len(history) - 1,
        "x0": data[0],
        "delta0": abs(data[1] - data[0]) if len(data) > 1 else 0,
        "sum": sum(data)
    }

# Prompt for raw digit string and chunk size
raw_input = input("Paste your decimal digit string (e.g. 33353337...):\n").strip()
chunk_size = int(input("Enter chunk size (e.g. 4, 8, etc): "))

# Convert string to list of integers (digit-by-digit)
data = [int(d) for d in raw_input if d.isdigit()]

# Chunking logic
def chunk_sequence(data, chunk_size):
    padded_len = ((len(data) + chunk_size - 1) // chunk_size) * chunk_size
    padded = data + [0] * (padded_len - len(data))
    return [padded[i:i+chunk_size] for i in range(0, padded_len, chunk_size)]

# Chunk and compress
chunks = chunk_sequence(data, chunk_size)
compression_table = []
for chunk in chunks:
    result = delta_collapse(chunk)
    compression_table.append({
        "Original Chunk": chunk,
        "Triplet": result["triplet"],
        "Depth": result["depth"],
        "x0": result["x0"],
        "Δ0": result["delta0"],
        "Σx": result["sum"]
    })

# Display the result
df = pd.DataFrame(compression_table)
display(df)
