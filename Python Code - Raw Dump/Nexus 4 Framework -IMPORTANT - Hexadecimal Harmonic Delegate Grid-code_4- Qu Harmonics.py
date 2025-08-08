import pandas as pd

# Define your exact symbolic 4-digit input chunks
selected_chunks = [
    [3, 3, 3, 3], [3, 3, 3, 2], [3, 3, 3, 2], [3, 3, 3, 2], [3, 3, 3, 2],
    [3, 3, 3, 4], [3, 3, 3, 4], [3, 3, 3, 4],
    [3, 3, 3, 1],
    [3, 3, 3, 5], [3, 3, 3, 5], [3, 3, 3, 5], [3, 3, 3, 5],
    [3, 3, 3, 0], [3, 3, 3, 0], [3, 3, 3, 0],
    [3, 3, 3, 6], [3, 3, 3, 6], [3, 3, 3, 6], [3, 3, 3, 6],
    [3, 3, 3, 6], [3, 3, 3, 6], [3, 3, 3, 6], [3, 3, 3, 6],
    [3, 3, 3, 7], [3, 3, 3, 7], [3, 3, 3, 7], [3, 3, 3, 7], [3, 3, 3, 7],
    [3, 3, 3, 9], [3, 3, 3, 9],
    [3, 4, 3, 3],
    [3, 4, 3, 5], [3, 4, 3, 5],
    [3, 4, 3, 6], [3, 4, 3, 6]
]

# Define delta collapse logic
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

# Build the compression table
compression_table = []
for chunk in selected_chunks:
    result = delta_collapse(chunk)
    compression_table.append({
        "Original Chunk": chunk,
        "Triplet": result["triplet"],
        "Depth": result["depth"],
        "x0": result["x0"],
        "Δ0": result["delta0"],
        "Σx": result["sum"]
    })

# Convert to DataFrame and display
df = pd.DataFrame(compression_table)
display(df)
