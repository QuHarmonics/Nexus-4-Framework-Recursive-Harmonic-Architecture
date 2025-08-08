import pandas as pd

# Function: Perform only the FIRST collapse transformation (not recursive)
def first_delta_step(data):
    return [abs(data[i + 1] - data[i]) for i in range(len(data) - 1)]

# Prompt for user input (raw digit string) and chunk size
raw_input = input("Paste your decimal digit string (e.g. 33353337...):\n").strip()
chunk_size = int(input("Enter chunk size (e.g. 4, 8, etc): "))

# Parse raw input string into digits
data = [int(d) for d in raw_input if d.isdigit()]

# Chunking function
def chunk_sequence(data, chunk_size):
    padded_len = ((len(data) + chunk_size - 1) // chunk_size) * chunk_size
    padded = data + [0] * (padded_len - len(data))
    return [padded[i:i+chunk_size] for i in range(0, padded_len, chunk_size)]

# Build table: show only the first delta transformation per chunk
chunks = chunk_sequence(data, chunk_size)
transformation_table = []

for chunk in chunks:
    first_transform = first_delta_step(chunk)
    transformation_table.append({
        "Original Chunk": chunk,
        "1st Collapse Step": first_transform,
        "x0": chunk[0],
        "Δ0": abs(chunk[1] - chunk[0]) if len(chunk) > 1 else 0,
        "Σx": sum(chunk)
    })

df = pd.DataFrame(transformation_table)
display(df)
