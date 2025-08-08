from typing import List
import pandas as pd

def delta_collapse(data: List[int]) -> List[List[int]]:
    """Apply recursive delta collapse until the list cannot shrink further."""
    history = [data[:]]  # Store all intermediate collapse stages
    current = data[:]
    
    while len(current) > 1:
        next_stage = [abs(current[i + 1] - current[i]) for i in range(len(current) - 1)]
        history.append(next_stage)
        if next_stage == current:
            break  # Reached steady state
        current = next_stage
    
    return history

# Example: replace with your own 8-digit sequence
input_digits = [ 1, 4, 1, 5, 9, 2, 6, 5]

# Run the collapse
collapsed = delta_collapse(input_digits)

# Pad results for display
max_len = max(len(row) for row in collapsed)
padded = [row + [None] * (max_len - len(row)) for row in collapsed]

# Create and display DataFrame
df = pd.DataFrame(padded)
df.columns = [f"Index {i}" for i in range(max_len)]
df.index.name = "Collapse Step"
df
