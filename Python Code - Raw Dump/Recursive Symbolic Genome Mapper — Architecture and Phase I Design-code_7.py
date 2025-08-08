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

# Replace with your custom input
input_digits = [3,4,4,4,3,5,4,1,3,5,3,0,3,0,3,0,3,0,3,2,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,4,3,0,3,0,3,0,4,6,3,0,3,0,4,6,4,6,4,6,4,6,3,0,3,0,3,0,3,0,4,2,3,8,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,4,3,0,3,0,3,0,3,1,4,1,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,0,3,1,3,0,3,0,3,0,3,0
]

# Run the collapse
collapsed = delta_collapse(input_digits)

# Normalize row lengths for display
max_len = max(len(row) for row in collapsed)
padded = [row + [None] * (max_len - len(row)) for row in collapsed]

# Build DataFrame
df = pd.DataFrame(padded)
df.columns = [f"Index {i}" for i in range(max_len)]
df.index.name = "Collapse Step"

# Display in full width in Jupyter
pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)

# Show the result
df
