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

# Input: Your full recursive byte-structured sequence (can be replaced with any digits)
input_digits = [
    3,3,3,1,6,6,3,8
]

# Run collapse
collapsed = delta_collapse(input_digits)

# Print each step
for step_num, stage in enumerate(collapsed):
    print(f"Step {step_num}: {stage}")
