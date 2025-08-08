# Full working example of the bidirectional collapse/expand codec system

from typing import List, Tuple

def delta_collapse(data: List[int]) -> Tuple[List[int], List[int]]:
    """
    Collapse the input data using left-to-right subtraction.
    Returns the delta array and the final rightmost value.
    """
    if len(data) < 2:
        return [], data[0] if data else 0

    delta = [data[i+1] - data[i] for i in range(len(data)-1)]
    final = data[-1]
    return delta, final

def delta_expand(delta: List[int], final_value: int) -> List[int]:
    """
    Expand the data using right-to-left subtraction.
    Reconstructs the original list from the delta and the final value.
    """
    result = [final_value]
    for d in reversed(delta):
        result.insert(0, result[0] - d)
    return result

# Example 1: Basic test
original = [3, 4, 4, 4]

# Collapse the sequence
delta, last = delta_collapse(original)
print("Collapsed Delta:", delta)
print("Stored Final Digit:", last)

# Expand it back
reconstructed = delta_expand(delta, last)
print("Reconstructed Sequence:", reconstructed)

# Example 2: Symbolic insight
# Simulate a rolling math or symbolic encoding with phase reference
def symbolic_expand(delta: List[int], final_value: int) -> List[int]:
    """
    A variant that treats the expansion as a stitched oscillating path
    where direction and subtraction simulate a phase-based reconstruction.
    """
    result = [final_value]
    for i in range(len(delta) - 1, -1, -1):
        phase = (-1)**i  # + or - depending on the index parity
        result.insert(0, result[0] - delta[i] * phase)
    return result

# Run the symbolic expansion for creative analysis
symbolic_result = symbolic_expand(delta, last)
print("Symbolic Oscillating Expansion:", symbolic_result)
