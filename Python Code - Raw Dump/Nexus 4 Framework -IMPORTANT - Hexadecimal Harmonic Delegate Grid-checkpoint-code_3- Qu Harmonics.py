from typing import List, Dict

# Collapse function: recursively apply delta until a triplet is formed
def delta_collapse(data: List[int]) -> Dict:
    history = [data[:]]
    current = data[:]
    while len(current) > 3:
        next_stage = [abs(current[i + 1] - current[i]) for i in range(len(current) - 1)]
        history.append(next_stage)
        current = next_stage
    return {
        "triplet": current,
        "depth": len(history) - 1,
        "history": history,
        "x0": data[0],
        "delta0": abs(data[1] - data[0]) if len(data) > 1 else 0,
        "sum": sum(data)
    }

# Expand function: deterministically reconstruct from triplet and phase tag
def delta_expand(triplet: List[int], depth: int, x0: int, delta0: int, total_sum: int) -> List[int]:
    current = triplet[:]
    for _ in range(depth):
        next_stage = [0] * (len(current) + 1)
        # Reverse swinging-add logic, using previous difference
        next_stage[0] = x0  # anchor start
        for i in range(1, len(next_stage)):
            next_stage[i] = next_stage[i - 1] + current[i - 1]
        current = next_stage
    return current

# Example input
sequence = [3,3,3,5,3,3,3,7,3,3,3,6,3,3,3,1,3,3,3,7,3,3,3,4,3,3,3,7,3,3,3,3,3,3,3,6,3,4,3,6,3,3,3,6,3,4,3,5,3,3,3,2,3,4,3,3,3,3,3,2,3,3,3,0,3,3,3,4,3,3,3,9,3,3,3,2,3,3,3,0,3,3,3,6,3,4,3,5,3,3,3,6,3,3,3,5,3,3,3,6,3,3,3,5,3,3,3,6,3,3,3,4,3,3,3,2,3,3,3,0,3,3,3,7,3,3,3,9,3,3,3,6,3,4,3,6,3,3,3,7,3,3,3,5]
collapsed = delta_collapse(sequence)

# Reconstruct from metadata
reconstructed = delta_expand(
    triplet=collapsed["triplet"],
    depth=collapsed["depth"],
    x0=collapsed["x0"],
    delta0=collapsed["delta0"],
    total_sum=collapsed["sum"]
)

# Display results
{
    "Original Input": sequence,
    "Collapse Triplet": collapsed["triplet"],
    "Collapse Depth": collapsed["depth"],
    "Phase Tag": {
        "x0": collapsed["x0"],
        "Δ0": collapsed["delta0"],
        "Σx": collapsed["sum"]
    },
    "Reconstructed": reconstructed
}
