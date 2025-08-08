def compute_recursive_stacks(initial_header, cycles):
    """
    Recursive framework for computing sequence with multiple stacks.

    :param initial_header: Initial header bit pair (e.g., [1, 4]).
    :param cycles: Number of cycles to compute.
    :return: Hierarchical structure of past, present, and future stacks.
    """
    past = []  # Past stack
    present = initial_header.copy()  # Present stack
    future = []  # Future stack
    results = []  # To store all cycles

    for cycle in range(cycles):
        # Store the current cycle's state
        results.append({"past": past.copy(), "present": present.copy(), "future": future.copy()})

        # Grow the present stack
        stack_len = len(present)
        for _ in range(stack_len):
            present.append(stack_len)

        # Compute new values based on addition/subtraction rules
        if len(present) >= 2:
            present[-1] = present[-2] - present[-1]  # Subtract last two values
        if len(present) >= 4:
            pointer = present[-1]
            ref_index = -pointer - 1
            if abs(ref_index) <= len(present):
                new_value = present[ref_index] + pointer
                present.append(new_value)

        # Compute the new header bits via XOR
        if len(present) >= 2:
            new_header_bit1 = present[-1]
            new_header_bit2 = present[-2]
            new_header = [new_header_bit1 ^ new_header_bit2]

        # Update stacks
        past = present.copy()
        present = new_header.copy()
        future = []

    return results


# Example Usage:
initial_header = [1, 4]
num_cycles = 13
sequence_structure = compute_recursive_stacks(initial_header, num_cycles)

# Display each cycle's stack states
for cycle, stacks in enumerate(sequence_structure):
    print(f"Cycle {cycle + 1}:")
    print(f"  Past: {stacks['past']}")
    print(f"  Present: {stacks['present']}")
    print(f"  Future: {stacks['future']}")
