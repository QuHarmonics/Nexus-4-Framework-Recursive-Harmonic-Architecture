; Current Stack State: [1, 4, 1, 5, 9, 2, 6]
; Pointer: At the last value (`6`)

; STEP 10: Add the first two bits to compute the next number
MOV R1, [Stack - 7] ; Load the first value (1)
MOV R2, [Stack - 6] ; Load the second value (4)
ADD R3, R1, R2      ; Compute R3 = First + Second (1 + 4 = 5)
PUSH R3             ; Push the result onto the stack

; def compute_sequence_recursive(stack, steps):
    """
    Recursive function to compute the sequence.

    :param stack: Current state of the stack.
    :param steps: Remaining number of steps to compute.
    :return: The updated stack after computation.
    """
    # Base case: No more steps to process
    if steps == 0:
        return stack

    # Step 1: Calculate LEN and add it LEN times
    stack_len = len(stack)
    for _ in range(stack_len):
        stack.append(stack_len)

    # Step 2: Update the pointer value (subtraction step)
    if len(stack) >= 2:
        pointer_value = stack[-1]
        prev_value = stack[-2]
        stack[-1] = pointer_value - prev_value

    # Step 3: Calculate the next value in the sequence
    if len(stack) >= 4:
        pointer = stack[-1]
        reference_index = -pointer - 1
        if abs(reference_index) <= len(stack):
            reference_value = stack[reference_index]
            next_value = reference_value + pointer
            stack.append(next_value)

    # Recursive call
    return compute_sequence_recursive(stack, steps - 1)


# Example Usage:
initial_seed = [1, 4]  # Initial stack
num_steps = 5  # Number of steps to compute
result = compute_sequence_recursive(initial_seed, num_steps)
print("Final sequence:", result)
