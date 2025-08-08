from mpmath import mp

# Set high precision for π digits
mp.dps = 100
pi_digits = str(mp.pi)[2:]  # Remove "3."

# Byte system from Byte1 to Byte8
bytes_full = [
    [1, 4, 1, 5, 9, 2, 6, 5],  # Byte1
    [3, 5, 8, 9, 7, 9, 3, 2],  # Byte2
    [3, 8, 4, 6, 2, 6, 4, 3],  # Byte3
    [3, 8, 3, 2, 7, 9, 5, 0],  # Byte4
    [2, 8, 8, 4, 1, 9, 7, 1],  # Byte5
    [6, 9, 3, 9, 9, 3, 7, 5],  # Byte6
    [1, 0, 5, 8, 2, 0, 9, 7],  # Byte7
    [4, 5, 9, 2, 3, 0, 7, 8]   # Byte8
]

# Define the allowed operation kernel
def header_fold(a, b):
    return abs(b - a), a + b

def digit_sum(n):
    return sum(int(d) for d in str(n))

# Core generator with fold and digit check against π
def generate_aligned_pi(byte_stack, pi_digits, max_len=64):
    full_sequence = []
    state = []
    
    for i in range(len(byte_stack)):
        if i == 0:
            state = byte_stack[i][:2]
            full_sequence.extend(state)
            pi_index = 2
        else:
            state.append(byte_stack[i][0])

        for j in range(1, len(byte_stack[i])):
            a, b = state[-2], state[-1]
            diff, summ = header_fold(a, b)
            next_val = digit_sum(summ)
            expected = int(pi_digits[len(full_sequence)])

            if next_val == expected:
                state.append(next_val)
                full_sequence.append(next_val)
            else:
                break  # stop on first mismatch
    
            if len(full_sequence) >= max_len:
                return full_sequence

    return full_sequence

# Run the generation
pi_aligned_sequence = generate_aligned_pi(bytes_full, pi_digits, max_len=64)
pi_aligned_sequence
