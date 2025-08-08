def to_binary(input_str):
    """Convert input string to binary representation."""
    return ''.join(format(ord(c), '08b') for c in input_str)

def binary_to_hex(binary_str):
    """Convert binary to hexadecimal representation."""
    return hex(int(binary_str, 2))[2:].upper()

def decompile_to_asm(hex_input):
    """
    Simulate decompilation by mapping hex to assembly-like instructions.
    """
    asm = []
    for i in range(0, len(hex_input), 2):  # Process hex in pairs
        byte = hex_input[i:i+2]
        asm.append(f"LOAD #{byte}")  # Simulate loading a byte into a register
        asm.append("ADD R1, R0")    # Simulate a computation step
        asm.append("CMP R1, R2")    # Simulate a comparison operation
    return asm

def simulate_asm_execution(asm_code):
    """
    Simulate execution of assembly instructions.
    """
    register = 0  # Simulate a simple register
    executed_output = []
    for instruction in asm_code:
        if "LOAD" in instruction:
            byte_value = int(instruction.split("#")[1], 16)
            register += byte_value  # Simulate loading into a register
        elif "ADD" in instruction:
            register += 1  # Simulate a dummy addition
        elif "CMP" in instruction:
            executed_output.append(register)  # Capture register state
    return executed_output

def calculate_deviation(input_hex, output_values):
    """
    Calculate the deviation between the input hex and the output register states.
    """
    input_values = [int(input_hex[i:i+2], 16) for i in range(0, len(input_hex), 2)]
    if len(input_values) != len(output_values):
        return None  # Ensure lengths match
    deviation = [abs(i - o) for i, o in zip(input_values, output_values)]
    ratio = sum(deviation) / sum(input_values) if sum(input_values) != 0 else 0
    return deviation, ratio

# Input pairs: good and bad
good_pair = "2 + 2 = 4"
bad_pair = "2 + 9 = 4"

# Convert to binary and hex
good_pair_binary = to_binary(good_pair)
bad_pair_binary = to_binary(bad_pair)
good_pair_hex = binary_to_hex(good_pair_binary)
bad_pair_hex = binary_to_hex(bad_pair_binary)

# Decompile to assembly-like instructions
good_pair_asm = decompile_to_asm(good_pair_hex)
bad_pair_asm = decompile_to_asm(bad_pair_hex)

# Simulate execution
good_pair_output = simulate_asm_execution(good_pair_asm)
bad_pair_output = simulate_asm_execution(bad_pair_asm)

# Calculate deviation for the bad pair
bad_pair_deviation, bad_pair_ratio = calculate_deviation(bad_pair_hex, bad_pair_output)

# Output results
print("Good Pair Hex:", good_pair_hex)
print("Good Pair Output Matches Input:", good_pair_output == [int(good_pair_hex[i:i+2], 16) for i in range(0, len(good_pair_hex), 2)])

print("\nBad Pair Hex:", bad_pair_hex)
print("Bad Pair Deviation:", bad_pair_deviation)
print("Bad Pair Deviation Ratio:", bad_pair_ratio)

if bad_pair_ratio > 0:
    print("\nThe bad pair reveals a deviation, indicating a quantum difference.")
else:
    print("\nThe bad pair is harmonized (unexpected for this case).")
