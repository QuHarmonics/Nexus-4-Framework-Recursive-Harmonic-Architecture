binary_string = "01100001 01100010 01100011 1"  # Original binary string
zeros_to_add = 423  # Number of zeros to add

# Remove spaces for continuous binary string
binary_string = binary_string.replace(" ", "")

# Add the zeros
padded_binary = binary_string + ("0" * zeros_to_add)

# Format the output for readability (optional)
formatted_output = " ".join(padded_binary[i:i+8] for i in range(0, len(padded_binary), 8))

# Print the final output
print("Final Padded Binary String:")
print(formatted_output)

# Check the length
print(f"Length of final binary string: {len(padded_binary)} bits")