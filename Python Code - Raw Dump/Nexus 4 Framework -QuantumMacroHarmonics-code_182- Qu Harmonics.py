# Input sequence
decimal_sequence = [
    7, 90, 67, 255, 221, 153, 87, 114, 45, 30, 185, 36, 34, 142, 148, 81, 231, 237,
    187, 68, 55, 253, 210, 79, 188, 191, 101, 39, 126, 114, 189, 252
]

# Convert to hexadecimal
hex_sequence = [hex(num)[2:].zfill(2) for num in decimal_sequence]

# Convert to binary
binary_sequence = [bin(num)[2:].zfill(8) for num in decimal_sequence]

# Convert to ASCII (if valid printable ASCII range, else show as \xNN)
ascii_sequence = [
    chr(num) if 32 <= num <= 126 else f"\\x{num:02x}" for num in decimal_sequence
]

# Display results
print("Decimal Sequence:")
print(decimal_sequence)

print("\nHexadecimal Sequence:")
print(hex_sequence)

print("\nBinary Sequence:")
print(binary_sequence)

print("\nASCII Sequence:")
print(ascii_sequence)

# Write results to a file for better visualization
with open("sequence_representations.txt", "w") as f:
    f.write("Decimal Sequence:\n")
    f.write(" ".join(map(str, decimal_sequence)) + "\n\n")

    f.write("Hexadecimal Sequence:\n")
    f.write(" ".join(hex_sequence) + "\n\n")

    f.write("Binary Sequence:\n")
    f.write(" ".join(binary_sequence) + "\n\n")

    f.write("ASCII Sequence:\n")
    f.write("".join(ascii_sequence) + "\n")

print("\nResults saved to 'sequence_representations.txt'")
