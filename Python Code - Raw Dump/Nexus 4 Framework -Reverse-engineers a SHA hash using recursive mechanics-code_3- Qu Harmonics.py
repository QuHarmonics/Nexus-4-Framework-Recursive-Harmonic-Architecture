from mpmath import mp

# Set precision to get sufficient digits of pi
mp.dps = 100002  # 10,000 digits of Pi plus the initial "3."

# Get the first 10,000 digits of Pi as a string
pi_digits = str(mp.pi)[2:100002]  # Slice to exclude the initial "3."

# Convert digits of pi into 8-bit hex
hex_output = [hex(int(pi_digits[i:i+2]))[2:].zfill(2) for i in range(0, len(pi_digits), 2)]

# Join hex values into byte format
hex_byte_string = " ".join(hex_output)

# Display a portion to verify correctness (first 500 characters)
print(hex_byte_string)
