# Initialize the character and base for "2+2=9"
input_text = "2+29=11"
base = 16

# Initialize lists for X, Y, Z values
x_values = []  # Gain (ratio of change for each bit length)
y_values = []  # Frequency (difference between initial and iterative bit lengths)
z_values = []  # Pythagorean result (a^2 + b^2 = c^2)

# Process each character in the input text
initial_bit_length = len(bin(ord(input_text[0]))[2:])  # Bit length of the first character
for i, char in enumerate(input_text):
    ascii_value = ord(char)  # ASCII value
    quotient = ascii_value // base  # Quotient in hex conversion
    remainder = ascii_value % base  # Remainder in hex conversion
    hex_value = f"{quotient:X}{remainder:X}"  # Hexadecimal representation

    # Gain: Ratio of change for bit lengths (current bit length vs previous bit length)
    current_bit_length = len(bin(ascii_value)[2:])  # Current bit length
    if i > 0:
        previous_bit_length = len(bin(ord(input_text[i - 1]))[2:])
        gain = current_bit_length / previous_bit_length if previous_bit_length != 0 else 0
    else:
        gain = 1  # No gain for the first character

    # Frequency: Difference between initial bit length and current bit length
    frequency = abs(current_bit_length - initial_bit_length)

    # Pythagorean result: sqrt(a^2 + b^2)
    a = quotient
    b = remainder
    c = math.sqrt(a**2 + b**2)

    # Append results
    x_values.append(gain)
    y_values.append(frequency)
    z_values.append(c)

# Create 3D plot
fig = go.Figure()

fig.add_trace(go.Scatter3d(
    x=x_values,
    y=y_values,
    z=z_values,
    mode='markers+lines',
    marker=dict(size=5, color=z_values, colorscale='Viridis', opacity=0.8),
    line=dict(color='blue', width=2),
    name="Hex Stack Gain-Frequency Analysis"
))

# Update layout for visualization
fig.update_layout(
    title=f"3D Analysis: Gain, Frequency, and Transition Mix for {input_text}",
    scene=dict(
        xaxis_title="Gain (Bit Length Ratio Change)",
        yaxis_title="Frequency (Difference in Bit Lengths)",
        zaxis_title="Transition Mix (a^2 + b^2 = c^2)"
    ),
    height=700,
    width=900
)

# Show the plot
fig.show()
