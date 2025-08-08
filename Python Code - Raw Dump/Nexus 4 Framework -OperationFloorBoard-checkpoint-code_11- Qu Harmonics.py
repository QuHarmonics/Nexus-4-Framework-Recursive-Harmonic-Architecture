# Initialize the character and base
input_text = "2+2=4"
base = 16

# Initialize lists for X, Y, Z values
x_values = []  # ASCII to Hex ratios
y_values = []  # Bit Length Ratios
z_values = []  # Cumulative Hex Length over iterations

# Process each character in the input text
cumulative_hex_length = 0  # Track cumulative length of hex stack
for char in input_text:
    ascii_value = ord(char)  # ASCII value
    quotient = ascii_value // base  # Quotient in hex conversion
    remainder = ascii_value % base  # Remainder in hex conversion
    hex_value = f"{quotient:X}{remainder:X}"  # Hexadecimal representation
    
    # Calculate ASCII to Hex ratio
    ascii_to_hex_ratio = ascii_value / int(hex_value, 16) if int(hex_value, 16) != 0 else 0
    
    # Calculate Bit Length Ratios
    ascii_bit_length = len(bin(ascii_value)[2:])  # Bit length of ASCII value
    hex_bit_length = len(bin(quotient)[2:]) + len(bin(remainder)[2:])  # Sum of bit lengths in hex
    bit_length_ratio = ascii_bit_length / hex_bit_length if hex_bit_length != 0 else 0
    
    # Update cumulative hex length
    cumulative_hex_length += len(hex_value)
    
    # Append to lists
    x_values.append(ascii_to_hex_ratio)
    y_values.append(bit_length_ratio)
    z_values.append(cumulative_hex_length)

# Create 3D plot
fig = go.Figure()

fig.add_trace(go.Scatter3d(
    x=x_values,
    y=y_values,
    z=z_values,
    mode='markers+lines',
    marker=dict(size=5, color=z_values, colorscale='Viridis', opacity=0.8),
    line=dict(color='blue', width=2),
    name="Hex Stack Analysis"
))

# Update layout for visualization
fig.update_layout(
    title="3D Analysis of ASCII to Hex Stack Growth",
    scene=dict(
        xaxis_title="ASCII to Hex Ratio",
        yaxis_title="Bit Length Ratio",
        zaxis_title="Cumulative Hex Length"
    ),
    height=700,
    width=900
)

# Show the plot
fig.show()
