# Input text: "Pla"
input_text = "Pla"

# Step 1: Convert characters to their ASCII decimal values
ascii_decimal_values = [ord(char) for char in input_text]

# Step 2: Calculate the ratio of change between consecutive ASCII values
ratios_of_change = [
    ascii_decimal_values[i + 1] / ascii_decimal_values[i] if ascii_decimal_values[i] != 0 else 0
    for i in range(len(ascii_decimal_values) - 1)
]

# Prepare x-axis for ratio of change (between each pair)
x_indices = list(range(len(ratios_of_change)))

# Step 3: Plot the ratio of change for each step
fig = go.Figure()

# Plot the ratio of change
fig.add_trace(go.Scatter(
    x=x_indices,
    y=ratios_of_change,
    mode='lines+markers',
    name="Ratios of Change",
    line=dict(color='orange', width=2)
))

# Update layout
fig.update_layout(
    title="Ratios of Change for ASCII Conversion (Text: 'Pla')",
    xaxis_title="Step Index (Pairwise Comparison)",
    yaxis_title="Ratio of Change",
    height=600,
    width=800
)

# Show the plot
fig.show()
