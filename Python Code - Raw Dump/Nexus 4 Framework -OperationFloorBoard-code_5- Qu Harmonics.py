import math
import plotly.graph_objects as go

def text_to_waveform_analysis(text):
    # Convert text to hexadecimal (ASCII values)
    hex_values = [ord(c) for c in text]  # ASCII to decimal equivalent
    
    # Step 1: Compute deltas (changes between values)
    deltas = [hex_values[i + 1] - hex_values[i] for i in range(len(hex_values) - 1)]

    # Step 2: Compute Pythagorean relationships between consecutive deltas
    pythagorean_results = []
    for i in range(len(deltas) - 1):
        a = deltas[i]
        b = deltas[i + 1]
        c = math.sqrt(a**2 + b**2)
        pythagorean_results.append(c)

    # Step 3: Nesting levels (simulate by grouping values, e.g., pairs)
    nesting_level_results = []
    for i in range(0, len(hex_values), 2):  # Process in pairs
        group = hex_values[i:i + 2]
        if len(group) == 2:
            a, b = group
            nesting_level_results.append(math.sqrt(a**2 + b**2))

    return {
        "hex_values": hex_values,
        "deltas": deltas,
        "pythagorean_results": pythagorean_results,
        "nesting_level_results": nesting_level_results,
    }


def plot_waveform_analysis(result, text):
    # Extract results for plotting
    hex_values = result["hex_values"]
    deltas = result["deltas"]
    pythagorean_results = result["pythagorean_results"]
    nesting_level_results = result["nesting_level_results"]

    # Create subplots for analysis
    fig = go.Figure()

    # Plot Hexadecimal Values
    fig.add_trace(go.Scatter(
        y=hex_values,
        mode='lines+markers',
        name='Hex Values',
        line=dict(color='blue')
    ))

    # Plot Deltas
    fig.add_trace(go.Scatter(
        y=deltas,
        mode='lines+markers',
        name='Deltas',
        line=dict(color='red')
    ))

    # Plot Pythagorean Results
    fig.add_trace(go.Scatter(
        y=pythagorean_results,
        mode='lines+markers',
        name='Pythagorean Results',
        line=dict(color='green')
    ))

    # Plot Nesting Levels
    fig.add_trace(go.Scatter(
        y=nesting_level_results,
        mode='lines+markers',
        name='Nesting Level Results',
        line=dict(color='orange')
    ))

    # Update layout
    fig.update_layout(
        title=f"Waveform Analysis of Text: '{text}'",
        xaxis_title="Index",
        yaxis_title="Value",
        legend=dict(orientation="h", y=-0.2),
        height=600,
        width=800
    )

    return fig


# Input text for analysis
text = "hello"

# Perform waveform analysis
result = text_to_waveform_analysis(text)

# Plot results
fig = plot_waveform_analysis(result, text)
fig.show()
