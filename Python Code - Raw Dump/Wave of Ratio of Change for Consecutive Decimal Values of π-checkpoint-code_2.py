import plotly.graph_objects as go
from mpmath import mp

# Step 1: Generate 1,000,000 Decimal Places of π
mp.dps = 1_000_000  # Set precision to 1,000,000
pi_str = str(mp.pi)[2:]  # Get π as a string and exclude "3."

# Convert to integers
def get_pi_decimals(precision):
    return [int(digit) for digit in pi_str[:precision]]

# Compute the ratio of change between consecutive digits
def compute_ratio_of_change(digits):
    changes = [abs(digits[i] - digits[i - 1]) for i in range(1, len(digits))]
    ratios = [changes[i] / changes[i - 1] if changes[i - 1] != 0 else 0 for i in range(1, len(changes))]
    return ratios

# Step 2: Apply Nexus 2 Framework Harmonization
def harmonize_data(data, harmonic_constant=0.35):
    return [(value - harmonic_constant) ** 2 for value in data]

# Generate data and compute ratios
decimals = get_pi_decimals(1_000_000)
ratios = compute_ratio_of_change(decimals)

# Harmonize the ratios
harmonized_ratios = harmonize_data(ratios)

# Step 3: Visualization
fig = go.Figure()

# Plot original ratios
fig.add_trace(go.Scatter(
    x=list(range(len(ratios[:10_000]))),  # Limit visualization to 10,000 points for clarity
    y=ratios[:10_000],
    mode='lines',
    name='Original Ratios of Change'
))

# Plot harmonized ratios
fig.add_trace(go.Scatter(
    x=list(range(len(harmonized_ratios[:10_000]))),  # Limit visualization to 10,000 points
    y=harmonized_ratios[:10_000],
    mode='lines',
    name='Harmonized Ratios of Change (H=0.35)'
))

# Update layout
fig.update_layout(
    title="Harmonized vs. Original Ratios of Change (First 10,000 Points)",
    xaxis_title="Index",
    yaxis_title="Ratio Value",
    template="plotly_dark",
    showlegend=True
)

# Show the plot
fig.show()
