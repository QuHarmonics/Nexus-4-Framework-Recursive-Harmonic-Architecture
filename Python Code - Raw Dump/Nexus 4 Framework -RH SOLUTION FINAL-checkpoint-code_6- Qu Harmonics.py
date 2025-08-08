import plotly.graph_objects as go

# Generate primes up to 10000
def sieve_of_eratosthenes(limit):
    is_prime = [True] * (limit + 1)
    is_prime[0], is_prime[1] = False, False  # 0 and 1 are not primes
    primes = []
    for i in range(2, limit + 1):
        if is_prime[i]:
            primes.append(i)
            for multiple in range(i * i, limit + 1, i):
                is_prime[multiple] = False
    return primes

# Calculate gaps and ratios
primes = sieve_of_eratosthenes(10000)
prime_gaps = [primes[i] - primes[i - 1] for i in range(1, len(primes))]
prime_gap_ratios = [prime_gaps[i] / primes[i] for i in range(len(prime_gaps))]

# Create a scatter plot with gaps as nodes and primes as locations
fig = go.Figure()

# Add prime locations as nodes
fig.add_trace(go.Scatter(
    x=primes[:-1],  # Exclude the last prime since gaps are one less in length
    y=prime_gaps,
    mode='markers',
    name='Prime Gaps',
    marker=dict(size=8, color='blue'),
    text=[f"Prime: {p}, Gap: {g}" for p, g in zip(primes[:-1], prime_gaps)]
))

# Add lines connecting consecutive primes to show relationships
fig.add_trace(go.Scatter(
    x=primes[:-1],
    y=prime_gaps,
    mode='lines',
    name='Gap Connections',
    line=dict(color='blue', dash='dot')
))

# Add ratio nodes
fig.add_trace(go.Scatter(
    x=primes[:-1],
    y=prime_gap_ratios,
    mode='markers',
    name='Gap Ratios',
    marker=dict(size=8, color='orange'),
    text=[f"Prime: {p}, Ratio: {r:.4f}" for p, r in zip(primes[:-1], prime_gap_ratios)]
))

# Layout
fig.update_layout(
    title="Prime Gaps and Gap Ratios with Prime Locations",
    xaxis_title="Prime Numbers",
    yaxis_title="Gap Values and Ratios",
    legend=dict(font=dict(size=12)),
    template="plotly_dark",
    xaxis=dict(showgrid=True),
    yaxis=dict(showgrid=True)
)

fig.show()
