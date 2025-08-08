import pandas as pd
import matplotlib.pyplot as plt

# Step 1: Generate mirrored expressions
mirror_pairs = []
for a in range(20):
    for b in range(20):
        if a != b and (a + b) < 10 and (a + b) % 2 == 1:
            expr1 = f"{a}+{b}="
            expr2 = f"{b}+{a}="
            hex1 = expr1.encode().hex()
            hex2 = expr2.encode().hex()
            dec1 = int(hex1, 16)
            dec2 = int(hex2, 16)
            delta = abs(dec1 - dec2)
            mirror_pairs.append({
                "Expression 1": expr1,
                "Expression 2": expr2,
                "Hex 1": hex1.upper(),
                "Hex 2": hex2.upper(),
                "Decimal 1": dec1,
                "Decimal 2": dec2,
                "Delta": delta
            })

# Step 2: Organize into DataFrame
df_mirrors = pd.DataFrame(mirror_pairs).drop_duplicates(subset=["Delta"]).sort_values("Delta")

# Step 3: Prepare data for plotting
expr_labels = [f"{row['Expression 1']} / {row['Expression 2']}" for _, row in df_mirrors.iterrows()]
delta_values = df_mirrors["Delta"].values

# Step 4: Digital root calculation
def digital_root(n):
    while n >= 10:
        n = sum(int(d) for d in str(n))
    return n

delta_droots = [digital_root(d) for d in delta_values]

# Step 5: Plot the results
plt.figure(figsize=(14, 6))
plt.plot(delta_values, marker='o', label="Delta (Decimal Difference)")
plt.plot(delta_droots, marker='x', linestyle='--', label="Digital Root of Delta")
plt.xticks(ticks=range(len(expr_labels)), labels=expr_labels, rotation=45, ha='right')
plt.title("Mirrored Expression Drift and Digital Roots")
plt.xlabel("Expression Pairs")
plt.ylabel("Value")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
