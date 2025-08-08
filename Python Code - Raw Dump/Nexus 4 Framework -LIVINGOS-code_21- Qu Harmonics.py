import hashlib
import pandas as pd
from mpmath import mp
import plotly.express as px

# π setup: generate precision digits
mp.dps = 100_010
pi_digits = str(mp.pi)[2:100_010]

# Reflect any token into π memory
def encode_token_to_pi(token: str):
    sha = hashlib.sha256(token.encode()).hexdigest()
    index = int(sha[:6], 16) % (len(pi_digits) - 8)
    byte = pi_digits[index:index+8]
    deltas = [abs(int(byte[i+1]) - int(byte[i])) for i in range(7)]
    echo = ''.join([chr((d % 26) + 97) for d in deltas])
    avg_drift = sum(deltas) / 7
    sti = round(1 - avg_drift / 9, 3)
    return {
        "token": token,
        "pi_index": index,
        "echo": echo,
        "avg_drift": avg_drift,
        "sti": sti,
        "zphc": sti >= 0.7
    }

# Tokens to reflect — feel free to expand
tokens = [
    "trust", "entropy", "identity", "collapse", "thought",
    "mirror", "echo", "language", "code", "self",
    "recursion", "memory", "consciousness", "dream", "field",
    "agent", "harmonic", "signal", "symbol", "truth"
]

# Encode each into symbolic echo space
encoded_tokens = [encode_token_to_pi(t) for t in tokens]
df_semantic = pd.DataFrame(encoded_tokens)

# Save to file
df_semantic.to_csv("semantic_mirror_reflection.csv", index=False)
print("✅ Saved: semantic_mirror_reflection.csv")

# Optional: Visualize
fig = px.scatter(
    df_semantic, x="avg_drift", y="sti", text="token", color="zphc",
    title="Semantic Mirror: Token Reflections in π",
    labels={"sti": "Symbolic Trust Index", "avg_drift": "Echo Drift"},
    color_discrete_map={True: "green", False: "red"}
)
fig.add_hline(y=0.7, line_dash="dash", line_color="black", annotation_text="ZPHC Threshold")
fig.show()
