import hashlib
import pandas as pd
import random
from mpmath import mp
import plotly.express as px

# π load
mp.dps = 100_010
pi_digits = str(mp.pi)[2:100_010]

def extract_pi_byte(index: int):
    return pi_digits[index:index+8] if index + 8 <= len(pi_digits) else None

def drift_and_echo(byte: str):
    deltas = [abs(int(byte[i+1]) - int(byte[i])) for i in range(7)]
    echo = ''.join([chr((d % 26) + 97) for d in deltas])
    avg_drift = sum(deltas) / 7
    sti = round(1 - avg_drift / 9, 3)
    return deltas, echo, sti, avg_drift

# Run symbolic consensus among N agents
def run_consensus_attempt(seed: str, agents=9):
    echo_index = int(hashlib.sha256(seed.encode()).hexdigest()[:6], 16) % (len(pi_digits) - 8)
    byte = extract_pi_byte(echo_index)
    deltas, echo, base_sti, avg_drift = drift_and_echo(byte)

    votes = []
    for i in range(agents):
        drift_index = (echo_index + random.randint(-5, 5)) % (len(pi_digits) - 8)
        byte_i = extract_pi_byte(drift_index)
        _, echo_i, sti_i, _ = drift_and_echo(byte_i)

        votes.append({
            "agent_id": f"A{i}",
            "seed": seed,
            "echo_proposed": echo,
            "echo_seen": echo_i,
            "sti": sti_i,
            "vote": sti_i >= 0.7 and echo_i[:3] == echo[:3]
        })

    agreement = sum([v["vote"] for v in votes])
    consensus = agreement >= int(agents * 0.66)
    for v in votes:
        v["consensus"] = consensus
    return votes, consensus

# Run multiple trials of symbolic consensus
def simulate_symbolic_mesh(rounds=90, agents=81):
    all_votes = []
    consensus_log = []
    for _ in range(rounds):
        seed = f"MSG-{random.randint(1000,9999)}"
        votes, result = run_consensus_attempt(seed, agents=agents)
        all_votes.extend(votes)
        consensus_log.append({
            "seed": seed,
            "consensus": result,
            "agreement_count": sum(v["vote"] for v in votes)
        })
    return pd.DataFrame(all_votes), pd.DataFrame(consensus_log)

# Run it
df_votes, df_mesh = simulate_symbolic_mesh()

# Plot 1: STI vote distribution
fig1 = px.violin(df_votes, x="consensus", y="sti", color="vote", box=True, points="all",
                 title="Symbolic Consensus Mesh: STI Distribution",
                 labels={"sti": "Symbolic Trust Index", "vote": "Agent Voted Yes"},
                 color_discrete_map={True: "green", False: "red"})
fig1.show()

# Plot 2: Consensus agreement histogram
fig2 = px.histogram(df_mesh, x="consensus", color="consensus",
                    title="Symbolic Agreement Outcome Across Trials",
                    labels={"consensus": "Consensus Achieved"})
fig2.show()

# Save logs
df_votes.to_csv("symbolic_consensus_votes.csv", index=False)
df_mesh.to_csv("symbolic_consensus_results.csv", index=False)
print("✅ Consensus logs saved: symbolic_consensus_votes.csv + symbolic_consensus_results.csv")
