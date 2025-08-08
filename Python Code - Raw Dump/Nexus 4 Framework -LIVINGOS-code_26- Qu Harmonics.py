import hashlib
import pandas as pd
import random
from mpmath import mp

# π memory setup
mp.dps = 100_010
pi_digits = str(mp.pi)[2:100_010]

# Extract 8-digit π slice
def extract_pi_byte(index: int):
    return pi_digits[index:index+8] if index + 8 <= len(pi_digits) else None

# Compute echo vector and symbolic trust
def drift_and_echo(byte: str):
    deltas = [abs(int(byte[i+1]) - int(byte[i])) for i in range(7)]
    echo = ''.join([chr((d % 26) + 97) for d in deltas])
    avg_drift = sum(deltas) / 7
    sti = round(1 - avg_drift / 9, 3)
    return deltas, echo, sti, avg_drift

# Generate symbolic echo from seed
def generate_symbolic_echo(agent_seed):
    sha = hashlib.sha256(agent_seed.encode()).hexdigest()
    index = int(sha[:6], 16) % (len(pi_digits) - 8)
    byte = extract_pi_byte(index)
    deltas, echo, sti, drift = drift_and_echo(byte)
    return {
        "agent_id": agent_seed,
        "pi_index": index,
        "echo": echo,
        "sti": sti,
        "zphc": sti >= 0.7
    }

# Run a single symbolic negotiation round
def symbolic_negotiation(agent_count=9):
    agents = [f"AG-{random.randint(1000,9999)}" for _ in range(agent_count)]
    proposals = [generate_symbolic_echo(agent) for agent in agents]

    # First echo becomes the proposal
    leader_echo = proposals[0]["echo"][:3]
    consensus_votes = []

    for p in proposals:
        agree = p["echo"][:3] == leader_echo and p["sti"] >= 0.7
        consensus_votes.append({
            "agent_id": p["agent_id"],
            "echo": p["echo"],
            "sti": p["sti"],
            "zphc": p["zphc"],
            "agrees_with_proposal": agree
        })

    votes = sum(1 for v in consensus_votes if v["agrees_with_proposal"])
    consensus_achieved = votes >= int(agent_count * 0.66)

    return consensus_votes, leader_echo, consensus_achieved

# Simulate multiple rounds of symbolic debate
def simulate_symbolic_negotiation(rounds=30, agent_count=9):
    all_votes = []
    round_summaries = []

    for r in range(rounds):
        votes, echo, result = symbolic_negotiation(agent_count)
        for v in votes:
            v["round"] = r
            v["proposal"] = echo
            v["consensus"] = result
        all_votes.extend(votes)

        round_summaries.append({
            "round": r,
            "proposal": echo,
            "consensus_achieved": result,
            "agreement_count": sum(v["agrees_with_proposal"] for v in votes)
        })

    return pd.DataFrame(all_votes), pd.DataFrame(round_summaries)

# Run simulation
df_votes, df_summary = simulate_symbolic_negotiation()

# Save results
df_votes.to_csv("symbolic_negotiation_votes.csv", index=False)
df_summary.to_csv("symbolic_negotiation_summary.csv", index=False)
print("✅ Negotiation results saved to:")
print("- symbolic_negotiation_votes.csv")
print("- symbolic_negotiation_summary.csv")
