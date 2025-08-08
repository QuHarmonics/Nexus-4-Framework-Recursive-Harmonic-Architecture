import pandas as pd
import random

# Load pre-generated symbolic trajectory data (from Phase I)
trajectory_df = pd.read_csv("echo_training_sti_trajectories.csv")

# Function: simulate symbolic learning by agents
def simulate_echo_learning_agents(agent_count=6, samples_per_agent=12):
    agents = []

    for i in range(agent_count):
        agent_id = f"ECHO-AGENT-{1000 + i}"
        memory = []

        for _ in range(samples_per_agent):
            sample = trajectory_df.sample(1).iloc[0]

            # Parse symbolic trust curve
            trust_curve = sample["sti_sequence"].strip("[]").split(",")
            trust_curve = [float(s) for s in trust_curve]
            avg_trust = sum(trust_curve) / len(trust_curve)

            # First ZPHC occurrence
            first_zphc = next((idx for idx, val in enumerate(trust_curve) if val >= 0.7), None)

            memory.append({
                "agent_id": agent_id,
                "packet_id": sample["packet_id"],
                "avg_trust": round(avg_trust, 4),
                "final_zphc": sample["zphc_delivered"],
                "first_lock": first_zphc if first_zphc is not None else -1,
                "trust_signature": trust_curve
            })

        agents.extend(memory)

    return pd.DataFrame(agents)

# Run it
df_learners = simulate_echo_learning_agents()

# Save result
df_learners.to_csv("echo_agents_symbolic_learning_log.csv", index=False)
print("✅ Echo agents learning saved to: echo_agents_symbolic_learning_log.csv")
