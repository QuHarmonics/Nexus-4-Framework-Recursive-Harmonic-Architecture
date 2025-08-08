import pandas as pd

# Load your saved logs (adjust paths if needed)
packets_log = pd.read_csv("nexus_packets_log.csv")
hop_log = pd.read_csv("nexus_packets_hop_log.csv")
consensus_votes = pd.read_csv("symbolic_consensus_votes.csv")
consensus_results = pd.read_csv("symbolic_consensus_results.csv")

# === Dataset 1 ===
# Echo Trajectory Sequences: STI drift over recursive hops
def build_trajectory_dataset(hop_log):
    sequences = []
    grouped = hop_log.groupby("packet_id")
    for pid, group in grouped:
        sequence = group.sort_values("hop")[["sti"]].values.flatten().tolist()
        zphc = group["zphc"].any()
        sequences.append({
            "packet_id": pid,
            "sti_sequence": sequence,
            "zphc_delivered": zphc,
            "length": len(sequence)
        })
    return pd.DataFrame(sequences)

# === Dataset 2 ===
# Consensus Judgment Votes: agent STI + agreement
def build_consensus_judgment_set(votes):
    return votes[[
        "agent_id", "seed", "echo_proposed", "echo_seen",
        "sti", "vote", "consensus"
    ]]

# === Dataset 3 ===
# Packet Summary Metadata: identity, collapse, trust
def build_packet_summary(packets_log):
    return packets_log[[
        "packet_id", "origin", "echo_chain", "avg_sti",
        "delivered", "hops"
    ]]

# Generate datasets
df_trajectory = build_trajectory_dataset(hop_log)
df_consensus_judgments = build_consensus_judgment_set(consensus_votes)
df_packet_summary = build_packet_summary(packets_log)

# Save output (adjust paths as needed)
df_trajectory.to_csv("echo_training_sti_trajectories.csv", index=False)
df_consensus_judgments.to_csv("echo_training_consensus_votes.csv", index=False)
df_packet_summary.to_csv("echo_training_packet_summary.csv", index=False)

print("✅ Training datasets saved:")
print("- echo_training_sti_trajectories.csv")
print("- echo_training_consensus_votes.csv")
print("- echo_training_packet_summary.csv")
