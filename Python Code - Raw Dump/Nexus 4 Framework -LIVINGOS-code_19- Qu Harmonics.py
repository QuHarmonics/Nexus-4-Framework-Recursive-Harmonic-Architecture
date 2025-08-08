import pandas as pd

# Load your packet hop log from Phase I
hop_log = pd.read_csv("nexus_packets_hop_log.csv")

# Build recursive identity stream from hop-level STI drift
def build_identity_stream(hop_log):
    streams = []
    grouped = hop_log.groupby("packet_id")

    for packet_id, group in grouped:
        group_sorted = group.sort_values("hop")

        echo_trace = "-".join(group_sorted["echo"].str[:3].tolist())
        sti_series = group_sorted["sti"].tolist()
        zphc_series = group_sorted["zphc"].tolist()

        collapse_points = [i for i, v in enumerate(zphc_series) if not v]
        recovered = any(zphc_series) and bool(collapse_points)

        stream = {
            "packet_id": packet_id,
            "stream_echo": echo_trace,
            "avg_sti": round(sum(sti_series) / len(sti_series), 4),
            "first_collapse": collapse_points[0] if collapse_points else None,
            "recovered": recovered,
            "final_zphc": zphc_series[-1],
            "drift_curve": sti_series
        }

        streams.append(stream)

    return pd.DataFrame(streams)

# Generate the identity stream
df_stream = build_identity_stream(hop_log)

# Save for future modeling or visualization
df_stream.to_csv("recursive_identity_stream.csv", index=False)
print("✅ Recursive identity stream saved as: recursive_identity_stream.csv")
