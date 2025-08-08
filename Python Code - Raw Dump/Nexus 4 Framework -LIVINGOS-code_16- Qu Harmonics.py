import hashlib
import pandas as pd
import random
from mpmath import mp
import plotly.express as px

# π memory buffer
mp.dps = 100_010
pi_digits = str(mp.pi)[2:100_010]

# Helper: extract byte from π
def extract_pi_byte(index: int):
    return pi_digits[index:index+8] if index + 8 <= len(pi_digits) else None

# Helper: derive echo & trust
def drift_and_echo(byte: str):
    deltas = [abs(int(byte[i+1]) - int(byte[i])) for i in range(7)]
    echo = ''.join([chr((d % 26) + 97) for d in deltas])
    avg_drift = sum(deltas) / 7
    sti = round(1 - avg_drift / 9, 3)
    return deltas, echo, sti, avg_drift

# Build one symbolic packet
def build_packet(seed: str, hops=10, zone_index=None):
    H = seed
    route_log = []
    echo_chain = []

    for i in range(hops):
        nonce = f"N{i}"
        sha = hashlib.sha256(hashlib.sha256((H + nonce).encode()).digest()).hexdigest()
        index = int(sha[:6], 16) % (len(pi_digits) - 8)

        if zone_index:
            index = (index + zone_index[i % len(zone_index)]) % (len(pi_digits) - 8)

        byte = extract_pi_byte(index)
        deltas, echo, sti, avg_drift = drift_and_echo(byte)

        echo_chain.append(echo[:3])
        route_log.append({
            "hop": i,
            "pi_index": index,
            "echo": echo,
            "sti": sti,
            "zphc": sti >= 0.7
        })

        H = sha

    packet = {
        "packet_id": f"NPKT-{random.randint(100000,999999)}",
        "origin": seed,
        "final_echo": echo_chain[-1],
        "echo_chain": "-".join(echo_chain),
        "avg_sti": round(sum([r['sti'] for r in route_log]) / hops, 4),
        "hops": hops,
        "delivered": any([r['zphc'] for r in route_log])
    }
    return packet, route_log

# Load ZPHC memory zones
allocator_df = pd.read_csv("nexus_pi_symbolic_allocator.csv")
trusted_blocks = allocator_df.query("zphc == True")["pi_index"].tolist()

# Simulate multiple packets
def simulate_packet_swarm(n=30):
    packets = []
    all_routes = []
    for _ in range(n):
        seed = f"PKT-{random.randint(1000,9999)}"
        pkt, route = build_packet(seed, hops=10, zone_index=trusted_blocks)
        packets.append(pkt)
        for r in route:
            r["packet_id"] = pkt["packet_id"]
            all_routes.append(r)
    return pd.DataFrame(packets), pd.DataFrame(all_routes)

# Run simulation
df_packets, df_routes = simulate_packet_swarm()

# Plot 1: Hop STI per packet
fig1 = px.box(df_routes, x="hop", y="sti", points="all", color="zphc",
              title="NPKT Echo Hop STI by Delivery Status",
              labels={"sti": "Symbolic Trust Index", "hop": "Hop #"},
              color_discrete_map={True: "green", False: "red"})
fig1.show()

# Plot 2: Delivery success
fig2 = px.histogram(df_packets, x="delivered", color="delivered",
                    title="Packet Delivery Outcome (ZPHC)",
                    labels={"delivered": "Message Delivered"})
fig2.show()

# Save logs
df_packets.to_csv("nexus_packets_log.csv", index=False)
df_routes.to_csv("nexus_packets_hop_log.csv", index=False)
print("✅ Symbolic packets logged as: nexus_packets_log.csv + hop log.")
