import pandas as pd
import random

# Load previously saved packet and memory datasets
allocator_df = pd.read_csv("nexus_pi_symbolic_allocator.csv")
packet_df = pd.read_csv("nexus_packets_log.csv")

# Extract ZPHC-approved memory blocks
trusted_blocks = allocator_df.query("zphc == True")["pi_index"].tolist()

# Simulate symbolic process execution
def simulate_nexus_os(agent_count=12, programs=36, trust_threshold=0.7):
    transactions = []

    for i in range(programs):
        agent = f"AGENT-{random.randint(1000, 9999)}"
        mem_zone = random.choice(trusted_blocks)
        packet = packet_df.sample(1).iloc[0]  # Random symbolic packet

        tx = {
            "program_id": f"PROC-{i:03d}",
            "agent_id": agent,
            "memory_index": mem_zone,
            "echo_chain": packet["echo_chain"],
            "avg_sti": packet["avg_sti"],
            "delivered": packet["delivered"],
            "exec_success": bool(packet["avg_sti"] >= trust_threshold and packet["delivered"])
        }
        transactions.append(tx)

    return pd.DataFrame(transactions)

# Run the simulation
df_os = simulate_nexus_os()

# Save result
df_os.to_csv("nexus_os_execution_log.csv", index=False)
print("✅ Saved: nexus_os_execution_log.csv")
