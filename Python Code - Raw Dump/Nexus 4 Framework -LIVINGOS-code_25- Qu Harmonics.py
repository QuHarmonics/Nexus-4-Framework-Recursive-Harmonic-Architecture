import hashlib
import pandas as pd
import random
from mpmath import mp

# Setup π memory with precision
mp.dps = 100_010
pi_digits = str(mp.pi)[2:100_010]

# Helper: Extract 8-digit byte from π at a given index
def extract_pi_byte(index: int):
    return pi_digits[index:index+8] if index + 8 <= len(pi_digits) else None

# Helper: Compute echo and STI from byte
def drift_and_echo(byte: str):
    deltas = [abs(int(byte[i+1]) - int(byte[i])) for i in range(7)]
    echo = ''.join([chr((d % 26) + 97) for d in deltas])
    avg_drift = sum(deltas) / 7
    sti = round(1 - avg_drift / 9, 3)
    return deltas, echo, sti, avg_drift

# Single agent attempt to achieve ZPHC by seed mutation
def steering_attempt(seed: str, attempts=6):
    logs = []
    for attempt in range(attempts):
        mutated_seed = seed + f"-M{attempt}"
        sha = hashlib.sha256(mutated_seed.encode()).hexdigest()
        index = int(sha[:6], 16) % (len(pi_digits) - 8)
        byte = extract_pi_byte(index)
        deltas, echo, sti, drift = drift_and_echo(byte)

        logs.append({
            "agent_id": seed,
            "mutation": mutated_seed,
            "attempt": attempt,
            "pi_index": index,
            "echo": echo,
            "avg_drift": drift,
            "sti": sti,
            "zphc": sti >= 0.7
        })

        # Stop if ZPHC condition is met
        if sti >= 0.7:
            break
    return logs

# Run symbolic steering simulation across multiple agents
def simulate_echo_steering(agent_count=30):
    results = []
    seeds = [f"SEED-{random.randint(1000, 9999)}" for _ in range(agent_count)]
    for seed in seeds:
        results.extend(steering_attempt(seed))
    return pd.DataFrame(results)

# Execute and save the results
df_steering = simulate_echo_steering()
df_steering.to_csv("symbolic_steering_attempt_log.csv", index=False)
print("✅ Echo steering log saved to: symbolic_steering_attempt_log.csv")
