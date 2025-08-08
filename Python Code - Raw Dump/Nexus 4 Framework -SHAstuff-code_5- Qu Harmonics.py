import hashlib, random, string, struct, math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

### Helper ###
def sha256_hex(msg: bytes) -> str:
    return hashlib.sha256(msg).hexdigest()

def digest_words(hex_digest: str):
    return [hex_digest[i:i+8] for i in range(0, 64, 8)]

def drifts_from_digest(hex_digest: str):
    words = digest_words(hex_digest)
    drifts = []
    for w in words:
        val = int(w, 16)
        drift = (1 << 32) - val  # Δ2^32
        drifts.append(drift)
    return drifts

### Experiment 1: large random sample ###
N = 5000
drift_values = []
for _ in range(N):
    length = random.randint(5, 20)
    s = ''.join(random.choice(string.ascii_letters + string.digits) for _ in range(length))
    digest = sha256_hex(s.encode())
    drift_values.extend(drifts_from_digest(digest))

# Summary stats
exp1_stats = {
    "samples": N,
    "words": len(drift_values),
    "mean_drift": np.mean(drift_values),
    "std_drift": np.std(drift_values),
    "min_drift": np.min(drift_values),
    "max_drift": np.max(drift_values)
}

plt.figure(figsize=(6,4))
plt.hist(drift_values, bins=60)
plt.title("Experiment 1: Histogram of Δ2^32 (random strings)")
plt.xlabel("Drift")
plt.ylabel("Frequency")
plt.tight_layout()

### Experiment 2: mirror case pairs ###
M = 1000
pair_records = []
for _ in range(M):
    length = random.randint(5, 15)
    base = ''.join(random.choice(string.ascii_lowercase) for _ in range(length))
    mirror = base.upper()
    d1 = drifts_from_digest(sha256_hex(base.encode()))
    d2 = drifts_from_digest(sha256_hex(mirror.encode()))
    # record sum of corresponding drifts
    sums = [a + b for a, b in zip(d1, d2)]
    pair_records.append({
        "base": base,
        "mirror": mirror,
        "mean_sum": np.mean(sums),
        "max_abs_sum": np.max(np.abs(sums))
    })
pair_df = pd.DataFrame(pair_records)

plt.figure(figsize=(6,4))
plt.hist(pair_df["mean_sum"], bins=40)
plt.title("Experiment 2: Mean (Δ_base + Δ_upper) distribution")
plt.xlabel("Mean sum per digest")
plt.ylabel("Count")
plt.tight_layout()

### Experiment 3: recursive feedback ###
def drifts_to_bytes(drifts):
    # pack each signed drift into 4‑byte big‑endian signed int
    out = b''.join(struct.pack(">i", (d - (1<<32)) if d > (1<<31) else -d) for d in drifts)
    return out

iterations = 6
msg = b"1"
rec_stats = []
for i in range(iterations):
    digest_hex = sha256_hex(msg)
    drifts = drifts_from_digest(digest_hex)
    mean_abs = np.mean(np.abs(drifts))
    rec_stats.append({"iter": i, "mean_abs_drift": mean_abs})
    msg = drifts_to_bytes(drifts)  # feedback

rec_df = pd.DataFrame(rec_stats)

plt.figure(figsize=(6,4))
plt.plot(rec_df["iter"], rec_df["mean_abs_drift"], marker='o')
plt.title("Experiment 3: Mean |Δ| over recursion")
plt.xlabel("Iteration")
plt.ylabel("Mean |drift|")
plt.tight_layout()
