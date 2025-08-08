import hashlib, random, string, struct, math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ---------- helper ----------
def sha256_hex(msg: bytes) -> str:
    return hashlib.sha256(msg).hexdigest()

def digest_words(digest_hex: str):
    """Return list of 8 unsigned 32‑bit ints and their drifts to 2^32."""
    words = [int(digest_hex[i:i+8], 16) for i in range(0, 64, 8)]
    drifts = [(1<<32) - w for w in words]          # unsigned drift
    signed = [w - (1<<32) for w in words]          # negative drift
    return words, drifts, signed

# ---------- 1. Large random sample ----------
SAMPLE_N = 2000
rand_strings = [''.join(random.choices(string.ascii_letters+string.digits, k=random.randint(5,20)))
                for _ in range(SAMPLE_N)]

all_drifts=[]
for s in rand_strings:
    _, drifts, _ = digest_words(sha256_hex(s.encode()))
    all_drifts.extend(drifts)

# histogram
plt.figure(figsize=(6,4))
plt.hist(all_drifts, bins=60)
plt.title("Δ2^32 distribution – 2000 random strings")
plt.xlabel("drift (target − value)")
plt.ylabel("frequency")
plt.tight_layout()

# ---------- 2. Mirror‑pair case‑flip test ----------
PAIR_N = 500
pairs=[]
correlations=[]
for _ in range(PAIR_N):
    word = ''.join(random.choices(string.ascii_lowercase, k=random.randint(5,10)))
    mirror = word.capitalize()  # first char upper
    for label, txt in [(0,word),(1,mirror)]:
        digest=sha256_hex(txt.encode())
        _, drifts, signed = digest_words(digest)
        pairs.append({"pair_id":_, "case":label, **{f"d{i}":d for i,d in enumerate(signed)}})
    # correlation per pair
    dA = np.array(digest_words(sha256_hex(word.encode()))[2])
    dB = np.array(digest_words(sha256_hex(mirror.encode()))[2])
    correlations.append(abs(np.corrcoef(dA, dB)[0,1]))

pair_df = pd.DataFrame(pairs)
tools.display_dataframe_to_user("First 10 case‑flip pairs (signed drifts)", pair_df.head(20))

plt.figure(figsize=(5,5))
plt.scatter([abs(s) for s in pair_df[pair_df['case']==0]['d0']][:100],
            [abs(s) for s in pair_df[pair_df['case']==1]['d0']][:100])
plt.title("Signed drift word 0: lower vs capitalised (first 100 pairs)")
plt.xlabel("|signed drift| lower")
plt.ylabel("|signed drift| capitalised")
plt.tight_layout()

# summary correlation
mean_corr = sum(correlations)/len(correlations)

# ---------- 3. Recursive feedback ----------
INIT_MSG = b"1"
ITER = 20
rms_list=[]
msg = INIT_MSG
for i in range(ITER):
    digest = sha256_hex(msg)
    _, drifts, signed = digest_words(digest)
    rms = math.sqrt(sum(d*d for d in signed)/len(signed))
    rms_list.append(rms)
    # feed signed drifts back as bytes (big‑endian int32)
    msg = b"".join(struct.pack(">i", s) for s in signed)

plt.figure(figsize=(6,4))
plt.plot(range(1,ITER+1), rms_list, marker='o')
plt.title("RMS(|signed drift|) over recursive SHA iteration")
plt.xlabel("iteration")
plt.ylabel("RMS drift")
plt.tight_layout()

# show numeric table for recursion
rec_df = pd.DataFrame({"iteration": range(1,ITER+1), "RMS_drift": rms_list})
tools.display_dataframe_to_user("Recursive SHA RMS drift", rec_df)

print(f"Average correlation (signed drift vector) between case‑flip pairs: {mean_corr:.3f}")
