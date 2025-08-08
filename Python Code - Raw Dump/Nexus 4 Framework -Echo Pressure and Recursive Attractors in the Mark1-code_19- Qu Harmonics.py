
import argparse, hashlib, os, random, sys, time
from pathlib import Path
from typing   import Iterable, List, Tuple

# ─────────────────────────────────────────────────────────────────────
class EchoScanner:
    """Detect SHA-256 digests whose head mirrors the tail."""

    FRONT_LEN = 16           # compare first 16 bytes ↔ reversed last 16 bytes

    @staticmethod
    def sha256(data: bytes) -> bytes:
        return hashlib.sha256(data).digest()

    @staticmethod
    def echo_score(digest: bytes, n: int = FRONT_LEN) -> float:
        front, back = digest[:n], digest[-n:]
        reflect = sum(1 for a, b in zip(front, reversed(back)) if a == b)
        return reflect / n                       # range 0 … 1

    @classmethod
    def scan_iter(cls,
                  seq: Iterable[bytes],
                  score_cutoff: float = 0.75,
                  show_progress: bool = True
                  ) -> List[Tuple[str, str, float]]:
        """
        Iterate over seeds, return [(seed, hex_digest, echo_score)]
        where echo_score ≥ score_cutoff.
        """
        hits, tic = [], time.time()
        for i, item in enumerate(seq, 1):
            raw = item if isinstance(item, (bytes, bytearray)) else str(item).encode()
            d   = cls.sha256(raw)
            s   = cls.echo_score(d)
            if s >= score_cutoff:
                hits.append((raw.decode(errors='replace'), d.hex(), s))
            if show_progress and i % 1_000_000 == 0:
                rate = i / (time.time() - tic + 1e-9)
                print(f"[{i:,}] seeds, {len(hits)} echoes, {rate:,.0f} H/s")
        return hits

# ─────────────────────────────────────────────────────────────────────
def rand_seeds(count: int, byte_len: int = 32) -> Iterable[bytes]:
    for _ in range(count):
        yield os.urandom(byte_len)

def file_seeds(path: Path) -> Iterable[bytes]:
    with path.open("rb") as f:
        for line in f:
            yield line.rstrip(b"\r\n")

# ─────────────────────────────────────────────────────────────────────
def main() -> None:
    ap = argparse.ArgumentParser(description="SHA echo-symmetry scanner (Map-35 helper)")
    g  = ap.add_mutually_exclusive_group(required=False)
    g.add_argument("--random", type=int, metavar="N",
                   help="probe N random 32-byte seeds")
    g.add_argument("--file", type=Path,
                   help="scan one seed per line from text/binary file")
    ap.add_argument("--cutoff", type=float, default=0.75,
                    help="echo-score threshold (default 0.75)")
    args = ap.parse_args()

    if args.random:
        source = rand_seeds(args.random)
        label  = f"{args.random:,} random seeds"
    elif args.file:
        source = file_seeds(args.file)
        label  = f"file: {args.file}"
    else:
        source = (ln.rstrip(b"\r\n") for ln in sys.stdin.buffer)
        label  = "stdin stream"

    print(f"▶ scanning {label} for echo-score ≥ {args.cutoff} …")
    hits = EchoScanner.scan_iter(source, score_cutoff=args.cutoff)

    print("\n— Echo-lock summary —")
    print(f"hits ≥ cutoff        : {len(hits):,}")
    if hits:
        print("top 10 strongest echoes:")
        for seed, hexd, score in sorted(hits, key=lambda t: -t[2])[:10]:
            print(f"  score={score:.2f}  seed={seed[:32]!r:<34}  →  {hexd[:12]}…{hexd[-12:]}")
    else:
        print("no qualifying echoes in this run.")

# ─────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    main()
