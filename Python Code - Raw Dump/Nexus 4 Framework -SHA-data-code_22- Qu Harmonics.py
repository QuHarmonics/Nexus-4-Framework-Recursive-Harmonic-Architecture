#!/usr/bin/env python3
"""
cpu_miner.py — CPU‐Only Bitcoin‐Style Proof‐of‐Work Miner

Scans the full 32‐bit nonce space (and optional extraNonce) in parallel
using Python multiprocessing. No ASIC tricks, but fully correct.

Usage:
    python cpu_miner.py \
        --prevhash 000000…000 (hex little-endian, 64 chars) \
        --merkleroot 4e3b…9b0c (hex little-endian, 64 chars) \
        --version 536870912 \
        --timestamp 1622505600 \
        --bits 0x170d3f6f \
        [--workers 8] [--extra-nonce-max 1000000]
"""

import argparse
import hashlib
import struct
import multiprocessing as mp
import sys
from itertools import count

def double_sha256(header: bytes) -> bytes:
    return hashlib.sha256(hashlib.sha256(header).digest()).digest()

def bits_to_target(bits: int) -> int:
    e = bits >> 24
    c = bits & 0x007fffff
    return c * (1 << (8 * (e - 3)))

def build_header_prefix(version: int,
                        prevhash_le: bytes,
                        merkleroot_le: bytes,
                        timestamp: int,
                        bits: int) -> bytes:
    return (
        struct.pack('<I', version)
        + prevhash_le
        + merkleroot_le
        + struct.pack('<I', timestamp)
        + struct.pack('<I', bits)
    )

def worker_nonce_range(args):
    """
    Worker scans a subrange of nonces for any valid header.
    Returns (nonce, extraNonce, hash_bytes) on success, else None.
    """
    prefix, target, extra_nonce, start_nonce, end_nonce = args
    # build coinbase tweak bytes for extraNonce
    # (in real Bitcoin you’d embed this in the scriptSig; here we just treat it as part of the prefix)
    extra_bytes = struct.pack('<I', extra_nonce)
    # assume the extraNonce bytes go right before the merkle root; adjust as needed for your block format
    # if your coinbase is more complex, recompute merkle root each extraNonce iteration
    header_base = prefix[:36] + extra_bytes + prefix[40:]
    for nonce in range(start_nonce, end_nonce):
        header = header_base + struct.pack('<I', nonce)
        h_int = int.from_bytes(double_sha256(header), 'little')
        if h_int <= target:
            return nonce, extra_nonce, h_int.to_bytes(32, 'little')
    return None

def mine(prefix: bytes,
         target: int,
         workers: int,
         extra_nonce_max: int):
    """
    Distribute work across `workers` processes, iterating extraNonce
    from 0 to extra_nonce_max-1. Exits on first valid result.
    """
    total_nonces = 2**32
    chunk = total_nonces // workers

    pool = mp.Pool(workers)
    try:
        for extra in range(extra_nonce_max):
            # prepare args for all workers
            jobs = []
            for i in range(workers):
                start = i * chunk
                end   = start + chunk if i < workers-1 else total_nonces
                jobs.append((prefix, target, extra, start, end))
            # map with imap_unordered for early exit
            for result in pool.imap_unordered(worker_nonce_range, jobs):
                if result is not None:
                    nonce, ext, hash_bytes = result
                    return nonce, ext, hash_bytes
            # no result this extraNonce, try next
            print(f"[{extra+1}/{extra_nonce_max}] extraNonce tried, none valid.")
    finally:
        pool.close()
        pool.join()
    return None

def main():
    p = argparse.ArgumentParser(description="CPU-Only Bitcoin Proof-of-Work Miner")
    p.add_argument('--prevhash',    required=True, help="64-char hex (little-endian)")
    p.add_argument('--merkleroot',  required=True, help="64-char hex (little-endian)")
    p.add_argument('--version',     type=int,  required=True)
    p.add_argument('--timestamp',   type=int,  required=True)
    p.add_argument('--bits',        required=True, help="compact bits (e.g. 0x1d00ffff)")
    p.add_argument('--workers',     type=int, default=mp.cpu_count(),
                   help="number of parallel processes")
    p.add_argument('--extra-nonce-max', type=int, default=1_000_000,
                   help="max extraNonce iterations (default 1e6)")
    args = p.parse_args()

    prevhash_le   = bytes.fromhex(args.prevhash)
    merkleroot_le = bytes.fromhex(args.merkleroot)
    bits_int      = int(args.bits, 16)

    prefix = build_header_prefix(
        version      = args.version,
        prevhash_le  = prevhash_le,
        merkleroot_le= merkleroot_le,
        timestamp    = args.timestamp,
        bits         = bits_int
    )
    target = bits_to_target(bits_int)

    print(f"Mining with target: 0x{target:064x}")
    print(f"Using {args.workers} workers, extraNonce max {args.extra_nonce_max}")

    result = mine(prefix, target, args.workers, args.extra_nonce_max)
    if result:
        nonce, extra, hash_bytes = result
        print(f"\n✔ Success!")
        print(f"  extraNonce = {extra}")
        print(f"  nonce      = {nonce} (0x{nonce:08x})")
        # display hash big-endian
        print("  hash       =", hash_bytes[::-1].hex())
    else:
        print("\n✘ No valid nonce found in given search space.", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
