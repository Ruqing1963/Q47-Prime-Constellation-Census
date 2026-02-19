#!/usr/bin/env python3
"""
titan_miner_multicore.py — Exhaustive Prime Constellation Sweeper

Surveys Q(n) = n^47 - (n-1)^47 for prime-generating integers and
classifies them into morphological types:
    C1 (Solitary), C2 (Pair), C3 (Triplet), C4 (Quadruplet)

This is the workhorse algorithm used to produce the complete census
for n ∈ [1, 2×10^9] described in:

    R. Chen, "Statistical Morphology and Geodesic Rigidity of Prime
    Constellations in Q(n) = n^47 - (n-1)^47", 2026.

Dependencies: gmpy2, numpy
Install:      pip install gmpy2 numpy

Author: Ruqing Chen, GUT Geoservice Inc., Montreal
"""

import sys
import os
import time
import argparse
import multiprocessing as mp
from functools import partial

try:
    import gmpy2
    from gmpy2 import mpz, is_strong_prp
except ImportError:
    print("ERROR: gmpy2 is required.  Install via: pip install gmpy2")
    sys.exit(1)

import numpy as np

# ═════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═════════════════════════════════════════════════════════════════════

EXPONENT = 47
SIEVE_LIMIT = 5_000_000       # Small-factor trial-division bound
MR_BASES = [2, 3, 5, 7, 11]   # Miller-Rabin witness bases


# ═════════════════════════════════════════════════════════════════════
# SIEVE: Precompute small primes for trial division
# ═════════════════════════════════════════════════════════════════════

def sieve_of_eratosthenes(limit):
    """Return list of primes up to `limit` via Sieve of Eratosthenes."""
    is_prime = bytearray(b'\x01') * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i in range(2, limit + 1) if is_prime[i]]


SMALL_PRIMES = None  # Lazy-initialized


def init_sieve():
    """Initialize the small-prime sieve (called once per process)."""
    global SMALL_PRIMES
    if SMALL_PRIMES is None:
        SMALL_PRIMES = sieve_of_eratosthenes(SIEVE_LIMIT)


# ═════════════════════════════════════════════════════════════════════
# CORE: Q(n) computation and primality test
# ═════════════════════════════════════════════════════════════════════

def Q(n):
    """Compute Q(n) = n^47 - (n-1)^47 using GMP arbitrary precision."""
    n_mpz = mpz(n)
    return n_mpz**EXPONENT - (n_mpz - 1)**EXPONENT


def is_Q_prime(n):
    """
    Test whether Q(n) is a probable prime.

    Stage 1: Trial division by small primes (fast rejection).
    Stage 2: Multi-base strong probable prime test (Miller-Rabin).

    Returns True if Q(n) is a probable prime, False otherwise.
    """
    init_sieve()
    q = Q(n)

    # Stage 1: Trial division
    for p in SMALL_PRIMES:
        if q % p == 0:
            return q == p  # Only prime if Q(n) IS the small prime itself

    # Stage 2: Strong probable prime test (multiple bases)
    for base in MR_BASES:
        if not is_strong_prp(q, base):
            return False

    return True


# ═════════════════════════════════════════════════════════════════════
# SCANNER: Process a range [n_start, n_end)
# ═════════════════════════════════════════════════════════════════════

def scan_range(n_start, n_end, verbose=False):
    """
    Scan Q(n) for n in [n_start, n_end).

    Returns:
        primes: sorted list of n values where Q(n) is probable prime
        stats:  dict with counts {solitary, pairs, triplets, quadruplets}
    """
    init_sieve()
    primes = []

    t0 = time.time()
    for n in range(n_start, n_end):
        if is_Q_prime(n):
            primes.append(n)

        if verbose and (n - n_start) % 100_000 == 0 and n > n_start:
            elapsed = time.time() - t0
            rate = (n - n_start) / elapsed
            print(f"  [{os.getpid()}] n={n:>12,}  "
                  f"primes={len(primes):>8,}  "
                  f"rate={rate:,.0f} n/s", flush=True)

    return primes


def classify_morphology(primes):
    """
    Classify a sorted list of prime-generating n values into
    constellation types.

    Returns dict:
        solitary:    count of isolated primes
        pairs:       count of pair events (2 consecutive)
        triplets:    count of triplet events (3 consecutive)
        quadruplets: count of quadruplet events (4 consecutive)
        quad_starts: list of starting n for each quadruplet
    """
    if not primes:
        return {'solitary': 0, 'pairs': 0, 'triplets': 0,
                'quadruplets': 0, 'quad_starts': []}

    prime_set = set(primes)
    visited = set()
    stats = {'solitary': 0, 'pairs': 0, 'triplets': 0,
             'quadruplets': 0, 'quad_starts': []}

    for n in sorted(primes):
        if n in visited:
            continue

        # Find maximal chain starting at n
        chain_len = 1
        while (n + chain_len) in prime_set:
            chain_len += 1

        # Mark all members as visited
        for i in range(chain_len):
            visited.add(n + i)

        # Classify
        if chain_len == 1:
            stats['solitary'] += 1
        elif chain_len == 2:
            stats['pairs'] += 1
        elif chain_len == 3:
            stats['triplets'] += 1
        elif chain_len >= 4:
            stats['quadruplets'] += 1
            stats['quad_starts'].append(n)
            # If chain > 4, count additional sub-structures
            if chain_len > 4:
                print(f"  *** CHAIN OF LENGTH {chain_len} at n={n}! ***")

    return stats


# ═════════════════════════════════════════════════════════════════════
# WORKER: Multiprocessing wrapper
# ═════════════════════════════════════════════════════════════════════

def worker(args):
    """Worker function for multiprocessing pool."""
    n_start, n_end, verbose = args
    primes = scan_range(n_start, n_end, verbose=verbose)
    return primes


# ═════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='Titan Miner: Exhaustive Q(n) prime constellation search')
    parser.add_argument('start', type=int, help='Start of range (inclusive)')
    parser.add_argument('end', type=int, help='End of range (exclusive)')
    parser.add_argument('--cores', type=int, default=mp.cpu_count(),
                        help=f'Number of CPU cores (default: {mp.cpu_count()})')
    parser.add_argument('--chunk', type=int, default=1_000_000,
                        help='Chunk size per worker (default: 1,000,000)')
    parser.add_argument('--output', type=str, default=None,
                        help='Output file for prime list (default: stdout summary)')
    parser.add_argument('--verbose', action='store_true',
                        help='Print progress per worker')
    args = parser.parse_args()

    n_start, n_end = args.start, args.end
    print(f"{'═'*60}")
    print(f"  Titan Miner v3.0 — Q(n) = n^{EXPONENT} - (n-1)^{EXPONENT}")
    print(f"  Range: [{n_start:,}, {n_end:,})")
    print(f"  Cores: {args.cores}  |  Chunk: {args.chunk:,}")
    print(f"  Sieve: primes up to {SIEVE_LIMIT:,}")
    print(f"{'═'*60}")

    # Build work chunks
    chunks = []
    pos = n_start
    while pos < n_end:
        chunk_end = min(pos + args.chunk, n_end)
        chunks.append((pos, chunk_end, args.verbose))
        pos = chunk_end

    print(f"  {len(chunks)} chunks queued\n")

    # Execute
    t0 = time.time()
    all_primes = []

    if args.cores == 1:
        for chunk in chunks:
            all_primes.extend(worker(chunk))
    else:
        with mp.Pool(args.cores) as pool:
            for result in pool.imap(worker, chunks):
                all_primes.extend(result)

    all_primes.sort()
    elapsed = time.time() - t0

    # Classify
    stats = classify_morphology(all_primes)

    # Report
    print(f"\n{'═'*60}")
    print(f"  RESULTS for [{n_start:,}, {n_end:,})")
    print(f"{'═'*60}")
    print(f"  Total primes:  {len(all_primes):>12,}")
    print(f"  Solitary:      {stats['solitary']:>12,}")
    print(f"  Pairs:         {stats['pairs']:>12,}")
    print(f"  Triplets:      {stats['triplets']:>12,}")
    print(f"  Quadruplets:   {stats['quadruplets']:>12,}")
    if stats['quad_starts']:
        print(f"\n  Quadruplet starting positions:")
        for qs in stats['quad_starts']:
            digits = len(str(Q(qs)))
            print(f"    n = {qs:>15,}  ({digits} digits)")
    print(f"\n  Time: {elapsed:.1f}s  |  Rate: {(n_end-n_start)/elapsed:,.0f} n/s")
    print(f"{'═'*60}")

    # Save
    if args.output:
        with open(args.output, 'w') as f:
            f.write(f"# Titan Miner v3.0 — Q(n) primes\n")
            f.write(f"# Range: [{n_start}, {n_end})\n")
            f.write(f"# Count: {len(all_primes)}\n")
            for n in all_primes:
                f.write(f"{n}\n")
        print(f"  Saved to {args.output}")


if __name__ == '__main__':
    main()
