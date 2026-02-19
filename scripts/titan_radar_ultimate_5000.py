#!/usr/bin/env python3
"""
titan_radar_ultimate_5000.py — Chirality-Aware Deep Space Radar

Targeted search for satellite primes around known Quadruplet seeds,
implementing the Conditional Neighborhood Sieve (Definition 5.1):

    S = { Q(n) - k  |  n ∈ C₄,  k ∈ [2, R],  k ≢ 2 (mod 3) }

Key principles from the paper:
  1. Modular-3 Exclusion (Lemma 3.1): Q(n) ≡ 1 (mod 3) for all n,
     so Q(n)+k is composite whenever k ≡ 2 (mod 3).
     → Right spectrum (k > 0, twin slot k=+2) is PROVABLY DEAD.
  2. Chirality: Only the LEFT spectrum Q(n)-k is searched.
  3. Conditional Density: Quadruplet neighborhoods have enhanced
     local primality due to favorable modular alignment.

Author: Ruqing Chen, GUT Geoservice Inc., Montreal
Repository: https://github.com/Ruqing1963/Q47-Prime-Constellation-Census

Dependencies: gmpy2
Install:      pip install gmpy2

Usage:
    python titan_radar_ultimate_5000.py                    # scan all 15 seeds
    python titan_radar_ultimate_5000.py --seeds 23159557   # scan specific seed
    python titan_radar_ultimate_5000.py --radius 5000      # search radius
"""

import sys
import time
import argparse

try:
    import gmpy2
    from gmpy2 import mpz, is_strong_prp
except ImportError:
    print("ERROR: gmpy2 is required.  Install via: pip install gmpy2")
    sys.exit(1)

# ═════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═════════════════════════════════════════════════════════════════════

EXPONENT = 47
SIEVE_LIMIT = 5_000_000
MR_BASES = [2, 3, 5, 7, 11, 13]  # Extra base for deep-space confidence

# Complete catalog of 15 quadruplet seeds (n ≤ 2×10⁹)
KNOWN_QUADRUPLETS = [
      23_159_557,   117_309_848,  136_584_738,  218_787_064,
     411_784_485,   423_600_750,  523_331_634,  640_399_031,
     987_980_498, 1_163_461_515, 1_370_439_187, 1_643_105_964,
   1_691_581_855, 1_975_860_550, 1_996_430_175,
]

DEFAULT_RADIUS = 5000  # Maximum offset |k| to search


# ═════════════════════════════════════════════════════════════════════
# SIEVE
# ═════════════════════════════════════════════════════════════════════

_small_primes = None

def get_small_primes():
    global _small_primes
    if _small_primes is None:
        sieve = bytearray(b'\x01') * (SIEVE_LIMIT + 1)
        sieve[0] = sieve[1] = 0
        for i in range(2, int(SIEVE_LIMIT**0.5) + 1):
            if sieve[i]:
                sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
        _small_primes = [i for i in range(2, SIEVE_LIMIT + 1) if sieve[i]]
    return _small_primes


# ═════════════════════════════════════════════════════════════════════
# CORE FUNCTIONS
# ═════════════════════════════════════════════════════════════════════

def Q(n):
    """Compute Q(n) = n^47 - (n-1)^47."""
    n_mpz = mpz(n)
    return n_mpz**EXPONENT - (n_mpz - 1)**EXPONENT


def is_probable_prime(value):
    """
    Test if a large integer is a probable prime.
    Stage 1: Trial division by small primes.
    Stage 2: Multi-base Miller-Rabin.
    """
    if value <= 1:
        return False
    if value <= 3:
        return True

    primes = get_small_primes()
    v = mpz(value)

    for p in primes:
        if v % p == 0:
            return v == p

    for base in MR_BASES:
        if not is_strong_prp(v, base):
            return False

    return True


def is_excluded_by_mod3(k):
    """
    Modular-3 Exclusion Principle (Lemma 3.1):
    Q(n) ≡ 1 (mod 3), so Q(n)-k ≡ (1-k) (mod 3).
    Q(n)-k is divisible by 3 iff k ≡ 1 (mod 3).

    Returns True if offset k is EXCLUDED (Q(n)-k ≡ 0 mod 3).
    """
    return k % 3 == 1


# ═════════════════════════════════════════════════════════════════════
# RADAR SCAN
# ═════════════════════════════════════════════════════════════════════

def radar_scan(seed_n, radius=DEFAULT_RADIUS, verbose=True):
    """
    Scan the left-spectrum neighborhood of a quadruplet seed.

    Implements the search space S (Definition 5.1):
        S = { Q(n) - k  |  k ∈ [2, R],  k ≢ 1 (mod 3) }

    Note: k ≡ 1 (mod 3) means Q(n)-k ≡ 0 (mod 3), so excluded.
          Positive offsets (Q(n)+k) are excluded by chirality.

    Args:
        seed_n:  Starting n of a known quadruplet
        radius:  Maximum offset k to search
        verbose: Print progress

    Returns:
        List of (k, digit_count) for each satellite prime found
    """
    q_base = Q(seed_n)
    q_digits = len(str(q_base))

    if verbose:
        print(f"\n{'─'*60}")
        print(f"  RADAR SCAN: Quadruplet seed n = {seed_n:,}")
        print(f"  Q(n) ≈ {q_digits} digits")
        print(f"  Search space: Q(n) - k,  k ∈ [2, {radius}],  k ≢ 1 (mod 3)")
        print(f"  Right spectrum: BLOCKED (Modular-3 Exclusion Principle)")
        print(f"{'─'*60}")

    satellites = []
    tested = 0
    skipped = 0
    t0 = time.time()

    for k in range(2, radius + 1):
        # Enforce Modular-3 Exclusion
        if is_excluded_by_mod3(k):
            skipped += 1
            continue

        candidate = q_base - k
        tested += 1

        if is_probable_prime(candidate):
            digits = len(str(candidate))
            satellites.append((k, digits))
            if verbose:
                status = "★ LEFT TWIN" if k == 2 else "● SATELLITE"
                print(f"  {status}  k={k:<6}  "
                      f"Q(n)-{k} is PRP  ({digits} digits)")

        if verbose and tested % 500 == 0:
            elapsed = time.time() - t0
            print(f"  ... tested {tested}/{radius - skipped}  "
                  f"found {len(satellites)}  "
                  f"({elapsed:.1f}s)", flush=True)

    elapsed = time.time() - t0
    if verbose:
        total_possible = radius - 1  # k=2..R
        print(f"\n  Summary: {len(satellites)} satellites found")
        print(f"  Tested: {tested}  |  Skipped (mod 3): {skipped}  "
              f"|  Time: {elapsed:.1f}s")
        print(f"  Efficiency: {skipped}/{total_possible} candidates "
              f"pruned by Chirality Lemma ({100*skipped/total_possible:.1f}%)")

    return satellites


# ═════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='Titan Radar: Chirality-aware deep-space satellite search')
    parser.add_argument('--seeds', type=int, nargs='+', default=None,
                        help='Quadruplet seed n values (default: all 15 known)')
    parser.add_argument('--radius', type=int, default=DEFAULT_RADIUS,
                        help=f'Search radius R (default: {DEFAULT_RADIUS})')
    parser.add_argument('--output', type=str, default=None,
                        help='Output file for results')
    parser.add_argument('--quiet', action='store_true',
                        help='Suppress per-candidate output')
    args = parser.parse_args()

    seeds = args.seeds if args.seeds else KNOWN_QUADRUPLETS

    print(f"{'═'*60}")
    print(f"  Titan Deep Space Radar v5000")
    print(f"  Q(n) = n^{EXPONENT} - (n-1)^{EXPONENT}")
    print(f"  Seeds: {len(seeds)} quadruplet(s)")
    print(f"  Radius: {args.radius}")
    print(f"  Strategy: Left-spectrum only (Mod-3 Exclusion)")
    print(f"{'═'*60}")

    all_results = {}
    total_satellites = 0

    for seed in seeds:
        sats = radar_scan(seed, radius=args.radius, verbose=not args.quiet)
        all_results[seed] = sats
        total_satellites += len(sats)

    # Final summary
    print(f"\n{'═'*60}")
    print(f"  GRAND SUMMARY")
    print(f"{'═'*60}")
    for seed, sats in all_results.items():
        twin = any(k == 2 for k, _ in sats)
        twin_str = " ★ HAS LEFT TWIN" if twin else ""
        print(f"  n={seed:>15,}  →  {len(sats):>3} satellites{twin_str}")
    print(f"{'─'*60}")
    print(f"  Total satellites across all seeds: {total_satellites}")
    print(f"{'═'*60}")

    # Save
    if args.output:
        with open(args.output, 'w') as f:
            f.write(f"# Titan Radar v5000 results\n")
            f.write(f"# Radius: {args.radius}\n")
            f.write(f"# seed_n, offset_k, satellite_digits\n")
            for seed, sats in all_results.items():
                for k, d in sats:
                    f.write(f"{seed},{k},{d}\n")
        print(f"  Results saved to {args.output}")


if __name__ == '__main__':
    main()
