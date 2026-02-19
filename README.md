# Q47-Prime-Constellation-Census

### Statistical Morphology and Geodesic Rigidity of Prime Constellations in Q(n) = n⁴⁷ − (n−1)⁴⁷

> Complete census of the first 2 billion cases · 18,473,571 primes · **15 quadruplets** · Modular-3 Exclusion Principle

---

## Abstract

This repository contains the data, algorithms, and paper source for an exhaustive morphological census of prime-generating integers for Q(n) = n⁴⁷ − (n−1)⁴⁷ over the range n ∈ \[1, 2 × 10⁹\]. The survey identifies **15 prime quadruplets** (four consecutive integers n, n+1, n+2, n+3 all generating probable primes of 341–430 digits), correcting the previously published count of 14. A new quadruplet at **n = 23,159,557** was discovered when the \[10⁶, 10⁸\] data gap was filled.

Key theoretical results:
- **Modular-3 Exclusion Principle** (Lemma 3.1): Q(n) ≡ 1 (mod 3) for all n, so Q(n)+k is composite whenever k ≡ 2 (mod 3). This provably kills the right twin-prime slot and establishes intrinsic chirality.
- **Geodesic Rigidity**: The pair-to-solitary ratio R₂(x) decays only ~7% over \[0.5B, 2B\], far slower than the ~11% predicted by random models.
- **Conditional Neighborhood Sieve**: A formally defined search space that exploits chirality and conditional density for deep-space surveys (n > 10¹¹).

## Repository Structure

```
Q47-Prime-Constellation-Census/
├── paper/
│   ├── Q47_Morphology_Paper.tex          # Full LaTeX source
│   ├── Q47_Morphology_Paper.pdf          # Compiled paper (11 pages)
│   ├── figure1_hierarchy.pdf             # Fig 1: Morphological hierarchy
│   ├── figure2_rigidity.pdf              # Fig 2: Geodesic rigidity
│   └── figure3_discovery_map.pdf         # Fig 3: 15 quadruplets map
├── data/
│   └── Q47_morphology_COMPLETE.csv       # Per-10⁸ bin statistics
├── scripts/
│   ├── titan_miner_multicore.py          # Exhaustive constellation sweeper
│   ├── titan_radar_ultimate_5000.py      # Chirality-aware deep-space radar
│   └── generate_figures.py               # Reproduce all paper figures
└── README.md
```

## Core Modules

### 1. `titan_miner_multicore.py` — The Sweeper

Exhaustive search for prime constellations in Q(n).

- **Sieve shield**: Pre-computed Sieve of Eratosthenes up to 5 × 10⁶ eliminates composites by trial division
- **Primality test**: `gmpy2.is_strong_prp()` multi-round Miller–Rabin (error probability < 2⁻¹²⁸)
- **Morphology classifier**: Detects maximal consecutive chains (solitary / pair / triplet / quadruplet)
- **Multi-core**: Parallelized over arbitrary range partitions

### 2. `titan_radar_ultimate_5000.py` — The Deep Space Radar

Targeted search around identified quadruplet seeds, implementing the Conditional Neighborhood Sieve.

- **Chirality enforcement**: Applies the Modular-3 Exclusion Principle — skips all offsets k ≡ 2 (mod 3)
- **Left-spectrum scan**: Searches Q(n) − k for k ∈ \[2, R\], k ≢ 2 (mod 3)
- **Conditional density exploitation**: Prioritizes quadruplet neighborhoods where the modular pre-filter guarantees elevated local primality

### 3. `generate_figures.py` — Figure Reproduction

Standalone script to reproduce all three publication figures from `Q47_morphology_COMPLETE.csv`.

```bash
pip install matplotlib numpy
python scripts/generate_figures.py
```

## Key Results

| Morphology | Count | Frequency (per 10⁹) | Percentage |
|:---|---:|---:|---:|
| Solitary (k=1) | 18,121,562 | 9,060,781 | 98.095% |
| Pair (k=2) | 173,351 | 86,676 | 0.938% |
| Triplet (k=3) | 1,749 | 875 | 0.009% |
| **Quadruplet (k=4)** | **15** | **7.5** | **0.00008%** |

### The 15 Quadruplets

| # | Starting n | Digits | n / 10⁹ |
|---:|---:|---:|---:|
| **1★** | **23,159,557** | **341** | **0.023** |
| 2 | 117,309,848 | 380 | 0.117 |
| 3 | 136,584,738 | 385 | 0.137 |
| 4 | 218,787,064 | 390 | 0.219 |
| 5 | 411,784,485 | 400 | 0.412 |
| 6 | 423,600,750 | 401 | 0.424 |
| 7 | 523,331,634 | 405 | 0.523 |
| 8 | 640,399,031 | 408 | 0.640 |
| 9 | 987,980,498 | 415 | 0.988 |
| 10 | 1,163,461,515 | 420 | 1.163 |
| 11 | 1,370,439,187 | 423 | 1.370 |
| 12 | 1,643,105,964 | 426 | 1.643 |
| 13 | 1,691,581,855 | 427 | 1.692 |
| 14 | 1,975,860,550 | 429 | 1.976 |
| 15 | 1,996,430,175 | 430 | 1.996 |

★ New discovery (this work). Previously absent from all published catalogs.

## Dependencies

```bash
pip install gmpy2 numpy matplotlib
```

## Citation

If you use this data or code, please cite:

```bibtex
@article{chen2026q47,
  title   = {Statistical Morphology and Geodesic Rigidity of Prime
             Constellations in {$Q(n) = n^{47} - (n-1)^{47}$}:
             A Complete Census of the First 2 Billion Cases},
  author  = {Chen, Ruqing},
  year    = {2026},
  note    = {Titan Project, GUT Geoservice Inc.},
  url     = {https://github.com/Ruqing1963/Q47-Prime-Constellation-Census}
}
```

## License

MIT License. See `LICENSE` for details.

---

*Titan Project — GUT Geoservice Inc., Montréal, Canada*
