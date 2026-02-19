#!/usr/bin/env python3
"""
generate_figures.py — Reproduce all publication figures for:

  "Statistical Morphology and Geodesic Rigidity of Prime Constellations
   in Q(n) = n^{47} - (n-1)^{47}: A Complete Census of the First
   2 Billion Cases"

  Author: Ruqing Chen, GUT Geoservice Inc., Montreal
  Repository: https://github.com/Ruqing1963/Q47-Prime-Constellation-Census

Usage:
    python generate_figures.py            # outputs to current directory
    python generate_figures.py --outdir paper/   # outputs to paper/

Outputs (PDF vector + PNG 600 DPI):
    figure1_hierarchy.pdf / .png      — Morphological hierarchy (log-scale)
    figure2_discovery_map.pdf / .png  — Spatial distribution of 15 quadruplets
    figure3_rigidity.pdf / .png       — Pair-to-Solitary ratio R₂(x)

Dependencies: numpy, matplotlib
"""

import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import matplotlib.patheffects as pe
import warnings
warnings.filterwarnings('ignore')

# ═════════════════════════════════════════════════════════════════════
# GLOBAL STYLE
# ═════════════════════════════════════════════════════════════════════
plt.rcParams.update({
    'font.family':      'serif',
    'font.serif':       ['DejaVu Serif', 'Times New Roman', 'Computer Modern Roman'],
    'mathtext.fontset': 'cm',
    'font.size':        10.5,
    'axes.labelsize':   12,
    'axes.titlesize':   11.5,
    'xtick.labelsize':  9.5,
    'ytick.labelsize':  9.5,
    'legend.fontsize':  9,
    'figure.dpi':       150,
    'savefig.dpi':      600,
    'savefig.bbox':     'tight',
    'savefig.pad_inches': 0.05,
    'axes.linewidth':   0.65,
    'axes.grid':        True,
    'axes.axisbelow':   True,
    'xtick.direction':  'in',
    'ytick.direction':  'in',
    'grid.alpha':       0.22,
    'grid.linewidth':   0.4,
    'lines.linewidth':  1.3,
})

# ═════════════════════════════════════════════════════════════════════
# DATA — From Q47_morphology_COMPLETE.csv (verified, gap-free)
# ═════════════════════════════════════════════════════════════════════

BIN_MIDS = np.arange(0, 2e9, 1e8) + 5e7

SOLITARY = np.array([
    1059531, 978760, 951065, 938075, 924841, 917482, 910318,
    902568, 898472, 893721, 888837, 884084, 881161, 878894,
    875527, 871338, 870327, 866375, 865520, 864666])

PAIRS = np.array([
    11790, 10107, 9848, 9311, 8893, 8761, 8679, 8616, 8421,
    8349, 8424, 8072, 8206, 8082, 8092, 8095, 8007, 7806,
    7914, 7878])

TRIPLETS = np.array([
    144, 129, 101, 92, 75, 89, 88, 85, 80, 75,
    70, 71, 76, 87, 74, 79, 83, 69, 107, 75])

# All 15 quadruplet starting positions (complete corrected catalog)
QUAD_STARTS = np.array([
      23_159_557,   117_309_848,  136_584_738,  218_787_064,
     411_784_485,   423_600_750,  523_331_634,  640_399_031,
     987_980_498, 1_163_461_515, 1_370_439_187, 1_643_105_964,
   1_691_581_855, 1_975_860_550, 1_996_430_175])

# Cumulative sums for rigidity analysis
CUM_X  = np.array([i * 1e8 for i in range(1, 21)])
CUM_P1 = np.cumsum(SOLITARY)
CUM_P2 = np.cumsum(PAIRS)
CUM_R2 = CUM_P2 / CUM_P1

# ═════════════════════════════════════════════════════════════════════
# PALETTE
# ═════════════════════════════════════════════════════════════════════
C_SOL   = '#6A7B94'   # Slate blue-grey
C_PAIR  = '#0077B6'   # Deep ocean blue
C_TRIP  = '#D45E00'   # Burnt sienna
C_QUAD  = '#B5121B'   # Deep crimson
C_NEW   = '#7B2D8E'   # Royal violet — new discovery
C_RIGID = '#1A7A3A'   # Forest green
C_HL    = '#AA3344'   # Hardy–Littlewood baseline


def fmt_n(x, _):
    """Axis formatter: 0→'0', 5e8→'500M', 1e9→'1.0B', etc."""
    if x == 0:       return '0'
    if x >= 1e9:     return f'{x/1e9:.1f}B'
    if x >= 1e8:     return f'{x/1e8:.0f}00M'
    if x >= 1e6:     return f'{x/1e6:.0f}M'
    return f'{x:.0f}'


# ═══════════════════════════════════════════════════════════════════
def make_figure1(outdir):
    """Figure 1 — Morphological hierarchy (log-scale decay)."""
    print('▸ Figure 1 — Morphological Hierarchy …')

    fig, ax = plt.subplots(figsize=(7.2, 4.8))

    ax.semilogy(BIN_MIDS, SOLITARY, '-', color=C_SOL, lw=2.4, zorder=3,
                label=r'$\pi_1\;$ Solitary  ($k\!=\!1$)')
    ax.semilogy(BIN_MIDS, PAIRS, '-o', color=C_PAIR, lw=1.7, ms=3.8,
                markerfacecolor='white', markeredgewidth=0.65, zorder=4,
                label=r'$\pi_2\;$ Pairs  ($k\!=\!2$)')
    trips_safe = np.maximum(TRIPLETS, 0.7)
    ax.semilogy(BIN_MIDS, trips_safe, '-s', color=C_TRIP, lw=1.4, ms=3.8,
                markerfacecolor='white', markeredgewidth=0.65, zorder=5,
                label=r'$\pi_3\;$ Triplets  ($k\!=\!3$)')

    # All 15 quadruplets individually at y=1
    for i, qs in enumerate(QUAD_STARTS):
        is_new = (i == 0)
        c  = C_NEW if is_new else C_QUAD
        ms = 15   if is_new else 12
        mec = 'white' if is_new else '#3A0000'
        label = (r'$\pi_4\;$ Quadruplets  ($k\!=\!4$): 15 total' if i == 1 else
                 r'$\pi_4\;$ #1$^{\bigstar}$ New Discovery' if is_new else None)
        ax.semilogy(qs, 1, marker='*', ms=ms, color=c, zorder=12 if is_new else 10,
                    markeredgecolor=mec, markeredgewidth=0.5 if is_new else 0.4,
                    linestyle='none', label=label)

    # Hierarchy arrows
    ax.annotate('', xy=(1.6e9, SOLITARY[16]), xytext=(1.6e9, PAIRS[16]),
                arrowprops=dict(arrowstyle='<->', color='#555', lw=0.7, shrinkA=3, shrinkB=3))
    ax.text(1.65e9, np.sqrt(float(SOLITARY[16]) * float(PAIRS[16])),
            r'$\times\,104.5$', fontsize=8.5, color=C_PAIR, fontweight='bold', va='center', ha='left')
    ax.annotate('', xy=(1.75e9, PAIRS[17]), xytext=(1.75e9, trips_safe[17]),
                arrowprops=dict(arrowstyle='<->', color='#555', lw=0.7, shrinkA=3, shrinkB=3))
    ax.text(1.80e9, np.sqrt(float(PAIRS[17]) * float(trips_safe[17])),
            r'$\times\,99.1$', fontsize=8.5, color=C_TRIP, fontweight='bold', va='center', ha='left')

    ax.text(0.50, 0.95, r'"Noise" — solitary background  ($\sim\!10^6$ per bin)',
            transform=ax.transAxes, fontsize=8.5, color=C_SOL, fontstyle='italic', va='top',
            path_effects=[pe.withStroke(linewidth=3, foreground='white')])
    ax.text(0.50, 0.12, r'"Signal" — 15 rare quadruplets across entire range',
            transform=ax.transAxes, fontsize=8.5, color=C_QUAD, fontstyle='italic', va='bottom',
            path_effects=[pe.withStroke(linewidth=3, foreground='white')])
    ax.annotate(r'#1$^{\bigstar}$ New', xy=(QUAD_STARTS[0], 1),
                xytext=(QUAD_STARTS[0]+6e7, 4), fontsize=7.5, color=C_NEW, fontweight='bold',
                arrowprops=dict(arrowstyle='->', color=C_NEW, lw=0.8),
                path_effects=[pe.withStroke(linewidth=2, foreground='white')])

    ax.set_xlabel(r'$n$', fontsize=12.5)
    ax.set_ylabel(r'Count per $10^8$ bin  (quadruplets: individual events)', fontsize=11)
    ax.set_title(r'Figure 1.  Morphological hierarchy of $Q(n)=n^{47}-(n{-}1)^{47}$'
                 r' prime constellations', fontsize=11, fontweight='bold', pad=12)
    ax.xaxis.set_major_formatter(FuncFormatter(fmt_n))
    ax.set_xlim(-5e7, 2.15e9); ax.set_ylim(0.5, 3.5e6)
    leg = ax.legend(loc='center right', framealpha=0.94, edgecolor='#bbb',
                    fancybox=False, borderpad=0.7, handlelength=2.2)
    leg.get_frame().set_linewidth(0.5)

    fig.tight_layout()
    fig.savefig(f'{outdir}/figure1_hierarchy.pdf')
    fig.savefig(f'{outdir}/figure1_hierarchy.png')
    plt.close(fig)
    print('  ✓ figure1_hierarchy.pdf / .png')


# ═══════════════════════════════════════════════════════════════════
def make_figure2(outdir):
    """Figure 2 — Discovery Map (spatial distribution of 15 quadruplets)."""
    print('▸ Figure 2 — Discovery Map …')

    fig, (ax_main, ax_rug, ax_gap) = plt.subplots(
        3, 1, figsize=(7.4, 5.0),
        height_ratios=[4.5, 1.2, 1.3],
        gridspec_kw={'hspace': 0.06})

    # Top panel
    sol_norm = SOLITARY / SOLITARY.max()
    ax_main.fill_between(BIN_MIDS, 0, sol_norm, color=C_SOL, alpha=0.10, zorder=1)
    ax_main.plot(BIN_MIDS, sol_norm, '-', color=C_SOL, lw=0.6, alpha=0.4, zorder=2)

    for i, qs in enumerate(QUAD_STARTS):
        is_new = (i == 0)
        c  = C_NEW if is_new else C_QUAD
        lw = 2.0  if is_new else 0.9
        ms = 14   if is_new else 10
        ax_main.axvline(qs, color=c, lw=lw, alpha=0.85 if is_new else 0.5, zorder=4+is_new)
        y_pos = 0.88 - 0.09 * (i % 3)
        ax_main.plot(qs, y_pos, '*', ms=ms, color=c, zorder=8,
                     markeredgecolor='white' if is_new else '#4a0000',
                     markeredgewidth=0.6 if is_new else 0.35)
        lbl = r'#1$^{\bigstar}$' if is_new else f'#{i+1}'
        ax_main.text(qs, y_pos+0.065, lbl, ha='center', fontsize=6.5,
                     color=c, fontweight='bold' if is_new else 'normal',
                     path_effects=[pe.withStroke(linewidth=2.5, foreground='white')])

    ax_main.annotate('New Discovery\n' r'$n = 23\,159\,557$' '\n'
                     r'$Q(n) \approx 341$ digits',
                     xy=(QUAD_STARTS[0], 0.88), xytext=(3.8e8, 0.45),
                     fontsize=9, color=C_NEW, fontweight='bold',
                     bbox=dict(boxstyle='round,pad=0.5', fc='#F5ECFF', ec=C_NEW, lw=1.3),
                     arrowprops=dict(arrowstyle='-|>', color=C_NEW, lw=1.6,
                                     connectionstyle='arc3,rad=0.2'))
    ax_main.set_xlim(-4e7, 2.08e9); ax_main.set_ylim(0, 1.15)
    ax_main.set_ylabel('Normalized density', fontsize=10)
    ax_main.set_title(r'Figure 2.  Spatial distribution of 15 prime quadruplets, '
                      r'$n \leq 2 \times 10^{9}$', fontsize=11, fontweight='bold', pad=12)
    ax_main.xaxis.set_major_formatter(FuncFormatter(fmt_n))
    ax_main.tick_params(axis='x', labelbottom=False)
    ax_main.set_yticks([0, 0.5, 1.0])

    # Rug panel
    for i, qs in enumerate(QUAD_STARTS):
        ax_rug.axvline(qs, color=C_NEW if i == 0 else C_QUAD,
                       lw=3.5 if i == 0 else 1.6, zorder=5)
    for k in range(1, 15):
        ax_rug.axvline(k * 2e9/15, color='#ccc', lw=0.45, ls=':', zorder=1)
    ax_rug.set_xlim(-4e7, 2.08e9); ax_rug.set_ylim(0, 1); ax_rug.set_yticks([])
    ax_rug.set_ylabel('Rug', fontsize=7.5, color='#888', rotation=0, labelpad=18, va='center')
    ax_rug.xaxis.set_major_formatter(FuncFormatter(fmt_n))
    ax_rug.tick_params(axis='x', labelbottom=False)
    ax_rug.text(0.005, 0.5, 'Dense', transform=ax_rug.transAxes,
                fontsize=7, color=C_QUAD, va='center', fontstyle='italic')
    ax_rug.text(0.96, 0.5, 'Sparse', transform=ax_rug.transAxes,
                fontsize=7, color=C_QUAD, va='center', fontstyle='italic', ha='right')

    # Gap panel
    gaps = np.diff(QUAD_STARTS) / 1e6
    gap_mids = (QUAD_STARTS[:-1] + QUAD_STARTS[1:]) / 2
    ax_gap.bar(gap_mids, gaps, width=5.5e7, color=C_QUAD, alpha=0.55,
               edgecolor='#600', linewidth=0.3, zorder=3)
    ax_gap.bar(gap_mids[0], gaps[0], width=5.5e7, color=C_NEW, alpha=0.7,
               edgecolor=C_NEW, linewidth=0.6, zorder=4)
    ax_gap.axhline(gaps.mean(), color='#333', ls='--', lw=0.7, alpha=0.5, zorder=2)
    ax_gap.text(2.0e9, gaps.mean()+8, f'mean = {gaps.mean():.0f}M',
                fontsize=7.5, color='#555', ha='right', va='bottom')
    ax_gap.set_xlim(-4e7, 2.08e9)
    ax_gap.set_ylabel(r'$\Delta n$ (M)', fontsize=9)
    ax_gap.set_xlabel(r'$n$', fontsize=12.5)
    ax_gap.xaxis.set_major_formatter(FuncFormatter(fmt_n))

    fig.tight_layout()
    fig.savefig(f'{outdir}/figure2_discovery_map.pdf')
    fig.savefig(f'{outdir}/figure2_discovery_map.png')
    plt.close(fig)
    print('  ✓ figure2_discovery_map.pdf / .png')


# ═══════════════════════════════════════════════════════════════════
def make_figure3(outdir):
    """Figure 3 — Rigidity Verification (R₂(x) ratio evolution)."""
    print('▸ Figure 3 — Rigidity Verification …')

    fig, ax = plt.subplots(figsize=(7.2, 4.2))

    R_bin = PAIRS / SOLITARY
    ax.plot(BIN_MIDS, R_bin*1000, 'o', ms=5.5, color=C_PAIR, alpha=0.55, zorder=4,
            markeredgecolor='#003050', markeredgewidth=0.4,
            label=r'Per-bin  $R_2 = \pi_2 / \pi_1$')
    ax.plot(CUM_X, CUM_R2*1000, '-', color=C_RIGID, lw=2.5, zorder=6,
            label=r'Cumulative  $R_2(x)$')

    coeffs = np.polyfit(CUM_X, CUM_R2*1000, 1)
    slope_per_B = coeffs[0] * 1e9
    ax.plot(CUM_X, np.polyval(coeffs, CUM_X), '--', color='#2B2B2B', lw=1.0, alpha=0.65,
            zorder=5, label=f'Linear trend  (slope = {slope_per_B:.4f} / $10^9$)')

    HL_x = np.linspace(1e8, 2e9, 300)
    HL_y = CUM_R2[0]*1000 * np.log(CUM_X[0]) / np.log(HL_x)
    ax.plot(HL_x, HL_y, ':', color=C_HL, lw=1.1, alpha=0.55, zorder=3,
            label=r'Random baseline  $\propto 1/\ln(n)$')

    ax.axvspan(5e8, 2.05e9, color=C_RIGID, alpha=0.045, zorder=0)
    ax.axvline(5e8, color=C_RIGID, lw=0.7, ls='--', alpha=0.4, zorder=2)

    pct_drop = (1 - CUM_R2[19]*1000 / (CUM_R2[4]*1000)) * 100
    sol_pct  = (1 - SOLITARY[-1] / SOLITARY[4]) * 100
    ax.annotate(f'Geometric Rigidity Zone\n'
                f'$R_2$ decays only {pct_drop:.1f}% over [0.5B, 2B]\n'
                f'vs. {sol_pct:.1f}% decay for per-bin $\\pi_1$',
                xy=(1.25e9, np.interp(1.25e9, CUM_X, CUM_R2*1000)),
                xytext=(2.5e8, 9.05), fontsize=8.5, color=C_RIGID,
                fontstyle='italic', fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.45', fc='white', ec=C_RIGID, lw=1.0, alpha=0.92),
                arrowprops=dict(arrowstyle='->', color=C_RIGID, lw=1.0,
                                connectionstyle='arc3,rad=-0.18'))

    ax.set_xlabel(r'$n$', fontsize=12.5)
    ax.set_ylabel(r'$R_2 \;\;(\times\, 10^{3})$', fontsize=12)
    ax.set_title(r'Figure 3.  Pair-to-Solitary ratio $R_2(x)$:  '
                 'geodesic rigidity verification', fontsize=11, fontweight='bold', pad=12)
    ax.xaxis.set_major_formatter(FuncFormatter(fmt_n))
    ax.set_xlim(-5e7, 2.15e9); ax.set_ylim(8.5, 12.0)
    leg = ax.legend(loc='upper right', framealpha=0.94, edgecolor='#bbb',
                    fancybox=False, borderpad=0.7, fontsize=8.5)
    leg.get_frame().set_linewidth(0.5)

    fig.tight_layout()
    fig.savefig(f'{outdir}/figure3_rigidity.pdf')
    fig.savefig(f'{outdir}/figure3_rigidity.png')
    plt.close(fig)
    print('  ✓ figure3_rigidity.pdf / .png')


# ═══════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate paper figures')
    parser.add_argument('--outdir', default='.', help='Output directory')
    args = parser.parse_args()

    make_figure1(args.outdir)
    make_figure2(args.outdir)
    make_figure3(args.outdir)

    print(f'\n{"═"*55}')
    print(f'  All 3 figures saved to {args.outdir}/')
    print(f'{"═"*55}')
