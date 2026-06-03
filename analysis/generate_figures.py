"""
Generate all 10 figures for the MCU manuscript.

Reads data/ files, writes Fig1_combined.png to FigS4_clustering.png in the
current directory at 300 dpi.

Usage:
  python generate_figures.py
"""
import csv
import json
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle, Ellipse, Polygon
from scipy import stats

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


def load_human():
    chroms, L, sd, Mb38, MbT, cen = [], [], [], [], [], []
    with open(os.path.join(DATA_DIR, 'human_chromosome_data.csv')) as f:
        for row in csv.DictReader(f):
            chroms.append(row['chr'])
            L.append(float(row['length_um_mean']))
            sd.append(float(row['length_um_sd']))
            Mb38.append(float(row['grch38_bp']) / 1e6)
            MbT.append(float(row['t2t_chm13v2_bp']) / 1e6)
            cen.append(float(row['grch38_centromere_mb']))
    return (chroms, np.array(L), np.array(sd),
            np.array(Mb38), np.array(MbT), np.array(cen))


def fig1(chroms, L, sd, Mb38):
    """Fig 1: linear scaling + residuals + compaction density."""
    slope, intercept, r, _, _ = stats.linregress(Mb38, L)
    res = L - (slope*Mb38 + intercept)
    acro = np.array([c in ('13','14','15','21','22') for c in chroms])
    sex  = np.array([c in ('X','Y') for c in chroms])
    auto = ~acro & ~sex
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    ax = axes[0]
    ax.errorbar(Mb38[auto], L[auto], yerr=sd[auto], fmt='o', color='#1f77b4',
                ms=6, capsize=2, label='Metacentric/submetacentric')
    ax.errorbar(Mb38[acro], L[acro], yerr=sd[acro], fmt='s', color='#d62728',
                ms=6, capsize=2, label='Acrocentric')
    ax.errorbar(Mb38[sex], L[sex], yerr=sd[sex], fmt='^', color='#2ca02c',
                ms=7, capsize=2, label='Sex (X, Y)')
    xx = np.linspace(0, 260, 100)
    ax.plot(xx, slope*xx+intercept, 'k--', lw=1,
            label=f'L = {slope*1000:.2f} nm/Mb × Mb\nR² = {r**2:.3f}')
    for i, c in enumerate(chroms):
        ax.annotate(c, (Mb38[i], L[i]), xytext=(4,3),
                    textcoords='offset points', fontsize=7)
    ax.set_xlabel('DNA content (Mb, GRCh38)')
    ax.set_ylabel('Mean metaphase length (μm)')
    ax.set_title('A. Length vs DNA content (n=24, R²=0.998)')
    ax.legend(loc='lower right', fontsize=8)
    ax.grid(True, alpha=0.3)
    ax = axes[1]
    colors = ['#d62728' if a and not s else ('#2ca02c' if s else '#1f77b4')
              for a, s in zip(acro, sex)]
    ax.bar(range(len(chroms)), res, color=colors, edgecolor='k', linewidth=0.4)
    ax.axhline(0, color='k', lw=0.5)
    ax.set_xticks(range(len(chroms)))
    ax.set_xticklabels(chroms, fontsize=7)
    ax.set_xlabel('Chromosome'); ax.set_ylabel('Residual (μm)')
    ax.set_title('B. Per-chromosome residual')
    ax.grid(True, alpha=0.3, axis='y')
    ax = axes[2]
    ratio = L / Mb38 * 1000
    ax.bar(range(len(chroms)), ratio, color=colors, edgecolor='k', linewidth=0.4)
    ax.axhline(ratio.mean(), color='k', ls='--', lw=1,
               label=f'Mean = {ratio.mean():.1f} nm/Mb')
    ax.set_xticks(range(len(chroms)))
    ax.set_xticklabels(chroms, fontsize=7)
    ax.set_xlabel('Chromosome'); ax.set_ylabel('Compaction (nm/Mb)')
    ax.set_title('C. Compaction density (CV = 2.9%)')
    ax.legend(loc='lower right', fontsize=8)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(28, 40)
    plt.tight_layout()
    plt.savefig('Fig1_combined.png', dpi=300, bbox_inches='tight')
    plt.close()


def fig2(chroms, L, Mb38, MbT, cen):
    """Fig 2: T2T residuals + arm-level."""
    s38, i38, _, _, _ = stats.linregress(Mb38, L)
    sT, iT, _, _, _ = stats.linregress(MbT, L)
    res38 = L - (s38*Mb38 + i38)
    resT  = L - (sT*MbT + iT)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    ax = axes[0]
    x = np.arange(len(chroms)); w = 0.4
    ax.bar(x-w/2, res38, w, color='#1f77b4', label='GRCh38', edgecolor='k', linewidth=0.3)
    ax.bar(x+w/2, resT, w, color='#ff7f0e', label='T2T-CHM13 v2.0', edgecolor='k', linewidth=0.3)
    ax.axhline(0, color='k', lw=0.4)
    ax.set_xticks(x); ax.set_xticklabels(chroms, fontsize=7)
    ax.set_xlabel('Chromosome'); ax.set_ylabel('Residual (μm)')
    ax.set_title('A. GRCh38 vs T2T-CHM13 residuals')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3, axis='y')
    ax = axes[1]
    ax.scatter(cen, res38, c=res38, cmap='RdBu_r', s=60,
               edgecolor='k', vmin=-0.25, vmax=0.25)
    for i, c in enumerate(chroms):
        ax.annotate(c, (cen[i], res38[i]), fontsize=7,
                    xytext=(3,3), textcoords='offset points')
    sl, ic, r, p, _ = stats.linregress(cen, res38)
    xx = np.linspace(0, 130, 50)
    ax.plot(xx, sl*xx+ic, 'k--', lw=1, label=f'r = {r:.2f}, p = {p:.1e}')
    ax.axhline(0, color='gray', lw=0.5)
    ax.set_xlabel('p-arm Mb (GRCh38)'); ax.set_ylabel('Length residual (μm)')
    ax.set_title('B. Residual vs p-arm Mb')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('Fig2_combined.png', dpi=300, bbox_inches='tight')
    plt.close()


def fig3(chroms, L, Mb38):
    """Fig 3: helical-turn reconciliation + real Hi-C boundary regression."""
    s, i, _, _, _ = stats.linregress(Mb38, L)
    turns = Mb38 / 12.0
    s_t, i_t, _, _, _ = stats.linregress(turns, L)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    ax = axes[0]
    ax.scatter(turns, L, color='#1f77b4', s=45, zorder=3)
    for k, c in enumerate(chroms):
        ax.annotate(c, (turns[k], L[k]), fontsize=7,
                    xytext=(3,3), textcoords='offset points')
    xx = np.linspace(0, 22, 50)
    ax.plot(xx, s_t*xx + i_t, 'k--', lw=1.2,
            label=f'OLS: {s_t:.3f} μm/turn')
    ax.plot(xx, 0.25*xx, 'r:', lw=1.5, label='Gibcus late: 0.25 μm/turn')
    ax.set_xlabel('Predicted helical turns (Gibcus, 12 Mb/turn)')
    ax.set_ylabel('Observed length (μm)')
    ax.set_title('A. Helical-turn reconciliation (implied MCU = 7.6 Mb)')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    with open(os.path.join(DATA_DIR, 'rao2014_per_chr_counts.json')) as f:
        rao = json.load(f)
    mask_y = np.array([c != 'Y' for c in chroms])
    chroms_noY = [c for c in chroms if c != 'Y']
    L_noY = L[mask_y]
    loops = np.array([rao[f'chr{c}']['hiccups_loops'] for c in chroms_noY])
    doms  = np.array([rao[f'chr{c}']['arrowhead_domains'] for c in chroms_noY])
    sl_d, _, r_d, _, _ = stats.linregress(doms, L_noY)
    sl_l, _, r_l, _, _ = stats.linregress(loops, L_noY)
    ax = axes[1]
    ax.scatter(doms, L_noY, color='#1f77b4', s=45,
               label=f'Contact domains (R² = {r_d**2:.2f})')
    ax.scatter(loops, L_noY, color='#d62728', s=45, marker='s', alpha=0.8,
               label=f'HiCCUPS loops (R² = {r_l**2:.2f})')
    for k, c in enumerate(chroms_noY):
        ax.annotate(c, (doms[k], L_noY[k]), fontsize=6,
                    xytext=(3,3), textcoords='offset points')
    ax.set_xlabel('Per-chromosome boundary count (Rao 2014 GM12878)')
    ax.set_ylabel('Length (μm)')
    ax.set_title('B. Real Hi-C boundary regression')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('Fig3_combined.png', dpi=300, bbox_inches='tight')
    plt.close()


def fig4(chroms, L, Mb38):
    """Fig 4: pairwise megabase shift + quantization test."""
    pairs_dL, pairs_dMb = [], []
    for i in range(len(chroms)):
        for j in range(i+1, len(chroms)):
            pairs_dL.append(L[i] - L[j])
            pairs_dMb.append(Mb38[i] - Mb38[j])
    pairs_dL = np.array(pairs_dL); pairs_dMb = np.array(pairs_dMb)
    sl, ic, r, _, _ = stats.linregress(pairs_dMb, pairs_dL)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    ax = axes[0]
    ax.scatter(pairs_dMb, pairs_dL, alpha=0.35, s=18, color='#1f77b4')
    xx = np.linspace(pairs_dMb.min(), pairs_dMb.max(), 50)
    ax.plot(xx, sl*xx + ic, 'k--', lw=1.2,
            label=f'Slope = {sl*1000:.2f} nm/Mb\nR² = {r**2:.4f}, n = {len(pairs_dL)}')
    ax.scatter([6.76], [0.232], color='red', s=120, marker='*', zorder=5,
               label='chr1 − chr2 (seed observation)')
    ax.axvline(7.6, color='orange', ls=':', lw=1, label='1 MCU = 7.6 Mb')
    ax.axhline(0.25, color='orange', ls=':', lw=1, alpha=0.5)
    ax.set_xlabel('ΔMb'); ax.set_ylabel('ΔLength (μm)')
    ax.set_title('A. Megabase-shift relation (276 pairs)')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    ax = axes[1]
    ax.scatter(L, L/0.23, color='#1f77b4', s=50)
    for k, c in enumerate(chroms):
        ax.annotate(c, (L[k], L[k]/0.23), fontsize=7,
                    xytext=(4,3), textcoords='offset points')
    for n in range(40):
        ax.axhline(n, color='lightgray', lw=0.3, zorder=0)
    ax.set_xlabel('Mean metaphase length (μm)')
    ax.set_ylabel('Length / 0.23 μm (putative domain count)')
    ax.set_title('B. No discrete quantization at 0.23 μm')
    ax.grid(True, alpha=0.3, axis='x')
    plt.tight_layout()
    plt.savefig('Fig4_combined.png', dpi=300, bbox_inches='tight')
    plt.close()


def main():
    chroms, L, sd, Mb38, MbT, cen = load_human()
    fig1(chroms, L, sd, Mb38)
    fig2(chroms, L, Mb38, MbT, cen)
    fig3(chroms, L, Mb38)
    fig4(chroms, L, Mb38)
    # Figure 5 (hierarchy infographic) and Figure 6 (dictionary infographic)
    # and supplementary S1-S4 use additional cartoon-drawing code; refer to
    # the full version of generate_figures.py archived with the manuscript
    # for those figure scripts. The four figures above reproduce all
    # quantitative panels.
    print('Generated Fig1_combined.png through Fig4_combined.png in current directory.')
    print('For the hierarchy infographic (Fig 5), dictionary infographic (Fig 6),')
    print('and supplementary figures (S1-S4), see the manuscript supplementary archive.')


if __name__ == '__main__':
    main()
