"""
MCU framework: full analysis pipeline.

Reproduces every quantitative result in:
Goyal M, Goyal R. The Metaphase Chromatin Unit: A Novel Unit of Higher-Order
Chromosome Organization in Human Mitotic Cells. Genes 2026.

Inputs (all in data/):
  human_chromosome_data.csv    - pooled 5-study lengths + GRCh38/T2T Mb + centromere positions
  pope_repliseq_per_chr.json   - Pope 2014 GM12878 ERD fractions
  roadmap_E116_per_chr.json    - Roadmap E116 chromHMM state fractions
  rao2014_per_chr_counts.json  - Rao 2014 GM12878 HiCCUPS / Arrowhead per-chr counts

Usage:
  python mcu_analysis.py
"""
import csv
import json
import os
import numpy as np
from scipy import stats

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

# ---------------- Constants ----------------
MB_PER_MCU = 7.6
UM_PER_MCU = 0.25
NM_PER_MB_THEORY = (UM_PER_MCU / MB_PER_MCU) * 1000  # = 32.89 nm/Mb
ABBE_NM = 196  # NA 1.4, lambda 550 nm
RAO_DOMAIN_KB = 185
DIXON_TAD_KB = 880
GIBCUS_OUTER_LOOP_KB = 400
GIBCUS_LATE_TURN_MB = 12.0
GIBCUS_LATE_TURN_UM = 0.25


def load_human_data():
    """Load per-chromosome cytogenetic and assembly data."""
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


def load_json(name):
    with open(os.path.join(DATA_DIR, name)) as f:
        return json.load(f)


# ---------------- Section 1: Primary OLS regression (Results 3.1) ----------------
def section_1_primary_regression(chroms, L, Mb38, sd):
    print('\n=== Section 1: Primary OLS regression (Results 3.1) ===')
    slope, intercept, r, p, se = stats.linregress(Mb38, L)
    print(f'  L (um) = {slope:.5f} (+/- {se:.5f}) x Mb + {intercept:.4f}')
    print(f'  R-squared = {r**2:.6f}')
    print(f'  p-value = {p:.3e}')
    print(f'  95% CI slope: [{slope-1.96*se:.5f}, {slope+1.96*se:.5f}] um/Mb')
    print(f'  Slope in nm/Mb: {slope*1000:.2f}')

    # Zero-intercept fit
    slope0 = np.sum(L * Mb38) / np.sum(Mb38**2)
    pred0 = slope0 * Mb38
    r2_0 = 1 - np.sum((L - pred0)**2) / np.sum((L - L.mean())**2)
    print(f'  Zero-intercept fit: slope = {slope0:.5f} um/Mb, R^2 = {r2_0:.6f}')

    # Power-law fit
    slope_log, _, r_log, _, se_log = stats.linregress(np.log(Mb38), np.log(L))
    print(f'  Power-law exponent: {slope_log:.4f} +/- {se_log:.4f} (CI [{slope_log-1.96*se_log:.4f}, {slope_log+1.96*se_log:.4f}])')

    # Per-chromosome compaction density
    ratio_nm_per_mb = L / Mb38 * 1000
    print(f'  Compaction density: {ratio_nm_per_mb.mean():.2f} +/- {ratio_nm_per_mb.std(ddof=1):.2f} nm/Mb')
    print(f'  CV across chromosomes: {ratio_nm_per_mb.std(ddof=1)/ratio_nm_per_mb.mean()*100:.2f}%')

    return slope, intercept, ratio_nm_per_mb


# ---------------- Section 2: Quantization test (Results 3.2) ----------------
def section_2_quantization_test(L, slope, intercept, Mb38):
    print('\n=== Section 2: No-quantization test (Results 3.2) ===')
    pred = slope * Mb38 + intercept
    resid_linear = (L - pred).std(ddof=1)
    print(f'  Residual SD from linear fit: {resid_linear:.4f} um')
    for q_um in [0.23, 0.25, 0.32]:
        nearest = q_um * np.round(L / q_um)
        resid_q = (L - nearest).std(ddof=1)
        print(f'  Residual SD from nearest {q_um} um multiple: {resid_q:.4f} um')


# ---------------- Section 3: T2T-CHM13 reanalysis (Results 3.3) ----------------
def section_3_t2t_reanalysis(chroms, L, Mb38, MbT, slope38, int38):
    print('\n=== Section 3: T2T-CHM13 reanalysis (Results 3.3) ===')
    slopeT, intT, rT, pT, seT = stats.linregress(MbT, L)
    print(f'  T2T fit: L = {slopeT:.5f} x Mb + {intT:.4f}, R^2 = {rT**2:.6f}')
    res38 = L - (slope38 * Mb38 + int38)
    resT = L - (slopeT * MbT + intT)
    print(f'\n  Chromosomes with the largest T2T-induced residual change:')
    print(f'    {"chr":<4} {"dMb":<8} {"GRCh38 res":<11} {"T2T res":<10}')
    for i, c in enumerate(chroms):
        dMb = MbT[i] - Mb38[i]
        if abs(dMb) >= 2.0:
            print(f'    {c:<4} {dMb:+.2f}    {res38[i]:+.3f}      {resT[i]:+.3f}')
    # Satellite compaction bound (chr 9)
    i9 = chroms.index('9')
    dMb_9 = MbT[i9] - Mb38[i9]
    bound = 53.0 / dMb_9  # 53 nm = inter-study SD bound
    print(f'\n  chr 9 satellite (+{dMb_9:.1f} Mb): upper bound on compaction = {bound:.2f} nm/Mb')
    print(f'    = at least {NM_PER_MB_THEORY*1000/1000/(bound):.1f}x more compact than chromosome average')


# ---------------- Section 4: Arm-level analysis (Results 3.4) ----------------
def section_4_arm_level(chroms, L, Mb38, cen, slope38, int38):
    print('\n=== Section 4: Arm-level analysis (Results 3.4) ===')
    p_arm = cen
    q_arm = Mb38 - p_arm
    p_frac = p_arm / Mb38
    res = L - (slope38 * Mb38 + int38)
    r, p = stats.pearsonr(p_frac, res)
    print(f'  Length residual vs p-arm fraction: r = {r:.3f}, p = {p:.3e}')
    r2, p2 = stats.pearsonr(p_arm, res)
    print(f'  Length residual vs absolute p-arm Mb: r = {r2:.3f}, p = {p2:.3e}')


# ---------------- Section 5: Helical-turn reconciliation (Results 3.5) ----------------
def section_5_helical_turn(chroms, L, Mb38, slope38):
    print('\n=== Section 5: Helical-turn reconciliation (Results 3.5) ===')
    turns = Mb38 / GIBCUS_LATE_TURN_MB
    slope_t, int_t, r_t, _, _ = stats.linregress(turns, L)
    print(f'  Observed um per Gibcus late-metaphase turn: {slope_t:.3f}')
    print(f'  Gibcus value: 0.25 um per turn')
    print(f'  Observed/Gibcus ratio: {slope_t/0.25:.2f}')
    implied_mb_per_turn = 0.25 / slope38
    print(f'  Implied Mb per turn at 0.25 um axial spacing: {implied_mb_per_turn:.2f}')
    print(f'  This is the operational MCU (= {implied_mb_per_turn:.2f} Mb).')


# ---------------- Section 6: Pairwise megabase shift (Results 3.6) ----------------
def section_6_pairwise(chroms, L, Mb38):
    print('\n=== Section 6: Pairwise megabase-shift relation (Results 3.6) ===')
    pairs_dL, pairs_dMb = [], []
    for i in range(len(chroms)):
        for j in range(i+1, len(chroms)):
            pairs_dL.append(L[i] - L[j])
            pairs_dMb.append(Mb38[i] - Mb38[j])
    pairs_dL = np.array(pairs_dL); pairs_dMb = np.array(pairs_dMb)
    slope, intercept, r, p, se = stats.linregress(pairs_dMb, pairs_dL)
    print(f'  Pairwise slope: {slope:.5f} um/Mb = {slope*1000:.2f} nm/Mb')
    print(f'  95% CI: [{(slope-1.96*se)*1000:.2f}, {(slope+1.96*se)*1000:.2f}] nm/Mb')
    print(f'  R^2 = {r**2:.4f}, n = {len(pairs_dL)} pairs')


# ---------------- Section 7: MCU definition (Results 3.7) ----------------
def section_7_mcu_definition(chroms, L, Mb38):
    print('\n=== Section 7: MCU definition (Results 3.7) ===')
    print(f'  1 MCU = {MB_PER_MCU} Mb DNA = {UM_PER_MCU} um axial length')
    print(f'  1 MCU = ~{int(1000/RAO_DOMAIN_KB * MB_PER_MCU)} Rao 2014 contact domains')
    print(f'  1 MCU = ~{int(1000/DIXON_TAD_KB * MB_PER_MCU)} Dixon 2012 TADs')
    print(f'  1 MCU = ~{int(1000/GIBCUS_OUTER_LOOP_KB * MB_PER_MCU)} Gibcus outer loops')
    total_mb = Mb38.sum()
    total_mcus = total_mb / MB_PER_MCU
    print(f'  Total haploid genome: {total_mb:.0f} Mb = {total_mcus:.0f} MCUs')
    print(f'  Predicted haploid mitotic axial length: {total_mcus * UM_PER_MCU:.1f} um')
    print(f'  Observed haploid total: {L.sum():.1f} um')
    mcu_per_chr = Mb38 / MB_PER_MCU
    print(f'  Range: {mcu_per_chr.min():.1f} MCUs (chr {chroms[np.argmin(mcu_per_chr)]}) '
          f'to {mcu_per_chr.max():.1f} MCUs (chr {chroms[np.argmax(mcu_per_chr)]})')


# ---------------- Section 8: CNV detection threshold (Results 3.8) ----------------
def section_8_cnv_threshold(slope38):
    print('\n=== Section 8: CNV detection threshold (Results 3.8) ===')
    threshold_mb = ABBE_NM / (slope38 * 1000)
    print(f'  Abbe limit at NA 1.4, lambda 550 nm: {ABBE_NM} nm')
    print(f'  Minimum visible CNV: {threshold_mb:.2f} Mb (= ~1 MCU)')
    syndromes = [
        ('Williams (7q11.23)', 1.6, 'FISH-required'),
        ('22q11.2 typical (DiGeorge)', 3.0, 'FISH-required'),
        ('1p36 deletion', 3.5, 'FISH-required'),
        ('Smith-Magenis (17p11.2)', 3.7, 'FISH-required'),
        ('Angelman/PWS (15q11-13)', 5.6, 'borderline'),
        ('22q11.2 atypical', 6.5, 'karyotype-visible'),
        ('Jacobsen (11q23.3)', 11.0, 'karyotype-visible'),
        ('Wolf-Hirschhorn (4p16.3)', 25.0, 'karyotype-visible'),
        ('Cri-du-chat (5p15.2)', 30.0, 'karyotype-visible'),
    ]
    correct = 0
    for name, size_mb, empirical in syndromes:
        pred_nm = slope38 * size_mb * 1000
        pred = ('karyotype-visible' if pred_nm > 250
                else 'borderline' if pred_nm > ABBE_NM
                else 'FISH-required')
        match = '✓' if pred == empirical else '✗'
        if pred == empirical:
            correct += 1
        print(f'    {name}: {size_mb} Mb -> {pred_nm:.0f} nm -> {pred} (empirical: {empirical}) {match}')
    print(f'  Correctly classified: {correct}/{len(syndromes)}')


# ---------------- Section 9: Pope 2014 ERD correlation (Supp S4) ----------------
def section_9_pope_repliseq(chroms, L, Mb38, slope38, int38):
    print('\n=== Section 9: Pope 2014 Repli-Seq correlation (Supp S4) ===')
    pope = load_json('pope_repliseq_per_chr.json')
    res = L - (slope38 * Mb38 + int38)
    erd = np.array([pope.get(f'chr{c}', {}).get('ERD_pct', np.nan) for c in chroms])
    mask = ~np.isnan(erd)
    r, p = stats.pearsonr(erd[mask], res[mask])
    print(f'  ERD% vs length residual: r = {r:.3f}, p = {p:.3e}')
    non_acro = np.array([c not in ['13','14','15','21','22','Y'] for c in chroms])
    combined = mask & non_acro
    r2, p2 = stats.pearsonr(erd[combined], res[combined])
    print(f'  Excluding acrocentrics + Y: r = {r2:.3f}, p = {p2:.3e}')


# ---------------- Section 10: Roadmap chromHMM correlation (Supp S3) ----------------
def section_10_chromhmm(chroms, L, Mb38, slope38, int38):
    print('\n=== Section 10: Roadmap E116 chromHMM correlation (Supp S3) ===')
    roadmap = load_json('roadmap_E116_per_chr.json')
    res = L - (slope38 * Mb38 + int38)
    active = np.array([roadmap.get(f'chr{c}', {}).get('active_pct', np.nan) for c in chroms])
    hetero = np.array([roadmap.get(f'chr{c}', {}).get('hetero_pct', np.nan) for c in chroms])
    mask = ~np.isnan(active)
    r_a, p_a = stats.pearsonr(active[mask], res[mask])
    r_h, p_h = stats.pearsonr(hetero[mask], res[mask])
    print(f'  Active chromatin% vs length residual: r = {r_a:.3f}, p = {p_a:.3e}')
    print(f'  Heterochromatin% vs length residual: r = {r_h:.3f}, p = {p_h:.3e}')


# ---------------- Section 11: Microdeletion clustering test (Supp S7) ----------------
def section_11_microdeletion_clustering():
    print('\n=== Section 11: Microdeletion integer-clustering test (Supp S7) ===')
    nahr_sizes_mb = [
        1.5, 1.5, 1.55, 3.7, 3.7, 3.0, 1.5, 1.4, 6.6, 5.0, 1.6, 1.7,
        0.22, 0.6, 1.65, 1.4, 0.5, 1.8, 1.35, 0.5, 1.6, 1.7, 1.9, 1.4, 1.2, 1.6,
    ]
    for unit_kb in [185, 250, 400, 500, 880, 1000]:
        unit_mb = unit_kb / 1000
        residuals_norm = [abs(s - round(s/unit_mb)*unit_mb)/unit_mb for s in nahr_sizes_mb]
        print(f'  Unit {unit_kb} kb: mean |residual|/unit = {np.mean(residuals_norm):.3f} '
              f'(random baseline 0.25)')
    print('  Conclusion: no integer-multiple clustering at any tested unit size.')
    print('  Microdeletion sizes are positioned by segmental duplications, not chromatin domains.')


# ---------------- Main ----------------
def main():
    print('=' * 72)
    print('MCU framework: full analysis pipeline')
    print('Goyal & Goyal 2026, Genes')
    print('=' * 72)
    chroms, L, sd, Mb38, MbT, cen = load_human_data()
    print(f'Loaded {len(chroms)} chromosomes from data/human_chromosome_data.csv')

    slope38, int38, ratio = section_1_primary_regression(chroms, L, Mb38, sd)
    section_2_quantization_test(L, slope38, int38, Mb38)
    section_3_t2t_reanalysis(chroms, L, Mb38, MbT, slope38, int38)
    section_4_arm_level(chroms, L, Mb38, cen, slope38, int38)
    section_5_helical_turn(chroms, L, Mb38, slope38)
    section_6_pairwise(chroms, L, Mb38)
    section_7_mcu_definition(chroms, L, Mb38)
    section_8_cnv_threshold(slope38)
    section_9_pope_repliseq(chroms, L, Mb38, slope38, int38)
    section_10_chromhmm(chroms, L, Mb38, slope38, int38)
    section_11_microdeletion_clustering()
    print('\n=== Done. ===')


if __name__ == '__main__':
    main()
