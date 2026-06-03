# MCU Analysis Code

Python analysis code for:

**Goyal M, Goyal R. The Metaphase Chromatin Unit: A Novel Unit of Higher-Order Chromosome Organization in Human Mitotic Cells.** *Genes* 2026 (submitted).

This folder reproduces every quantitative result, table, and figure in the manuscript and Supplementary Material from the primary inputs.

## Quick start

```bash
# Install dependencies
pip install -r requirements.txt

# Run the full analysis (prints all results to stdout)
python mcu_analysis.py

# Generate all 10 figures (Fig 1-6 main, Fig S1-S4 supplementary)
python generate_figures.py
```

Both scripts read from `data/` and require no internet access at run time. Raw public datasets used (Pope 2014 Repli-Seq, Roadmap E116 chromHMM, Rao 2014 Hi-C) were pre-processed once into `data/*.json` files for reproducibility; the original processing commands are documented at the bottom of `mcu_analysis.py`.

## Contents

| File | Purpose |
|---|---|
| `mcu_analysis.py` | Single consolidated analysis script with 11 numbered sections matching the manuscript Results. Prints all reported numbers. |
| `generate_figures.py` | Generates all 10 figures (Fig 1-6 main + Fig S1-S4 supplementary) as 300 dpi PNG. |
| `data/human_chromosome_data.csv` | Per-chromosome pooled cytogenetic length data from 5 studies (Penrose 1964, Cohen 1966, Mayall 1984, Van Dyke 1986, Naumova 2013) and Mb counts from GRCh38 and T2T-CHM13 v2.0. |
| `data/pope_repliseq_per_chr.json` | Per-chromosome Early Replication Domain (ERD) fractions, derived from Pope 2014 GM12878 (GEO GSE53984) segmented BED file. |
| `data/roadmap_E116_per_chr.json` | Per-chromosome chromatin state fractions, derived from Roadmap Epigenomics E116 GM12878 15-state chromHMM. |
| `data/rao2014_per_chr_counts.json` | Per-chromosome HiCCUPS loop and Arrowhead contact-domain counts, derived from Rao 2014 GM12878 (GEO GSE63525). |
| `requirements.txt` | Python dependencies (pinned versions). |

## Reproducing each manuscript result

`mcu_analysis.py` is divided into 11 numbered sections corresponding to the manuscript:

| Section in analysis script | Manuscript section | Key result |
|---|---|---|
| 1 | Results 3.1 | OLS regression L = 0.0329 x Mb + 0.043, R² = 0.998 |
| 2 | Results 3.2 | Quantization test: no improvement over linear fit |
| 3 | Results 3.3 | T2T-CHM13 satellite over-condensation (chr 9, 16, Y) |
| 4 | Results 3.4 | Arm-level: residual vs p-arm Mb, r = 0.65, p = 5.6e-4 |
| 5 | Results 3.5 | Helical-turn reconciliation; implied 7.6 Mb per turn |
| 6 | Results 3.6 | Pairwise (276 pairs) megabase-shift slope = 32.75 nm/Mb |
| 7 | Results 3.7 | MCU definition; 406 MCUs per haploid genome |
| 8 | Results 3.8 | CNV detection threshold; 9-of-9 microdeletion classification |
| 9 | Supp S4 | Pope 2014 ERD correlation with residual |
| 10 | Supp S3 | Roadmap E116 chromHMM correlation |
| 11 | Supp S7 | Integer-domain clustering test (negative result) |

## Data sources and pre-processing

The `data/*.json` files were pre-processed once from the original public datasets. To regenerate them from scratch:

### Pope 2014 Repli-Seq (GM12878)
```bash
# Source: GEO GSE53984
curl -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE53nnn/GSE53984/suppl/GSE53984_GSM923451_Gm12878_Rep1_segments.bed.gz" \
  | gunzip > /tmp/pope.bed
python preprocess_pope.py /tmp/pope.bed > data/pope_repliseq_per_chr.json
```

### Roadmap Epigenomics E116 chromHMM
```bash
# Source: Roadmap Epigenomics Project
curl -L "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E116_15_coreMarks_dense.bed.gz" \
  | gunzip > /tmp/E116.bed
python preprocess_roadmap.py /tmp/E116.bed > data/roadmap_E116_per_chr.json
```

### Rao 2014 GM12878 Hi-C
```bash
# Source: GEO GSE63525
curl -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.txt.gz" \
  | gunzip > /tmp/rao_loops.txt
curl -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt.gz" \
  | gunzip > /tmp/rao_domains.txt
python preprocess_rao.py > data/rao2014_per_chr_counts.json
```

## License

MIT License. See the repo-root LICENSE file.

## Citation

Please cite the manuscript and the Zenodo archive:

- Goyal M, Goyal R. The Metaphase Chromatin Unit: A Novel Unit of Higher-Order Chromosome Organization in Human Mitotic Cells. *Genes* 2026.
- Zenodo DOI: [10.5281/zenodo.20532702](https://doi.org/10.5281/zenodo.20532702)
