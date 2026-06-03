# Metaphase Chromatin Unit (MCU) Calculator

Interactive tool accompanying:

**Goyal M, Goyal R. The Metaphase Chromatin Unit (MCU): a 7.6 Mb / 0.25 μm unit bridging cytogenetics, Hi-C, and clinical CNV resolution in human metaphase chromosomes.** *Genes* 2026 (submitted).

## What this is

A single-page web tool that translates between five units of chromosome scale:

- DNA content (Mb, kb)
- Metaphase axial length (μm, nm)
- Metaphase Chromatin Units (MCU)
- Hi-C boundary counts (Rao contact domains, Dixon TADs, Gibcus mitotic loops)
- Cytogenetic CNV detection verdict (karyotype-visible, borderline, FISH-required)

Type any value into any field; the other six update live. Reference tables for all 24 human chromosomes, common aneuploidies (T13, T18, T21, XXX, XXY, XYY, X0), and major microdeletion syndromes (Williams, 22q11.2, Smith-Magenis, Angelman/PWS, Jacobsen, Wolf-Hirschhorn, Cri-du-chat) are included.

## Conversion constants

| Quantity | Value | Source |
|---|---|---|
| 1 MCU | 7.6 Mb DNA | this study |
| 1 MCU | 0.25 μm axial length | Gibcus 2018 |
| Linear compaction | 33.4 nm/Mb (95% CI [32.3, 33.5]) | this study, R² = 0.998 |
| Karyotype detection threshold | ~6.0 Mb (Abbe limit at 200 nm) | this study |
| 1 MCU = | ~41 Rao 2014 contact domains | Rao 2014 |
| 1 MCU = | ~9 Dixon 2012 TADs | Dixon 2012 |
| 1 MCU = | ~150 Gibcus mitotic loops | Gibcus 2018 |

## Use it online

Live at: **https://goyal74.github.io/mcu-calculator/**

## Use it offline

Download `index.html`, open in any modern browser. No installation, no dependencies.

## Deployment to GitHub Pages

1. Create a new public GitHub repository (e.g. `mcu-calculator`).
2. Upload the two files in this folder (`index.html` and `README.md`) to the repo root.
3. In repo Settings → Pages, set Source = "Deploy from a branch", Branch = `main`, Folder = `/ (root)`. Save.
4. Wait ~1 minute. The tool will be live at `https://goyal74.github.io/mcu-calculator/`.
5. (Optional) Connect to Zenodo for a permanent DOI:
   - Go to https://zenodo.org/account/settings/github/
   - Toggle the `mcu-calculator` repo on.
   - Create a GitHub release (e.g. v1.0). Zenodo will mint a DOI.
   - Add the DOI badge to this README and to the manuscript's Data Availability statement.

## Citing the calculator

Please cite the manuscript:

> Goyal M, Goyal R. The Metaphase Chromatin Unit (MCU): a 7.6 Mb / 0.25 μm unit bridging cytogenetics, Hi-C, and clinical CNV resolution in human metaphase chromosomes. *Genes* 2026.

If using the tool itself in clinical or research workflows, additionally cite the Zenodo DOI (once minted).

## License

MIT License. Free to use, modify, redistribute. See [LICENSE](LICENSE).

## Contact

Correspondence: Ravi Goyal, Department of Obstetrics and Gynecology, University of Arizona College of Medicine. Email: goyalr@arizona.edu.

## Acknowledgements

The 33 nm/Mb empirical scaling is derived from five pooled cytogenetic studies of human peripheral lymphocytes (Penrose 1964, Cohen 1966, Mayall 1984, Van Dyke 1986, Naumova 2013), the GRCh38 and T2T-CHM13 v2.0 reference assemblies, the Rao 2014 GM12878 Hi-C dataset, the Pope 2014 Repli-Seq dataset, and the Roadmap Epigenomics E116 chromHMM segmentation. See the manuscript for full attribution.
