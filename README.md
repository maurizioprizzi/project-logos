# Project LOGOS: Algorithmic Complexity Analysis of the FOXP2 Gene

## Overview

This project applies **Information Theory** (Shannon Entropy), **Algorithmic Complexity** (Kolmogorov via Zlib), and **statistical inference** to compare the FOXP2 coding sequence (CDS) across four mammalian species: *Homo sapiens*, *Pan troglodytes*, *Mus musculus*, and *Gorilla gorilla*.

FOXP2 encodes a transcription factor critical for speech and language development. Its protein is among the 5% most conserved in mammals — differing by only 2 amino acids between humans and chimpanzees, and 3 between humans and mice.

### Version History

**v2.0 (February 2026) — Revised methodology.** The original analysis (v1) compared non-equivalent mRNA isoforms, producing a spurious 64% size difference and artificial information-theoretic divergence. This version corrects three critical flaws:

1. **Isoform control:** CDS extracted from canonical isoforms via fixed RefSeq accession numbers, eliminating UTR-driven size artifacts.
2. **Size-bias correction:** Non-overlapping windows (200 bp) for compression analysis, eliminating DEFLATE dictionary bias.
3. **Effect size reporting:** Cohen's d alongside p-values, preventing statistical significance inflation from autocorrelated sliding windows.

## Key Findings

When comparing equivalent coding sequences, **all four species are informationally indistinguishable** in FOXP2:

| Comparison | Δ Entropy (bits) | p-value | Cohen's d | Δ Compression | p-value | Cohen's d |
|---|---|---|---|---|---|---|
| Human vs Chimpanzee | −0.0004 | 0.988 | 0.007 | 0.0020 | 0.909 | 0.052 |
| Human vs Mouse | +0.0006 | 0.979 | 0.012 | 0.0005 | 0.976 | 0.013 |
| Human vs Gorilla | −0.0012 | 0.961 | 0.022 | 0.0020 | 0.910 | 0.051 |

All Cohen's d values are below 0.07 (negligible effect). The CDS sizes differ by at most 0.6% across species, compared to the 184% difference in full mRNA — confirming that the original findings were artifacts of isoform selection.

### What the original analysis got wrong

| Metric | Original (v1) | Corrected (v2) |
|---|---|---|
| Sequences compared | mRNA (different isoforms) | CDS (canonical isoforms) |
| Human size | 6,618 bp | 2,148 bp |
| Chimpanzee size | 2,380 bp | 2,151 bp |
| Size difference | 64% | 0.1% |
| Entropy difference | 0.0225 bits (p < 0.000001) | 0.0004 bits (p = 0.988) |
| Windows | Overlapping (step=1) | Non-overlapping (200 bp) |
| Effect size | Not reported | Cohen's d = 0.007 |

## Data Sources

Sequences obtained from NCBI RefSeq with controlled accession numbers:

| Species | Accession | CDS (bp) | Protein (aa) |
|---|---|---|---|
| *Homo sapiens* | NM_014491 | 2,148 | 715 |
| *Pan troglodytes* | NM_001009020 | 2,151 | 716 |
| *Mus musculus* | NM_053242 | 2,145 | 714 |
| *Gorilla gorilla* | XM_063708043 | 2,139 | 712 |

## Pipeline

### Requirements

```bash
git clone https://github.com/maurizioprizzi/project-logos.git
cd project-logos
python3 -m venv .venv
source .venv/bin/activate
pip install biopython matplotlib numpy scipy
```

### Execution

```bash
# Step 1: Download sequences (requires internet)
python3 coletas_dados.py

# Step 2: Basic statistics and comparability validation
python3 analise_basica.py

# Step 3: Zlib compression analysis (non-overlapping windows)
python3 analise_avancada.py

# Step 4: Shannon entropy with effect sizes
python3 calculo_entropia.py
```

Each script is self-contained and produces both terminal output and publication-quality figures (PNG, 300 dpi).

### Output Files

| File | Description |
|---|---|
| `*_foxp2_cds.fasta` | Coding sequences per species |
| `*_foxp2_mrna.fasta` | Full mRNA per species (for reference) |
| `metadados.json` | Collection metadata |
| `compressao_foxp2.png` | Compression analysis figure |
| `entropia_foxp2.png` | Entropy analysis figure |

## Tech Stack

- **Python 3.12**
- **BioPython** — GenBank fetching and CDS extraction
- **SciPy** — Welch's t-test, Mann-Whitney U, Levene's test
- **NumPy / Matplotlib** — Numerical analysis and visualization
- **Zlib** — DEFLATE compression as Kolmogorov complexity proxy

## Lessons Learned

This project serves as a case study in computational bioinformatics methodology:

1. **Data selection determines results.** Comparing non-equivalent isoforms created a 64% size difference that does not exist in biology. Controlled accession numbers are essential.
2. **Statistical significance ≠ practical significance.** With thousands of autocorrelated windows, any trivial difference yields p < 0.000001. Effect size (Cohen's d) is indispensable.
3. **Compression algorithms have size bias.** DEFLATE builds more efficient dictionaries for longer sequences. Fixed-size non-overlapping windows eliminate this confound.

## Author

**Maurizio Prizzi**
*Professor of Optimization & Data Science*

## License

MIT