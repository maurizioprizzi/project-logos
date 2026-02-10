# Project LOGOS: Algorithmic Complexity in Human FOXP2 Gene

## 🔬 Overview
This project applies **Information Theory (Shannon Entropy)**, **Algorithmic Complexity (Kolmogorov/Zlib)**, and **Bioinformatics** to analyze the structural differences between the Human and Chimpanzee *FOXP2* gene (associated with language and cognition).

While evolutionary biology often focuses on protein-coding similarity (~98%), this project investigates the **informational density** of the mRNA transcripts, revealing significant structural divergence in non-coding regulatory regions.

## 📊 Key Findings

### 1. Shannon Entropy (Information Density)
![Entropy Analysis](comparacao_entropia_foxp2.png)
* **Length Discrepancy:** The Human transcript (~6.6kb) is nearly **3x larger** than the Chimpanzee transcript (~2.4kb).
* **Entropy Valleys:** The human-exclusive region (3' UTR) exhibits specific areas of **low entropy (<1.8 bits)**, suggesting high algorithmic structure and regulatory syntax rather than random evolutionary "noise".

### 2. Algorithmic Compressibility (Zlib/DEFLATE)
![Compression Analysis](comparacao_zlib.png)
* **Syntactic Density:** Despite being significantly larger, the Human gene shows a **lower compression ratio (0.30)** compared to the Chimpanzee (0.32).
* **Conclusion:** In Information Theory, higher compressibility implies higher redundancy and order. This suggests the human "junk DNA" regions are actually highly structured regulatory code, not random insertions.

## 🛠 Tech Stack
* **Language:** Python 3.12
* **Libraries:** BioPython, NumPy, Matplotlib, Zlib
* **Data Source:** NCBI GenBank (Dynamic Fetching via Entrez)
* **Concepts:** Shannon Entropy, Sliding Window Analysis, GC Content, Lossless Compression Ratios.

## 🚀 How to Run
1.  Clone the repository:
    ```bash
    git clone [https://github.com/maurizioprizzi/project-logos.git](https://github.com/maurizioprizzi/project-logos.git)
    ```
2.  Install dependencies:
    ```bash
    pip install biopython matplotlib numpy pandas
    ```
3.  Run the analysis pipeline:
    ```bash
    # Step 1: Fetch latest data from NCBI
    python3 coleta_dados.py

    # Step 2: Generate Entropy Graph
    python3 visualizacao_entropia.py

    # Step 3: Calculate Compression Ratios
    python3 analise_avancada.py
    ```

## 👨‍🏫 Author
**Maurizio Prizzi**
*Professor of Optimization & Data Science*
*Expert in Mathematical Modeling & Logic*