# Project LOGOS: Algorithmic Complexity in Human FOXP2 Gene

## 🔬 Overview
This project applies **Information Theory (Shannon Entropy)**, **Algorithmic Complexity (Kolmogorov/Zlib)**, and **Bioinformatics** to analyze the structural differences between the Human (*Homo sapiens*) and Chimpanzee (*Pan troglodytes*) *FOXP2* gene.

While evolutionary biology often focuses on protein-coding similarity (~98%), this project investigates the **informational density** of the mRNA transcripts. As a Professor of Optimization, I apply mathematical modeling to demonstrate that the massive expansion of the human gene (3x larger) exhibits a significant increase in structural order and syntactic density.

## 📊 Key Findings

### 1. Statistical Significance & Bio-Data
* **Length Discrepancy:** The Human transcript (6,618 bp) is **64% larger** than the Chimpanzee (2,380 bp).
* **The "RefSeq" Delta:** Analysis of curated RefSeq data shows a **7.49% drop in GC content** in humans, suggesting a more dynamic and accessible regulatory architecture for gene expression.
* **P-Value Evidence:** Using an Independent T-Test on sliding windows, we found **t = -13.12** and **p < 0.000001**, proving that the difference in information density is not due to chance.

### 2. Shannon Entropy (Order vs. Noise)
![Entropy Analysis](entropia_perfil_janela100_passo20.png)
* **Lower Entropy = Higher Order:** Despite its larger size, the Human gene has a **lower average entropy (1.91 bits)** than the Chimpanzee (1.93 bits).
* **Entropy Valleys:** The human-exclusive 3' UTR region contains "vales" of low entropy, indicating highly structured regulatory "modules" rather than random evolutionary insertions.

### 3. Algorithmic Complexity (Zlib Compression)
![Compression Analysis](comparacao_zlib.png)
* **Compression Ratio:** Human (**0.3049**) vs. Chimpanzee (**0.3189**).
* **Information Density:** The human gene is more compressible, achieving a lower **Bits per Base (2.44)**.
* **Conclusion:** Higher compressibility in a larger sequence is a mathematical hallmark of **modular design and hierarchical syntax**. This suggests that human-specific "non-coding" regions are actually high-level regulatory software.

## 🛠 Tech Stack
* **Language:** Python 3.12
* **Core Libraries:** `BioPython` (GenBank ETL), `SciPy` (Inferential Statistics), `NumPy` & `Matplotlib` (Data Visualization), `Zlib` (Algorithmic Complexity).
* **Environment:** Linux (Ubuntu/Debian), VS Code, Python Virtual Environments (.venv).

## 🚀 How to Run the Pipeline
1.  **Clone & Setup:**
    ```bash
    git clone [https://github.com/maurizioprizzi/project-logos.git](https://github.com/maurizioprizzi/project-logos.git)
    cd project-logos
    python3 -m venv .venv
    source .venv/bin/activate
    pip install biopython matplotlib numpy scipy
    ```
2.  **Execute Analysis:**
    ```bash
    # 1. Fetch data from NCBI RefSeq
    python3 coleta_dados.py

    # 2. Basic EDA (Length, GC%, Dinucleotides)
    python3 analise_basica.py

    # 3. Statistical Analysis (Shannon Entropy & T-Test)
    python3 calculo_entropia.py

    # 4. Complexity Analysis (Zlib/Compression)
    python3 analise_avancada.py

    # 5. Visual Profile Generation
    python3 visualizacao_entropia.py
    ```

## 👨‍🏫 Author
**Maurizio Prizzi**
*Professor of Optimization & Data Science*
*Expert in Mathematical Modeling, Logic, and Complex Systems Analysis*