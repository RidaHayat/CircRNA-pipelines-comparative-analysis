# Comparative Evaluation of Bioinformatics Pipelines for Circular RNA Detection and Reconstruction

[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![University of Verona](https://img.shields.io/badge/University-Verona-red.svg)](https://www.univr.it/)
[![Lab](https://img.shields.io/badge/Lab-InfOmics-teal.svg)](https://github.com/InfOmics)

> **Master's Degree in Medical Bioinformatics (Academic Year 2025–2026)**  
> **Student:** Rida Hayat (Matricola: VR489060)  
> **Supervisor:** Prof. Rosalba Giugno | **Tutor:** Simone Avesani  
> *Department of Computer Science, Università degli Studi di Verona*

---

## 📌 Overview & Research Challenge

Circular RNAs (circRNAs) are covalently closed non-coding RNAs formed via non-canonical back-splicing. Because they lack 5' caps and 3' poly(A) tails, circRNAs are resistant to exonuclease degradation and exceptionally stable in bodily fluids, establishing them as key biomarkers and therapeutic candidates.

However, computational identification from RNA-seq data suffers from **striking discrepancies across algorithms** (due to differing alignment strategies, junction filters, and reconstruction models). Benchmarking reveals that individual tools often report less than 20% overlap, making single-tool discovery prone to false positives or missed candidates.

---

## 🎯 Key Contributions

1. **Systematic 4-Pipeline Benchmark:** A head-to-head evaluation of **CIRCexplorer2**, **CIRI-full**, **CyCLeR**, and **circtools** on high-throughput paired-end RNA-seq data from the **TMD8 diffuse large B-cell lymphoma (DLBCL)** cell line (`TMD8-8h-DMSO`).
2. **Automated R Wrappers:** Developed end-to-end R automation scripts that streamline trimming (Flexbar), alignment (BWA-MEM, STAR), and caller execution.
3. **Harmonized Output Schema:** Standardized all caller outputs into uniform coordinates (`chr:start|end:strand`) and normalized expression matrices (TPM).
4. **EasyCircR Expansion (`EasyCircR-extended`):** Extended the [EasyCircR](https://github.com/InfOmics/EasyCircR) framework to import multi-caller outputs, calculate overlap metrics (Venn, UpSet), and construct high-confidence consensus sets.

---

## 📊 Summary of Benchmark Results

| Pipeline | Detection Strategy | Detections (`TMD8-8h-DMSO`) | Median Expression (log2 TPM+1) | Methodological Strength |
| :--- | :--- | :---: | :---: | :--- |
| **CIRCexplorer2** | Splice-aware chimeric parsing (STAR/TopHat2) | **1,373** | 8.35 | Broad sensitivity; comprehensive alternative splicing annotations. |
| **CIRI-full** | Reverse Overlap (RO) + BSJ integration | **97** | 9.02 | High precision; full-length circular isoform reconstruction. |
| **CyCLeR** | Greedy graph traversal + EM quantification | **149** | 0.89 | Isoform-level resolution; co-quantifies circular and linear transcripts. |
| **circtools** | DCC junction detection + FUCHS reconstruction | **6** | Bimodal (14,080 outlier) | High specificity with strict replicate support filters. |

### Major Findings:
- **Core Consensus:** **108 circRNAs** were consistently detected across the three quantification-capable tools (CIRCexplorer2, CIRI-full, CyCLeR).
- **Quadruple Overlap:** Exactly **1 circRNA** (`chr10:133491216:133498010`) was captured across all four tools, representing the highest-confidence candidate in the dataset.
- **Quantification Discordance:** Up to a **1,618-fold difference** in reported expression across tools, underscoring the critical need for multi-pipeline harmonization before biological interpretation.

### 📈 Benchmark Visualizations

#### Multi-Tool Detection Overlap & Consensus
<p align="center">
  <img src="EasyCircR-extended/img/Venn_3tools.png" width="48%" alt="Triple Overlap Venn Diagram (108 consensus circRNAs)" />
  <img src="EasyCircR-extended/img/UpSet_4tools.png" width="48%" alt="4-Tool Intersection UpSet Plot" />
</p>
<p align="center">
  <em><strong>Figure 1:</strong> Inter-pipeline detection overlap. (Left) 108 consensus circRNAs identified across CIRCexplorer2, CIRI-full, and CyCLeR. (Right) 4-way UpSet intersection revealing exactly 1 circRNA (<code>chr10:133491216:133498010</code>) shared across all four callers.</em>
</p>

#### Quantification Concordance & Discordance
<p align="center">
  <img src="EasyCircR-extended/img/correlation_combined_3tools.png" width="48%" alt="Pairwise Expression Correlation" />
  <img src="EasyCircR-extended/img/Boxplot_global_TPM.png" width="48%" alt="Global circRNA Expression Across Tools (TPM)" />
</p>
<p align="center">
  <em><strong>Figure 2:</strong> Quantitative evaluation. (Left) Pairwise expression correlations (log2(TPM+1)) showing strong correlation between CIRCexplorer2 and CIRI-full (r = 0.73, p = 2.58e-58) and divergence with CyCLeR. (Right) Global expression distributions highlighting severe caller-dependent quantification outliers.</em>
</p>

---

## 🗂️ Repository Structure

```text
CircRNA-pipelines-comparative-analysis/
├── README.md                          # Repository documentation & benchmarking summary
└── EasyCircR-extended/                # Multi-pipeline extension module
    ├── DESCRIPTION                    # R package metadata
    ├── NAMESPACE                      # Exported functions and dependencies
    ├── img/                           # Publication-ready benchmark figures
    │   ├── Venn_3tools.png            # 3-way Venn diagram (108 consensus)
    │   ├── UpSet_4tools.png           # 4-way UpSet intersection plot
    │   ├── correlation_combined_3tools.png # Pairwise correlation scatter plots
    │   ├── Boxplot_global_TPM.png     # Global TPM distribution
    │   └── Spaghetti_triple_all.png   # Multi-tool expression discordance
    └── R/                             # Automated wrappers & parsers
        ├── CIRCexplorer2.R            # CIRCexplorer2 execution driver
        ├── circexplorer2_output.R     # Standardized parser for CIRCexplorer2
        ├── circtools.R                # Flexbar + STAR + circtools driver
        ├── circtools-output.R         # Standardized parser for circtools
        ├── ciri-full.R                # CIRI-full execution driver
        └── ciri-output.R              # Standardized parser for CIRI-full
```

---

## 🚀 Quick Usage

```R
# Source the automated wrappers
source("EasyCircR-extended/R/circtools.R")
source("EasyCircR-extended/R/circtools-output.R")

# Run circtools automated pipeline
run_circtools(
  samples_file = "data/samples.tsv",
  reference_genome = "path/to/hg38.fa",
  annotation_gtf = "path/to/genes.gtf"
)

# Parse output into standardized format
circ_df <- read_circtools_output("results/circtools/circRNA_counts.tsv")
head(circ_df)
```

---

## 📖 Citation

```bibtex
@mastersthesis{hayat2025circrna,
  author  = {Rida Hayat},
  title   = {Comparative Evaluation of Bioinformatics Pipelines for Circular RNA Detection and Reconstruction},
  school  = {Department of Computer Science, Universit{\`a} degli Studi di Verona},
  year    = {2026},
  note    = {Supervised by Prof. Rosalba Giugno and Simone Avesani}
}
```
