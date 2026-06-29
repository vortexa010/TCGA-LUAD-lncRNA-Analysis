# 🧬 Integrated Transcriptomic and Co-expression Network Analysis for the Identification of Candidate Hub Long Non-coding RNAs in Lung Adenocarcinoma (TCGA-LUAD)

![R](https://img.shields.io/badge/R-4.6.0-blue)
![Bioconductor](https://img.shields.io/badge/Bioconductor-DESeq2-green)
![WGCNA](https://img.shields.io/badge/WGCNA-Network%20Analysis-orange)
![Cancer](https://img.shields.io/badge/Disease-LUAD-red)
![Status](https://img.shields.io/badge/Status-Ongoing-success)
![License](https://img.shields.io/badge/License-MIT-yellow)

---

# Project Overview

Lung adenocarcinoma (LUAD) is the most common subtype of non-small cell lung cancer (NSCLC) and remains one of the leading causes of cancer-related mortality worldwide. Long non-coding RNAs (lncRNAs) have emerged as important regulators of gene expression and are increasingly recognized as potential biomarkers and therapeutic targets in cancer.

This repository presents an integrated computational transcriptomic workflow for identifying biologically relevant candidate hub lncRNAs associated with LUAD using RNA sequencing data obtained from The Cancer Genome Atlas (TCGA).

The project integrates differential expression analysis, weighted gene co-expression network analysis (WGCNA), functional enrichment analysis, and external validation to prioritize candidate lncRNAs for future biological and clinical investigation.

---

# 📌 Current Project Status

**Status:** 🟢 Active Research Project

This project is currently being expanded toward a research publication.

## ✅ Completed Analyses

- TCGA-LUAD RNA-seq data acquisition from Genomic Data Commons (GDC)
- RNA-seq count matrix construction
- Data preprocessing and quality control
- Differential expression analysis using DESeq2
- Principal Component Analysis (PCA)
- Volcano plot visualization
- Heatmap visualization
- Differentially expressed lncRNA identification
- High-confidence lncRNA filtering
- Weighted Gene Co-expression Network Analysis (WGCNA)
- Disease-associated module identification
- Hub lncRNA prioritization
- Selection of candidate hub lncRNAs
- Gene Ontology (GO) enrichment analysis
- KEGG pathway enrichment analysis
- External validation using GEPIA2

## 🔄 Ongoing Analyses

- Independent GEO dataset validation
- Gene Set Enrichment Analysis (GSEA)
- Immune infiltration analysis
- RNA-binding protein (RBP) interaction analysis
- Functional characterization
- Experimental validation
- Research manuscript preparation

---

# Objectives

- Download and preprocess TCGA-LUAD RNA sequencing data.
- Construct a unified RNA-seq count matrix.
- Perform differential expression analysis.
- Identify differentially expressed lncRNAs.
- Construct weighted co-expression networks.
- Identify disease-associated modules.
- Prioritize candidate hub lncRNAs.
- Perform functional enrichment analysis.
- Validate candidate lncRNAs using GEPIA2.
- Build a computational framework for future experimental validation.

---

# Dataset

**Database:** The Cancer Genome Atlas (TCGA)

**Project:** TCGA-LUAD

**Portal:** https://portal.gdc.cancer.gov/

### Dataset Summary

| Sample Type | Number |
|-------------|-------:|
| Tumor | 542 |
| Normal | 59 |
| Total | 601 |

### Data Type

- RNA-Seq Gene Expression Quantification
- STAR Counts
- Clinical Metadata

---

# Computational Workflow

```text
TCGA-LUAD RNA-seq
        │
        ▼
Data Download
        │
        ▼
Quality Control & Preprocessing
        │
        ▼
Count Matrix Construction
        │
        ▼
Differential Expression Analysis (DESeq2)
        │
        ▼
Differentially Expressed lncRNA Identification
        │
        ▼
High-confidence lncRNA Filtering
        │
        ▼
Weighted Gene Co-expression Network Analysis (WGCNA)
        │
        ▼
Disease-associated Module Identification
        │
        ▼
Hub lncRNA Prioritization
        │
        ▼
Selection of Candidate Hub lncRNAs
        │
        ▼
GO Enrichment Analysis
        │
        ▼
KEGG Pathway Analysis
        │
        ▼
GEPIA2 Validation
        │
        ▼
Future Downstream Analyses
```

---

# Software

- R (v4.6.0)
- RStudio
- DESeq2
- WGCNA
- clusterProfiler
- org.Hs.eg.db
- enrichplot
- ggplot2
- pheatmap
- dplyr
- tidyverse

---

# Repository Structure

```
TCGA-LUAD-lncRNA
│
├── data/
│   └── README_data_source.txt
│
├── scripts/
│   └── TCGA_LUAD_pipeline.R
│
├── results/
│   ├── PCA.png
│   ├── Volcano.png
│   ├── Heatmap.png
│   ├── WGCNA_SampleClustering.png
│   ├── SoftThreshold.png
│   ├── GO_Enrichment.png
│   ├── KEGG_Enrichment.png
│   └── Candidate_Hub_lncRNAs.csv
│
├── figures/
│
├── README.md
│
└── LICENSE
```

---

# Key Results

### Differential Expression

- Total genes analyzed: **60,616**
- Significant DEGs: **27,660**
- Differentially expressed lncRNAs: **7,965**
- High-confidence lncRNAs: **1,924**

### WGCNA

- Disease-associated blue module identified
- Blue module contained **202 hub lncRNAs**
- Candidate hub lncRNAs prioritized using intramodular connectivity

### Final Candidate Hub lncRNAs

- FENDRR
- LINC00968
- LINC02016
- HID1-AS1
- LANCL1-AS1
- AL606469.1
- AC135012.3
- AL590226.1
- LINC01996
- AC093110.1

---

# Functional Interpretation

Gene Ontology (GO) and KEGG pathway analyses were performed to investigate the biological processes and signaling pathways associated with the prioritized candidate hub lncRNAs. These analyses provide preliminary functional insight into the potential roles of the identified candidates in lung adenocarcinoma.

---

# Future Directions

The current computational framework provides the foundation for additional downstream analyses, including:

- Independent validation using GEO datasets
- Gene Set Enrichment Analysis (GSEA)
- Immune infiltration analysis
- RNA-binding protein (RBP) interaction analysis
- Experimental validation
- Functional characterization
- Biomarker evaluation
- Research publication

---

# Data Availability

The RNA sequencing data used in this study were obtained from the **Genomic Data Commons (GDC)**.

https://portal.gdc.cancer.gov/

Raw TCGA data are **not redistributed** in this repository. Users should download the data directly from the GDC portal in accordance with the GDC data usage policy.

---

# Citation

If you use this repository, please cite:

**Vaishnavi Mangam**

*Integrated Transcriptomic and Co-expression Network Analysis for the Identification of Candidate Hub Long Non-coding RNAs in Lung Adenocarcinoma (TCGA-LUAD).*

GitHub Repository.

---

# Author

**Vaishnavi Mangam**

M.Sc. Bioinformatics  
Department of Bioinformatics  
Savitribai Phule Pune University (SPPU), Pune

**Internship Research Project**

National Centre for Cell Science (NCCS), Pune


# Acknowledgements

This work was carried out during an internship at the **National Centre for Cell Science (NCCS), Pune**, under the guidance of **Dr. Vidisha Tripathi**. Publicly available RNA sequencing and clinical data were obtained from **The Cancer Genome Atlas (TCGA)** through the **Genomic Data Commons (GDC)**. The authors gratefully acknowledge the developers of the open-source R packages and Bioconductor tools used throughout this study.

---

# License

This project is released under the **MIT License**.
