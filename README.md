# 🧬 Cocaine-induced epigenetic remodeling of mouse sperm

Repository containing the bioinformatic analyses used in the manuscript:

> **Cocaine reshapes the sperm epigenome through coordinated DNA methylation and RNA remodeling in mice**

---

# 📂 Repository contents

This repository contains all custom scripts used to generate the analyses presented in the manuscript.

## Included analyses

- 🧬 RRBS hotspot analysis
- 🎯 Transcription factor motif enrichment
- 🧪 Histone ChIP-seq integration
- 📈 RNA-seq differential expression analysis
- 📊 Figure generation
- 📑 Supplementary analyses

---

# 💾 Sequencing data

The raw sequencing datasets are available in the NCBI Gene Expression Omnibus (GEO).

| Dataset | GEO accession |
|---------|---------------|
| RRBS | **GSE341998** |
| RNA-seq | **GSE341997** |

---

# ⚠ Required external files

Some files are not included in this repository because of their size.

## 🧬 RRBS CX reports

The following scripts require the original Bismark **CX_report** files:

- `TF_analysis.R`
- `TE_enrichment.R`

Download the corresponding

```
*_CX_report.txt.gz
```

files from

**GEO: GSE341998**

and place them in the working directory before running these analyses.

---

## 🧪 Histone ChIP-seq data

The histone analysis requires processed **BigWig** files from:

**GEO: GSE79227**

Place the downloaded files inside

```text
histone data/
```

before running

```
Histone_analysis.R
```

---

# 💻 Software

Analyses were performed using

- R 4.5.1
- Bioconductor
- DESeq2
- GenomicRanges
- rtracklayer
- motifmatchr
- JASPAR2024
- tidyverse
- pheatmap
- EnhancedVolcano

---

# 🔬 Genome assembly

All analyses were performed using

**GRCm39 / mm39**

---

# 📜 Citation

If you use these scripts, please cite:

> Gonzalez *et al.*  
> *Cocaine reshapes the sperm epigenome through coordinated DNA methylation and RNA remodeling in mice*

---

# 📧 Contact

**Betina Gonzalez**

Laboratory of Epigenetics and Functional Genomics

Instituto Tecnológico de Buenos Aires (ITBA)

Argentina

<img width="672" height="480" alt="fractal-large" src="https://github.com/user-attachments/assets/fc21c83d-6d40-4465-851a-c80312ebd654" />
