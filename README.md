# 🧬 Cocaine-induced epigenetic remodeling of mouse sperm

Repository containing the bioinformatic analyses used in the manuscript:

> **Cocaine reshapes the sperm epigenome through coordinated DNA methylation and RNA remodeling in mice**

---

# 📂 Repository contents

This repository contains the custom scripts and input files used to generate the analyses presented in the manuscript.

## Included analyses

- 🧬 RRBS hotspot analysis
- 🎯 Transcription factor motif analysis
- 🧪 Histone ChIP-seq integration
- 📈 RNA-seq differential expression analysis
- 📊 Figure generation
- 📑 Supplementary analyses

---

# 💾 Sequencing data

The raw sequencing datasets generated in this study are available in the NCBI Gene Expression Omnibus (GEO).

| Dataset | GEO accession |
| --- | --- |
| RRBS | **GSE341998** |
| RNA-seq | **GSE341997** |

---

# 📁 Required input files

## Cocaine-associated DMRs

`hotspots.csv` contains the genomic coordinates and annotations of the 24
cocaine-associated sperm DMRs analyzed in the manuscript and corresponds to
the regions reported in Table S1.

This file is included in the repository and is used as input by analyses
requiring the DMR coordinates, including the histone ChIP-seq randomization
analyses.

---

## 🧬 RRBS CX reports

The original Bismark `CX_report` files are used in analyses requiring the
RRBS-accessible genomic space, including the transcription factor motif
analysis and the RRBS-restricted histone ChIP-seq randomization analysis.

Download the corresponding:

`*_CX_report.txt.gz`

files from **GEO: GSE341998** and place them in the working directory before
running these analyses.

---

## 🧪 Histone ChIP-seq data

The histone analyses use publicly available sperm histone ChIP-seq signal
tracks from **GEO: GSE79227**.

The following seven histone modifications are analyzed:

- H3K9ac
- H3K27ac
- H3K4me1
- H3K4me3
- H3K36me3
- H3K27me3
- H3K9me3

Two replicate WIG tracks are analyzed for each histone modification.

Place the downloaded WIG files inside:

`histone data/`

before running the histone-analysis scripts.

---

---

# 💻 Software

Analyses were performed using:

- R 4.5.1
- Bioconductor
- DESeq2
- GenomicRanges
- GenomeInfoDb
- rtracklayer
- motifmatchr
- JASPAR2024
- tidyverse
- pheatmap
- EnhancedVolcano

---

# 🔬 Genome assemblies

The cocaine-associated DMRs and primary RRBS analyses use the
**GRCm39/mm39** mouse genome assembly.

Public sperm histone ChIP-seq tracks from GSE79227 are provided in **mm9**.

For the histone analyses, DMRs and random background regions are converted
from mm39 to mm9 using liftOver. Only uniquely mapped regions on standard
chromosomes are retained. Region-length preservation is required during the
genome-wide randomization procedure.

Genome assembly conversions are therefore performed explicitly within the
histone-analysis pipeline.

---uto Tecnológico de Buenos Aires (ITBA)  
Argentina
