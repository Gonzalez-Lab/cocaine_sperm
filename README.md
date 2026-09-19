# 🧬 Cocaine-induced epigenetic remodeling of mouse sperm

Repository containing the bioinformatic analyses used in the manuscript:

> **Cocaine reshapes the sperm epigenome through coordinated DNA methylation and RNA remodeling in mice**

---

# 📂 Repository contents

This repository contains the custom scripts and input files used to generate the analyses presented in the manuscript.

## Included analyses

- 🧬 RRBS hotspot analysis
- 🎯 Transcription factor motif analysis
- 🧪 Histone ChIP-seq integration and genomic background randomization
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
requiring the DMR coordinates, including the histone ChIP-seq analysis.

---

## 🧬 RRBS CX reports

The transcription factor motif analysis requires the original Bismark
`CX_report` files.

Download the corresponding:

*_CX_report.txt.gz

files from **GEO: GSE341998** and place them in the working directory before
running the analysis.

---

## 🧪 Histone ChIP-seq data

The histone analysis uses publicly available sperm histone ChIP-seq signal
tracks from **GEO: GSE79227**.

The analysis includes the following seven histone modifications:

- H3K9ac
- H3K27ac
- H3K4me1
- H3K4me3
- H3K36me3
- H3K27me3
- H3K9me3

Two replicate WIG tracks are analyzed for each histone modification.

Place the downloaded WIG files inside:

histone data/

before running:

Histone_analysis.R

The histone analysis compares continuous ChIP-seq signal across the 24
cocaine-associated DMRs with 10,000 sets of length-matched random genomic
regions. The 75th percentile of the independent genomic background is used
only to generate the descriptive binary heatmap of high relative histone
ChIP-seq signal.

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
For the histone analysis, DMRs and random genomic regions are converted from
mm39 to mm9 using liftOver. Only uniquely mapped regions on standard
chromosomes that preserve the original region length are retained.

Genome assembly conversions are therefore performed explicitly within the
histone-analysis pipeline.

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
