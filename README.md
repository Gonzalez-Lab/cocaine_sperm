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
- 🎲 Genome-wide and RRBS-restricted histone randomization analyses
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

# 🧪 Histone ChIP-seq randomization analyses

Histone-modification ChIP-seq signal across the 24 cocaine-associated DMRs is
evaluated using two complementary empirical background models.

## Genome-wide background

The primary randomization analysis compares the continuous ChIP-seq signal
observed across the 24 DMRs with 10,000 independent random genomic sets.

Each random set contains 24 genomic regions matched exactly to the lengths of
the observed DMRs. Cocaine-associated DMRs themselves are excluded from the
background.

For each histone modification, the mean `log2(signal + 1)` across the 24 DMRs
is compared with the empirical distribution of the same statistic obtained
from the 10,000 random genomic sets.

An omnibus test is additionally performed across the seven histone
modifications. The statistic for each modification is standardized relative
to its own empirical random distribution, and the seven resulting Z scores
are averaged to obtain a global statistic. The same random genomic set is
retained across all seven histone modifications when constructing each global
null statistic, preserving the dependence among histone signals.

## RRBS-restricted background

A complementary sensitivity analysis repeats the same ChIP-seq signal
quantification and statistical framework using random regions sampled from
the genomic space represented by RRBS.

This analysis provides a more conditional null model by asking whether the
histone-modification signal observed across the cocaine-associated DMRs is
elevated relative to regions accessible to the methylation assay.

The genome-wide and RRBS-restricted analyses therefore address related but
distinct null hypotheses.

## Descriptive histone-signal heatmap

For visualization, high relative histone ChIP-seq signal is defined using the
75th percentile of the independent genome-wide background distribution for
each histone modification.

This binary classification is used only for the descriptive heatmap and is
not used for statistical inference.

All statistical inference is based on the continuous ChIP-seq signal and the
empirical randomization analyses described above.

The descriptive heatmap shown in Figure 1D is generated from the genome-wide analysis 
using the 75th percentile of the independent genome-wide background. The RRBS-restricted analysis 
is used exclusively as a sensitivity analysis of the continuous ChIP-seq signal.

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

---

# ♻️ Reproducibility

Random seeds are fixed in the randomization analyses to ensure reproducible
generation of empirical background distributions.

Intermediate histone-signal calculations can be stored as RDS cache files so
that the WIG tracks do not need to be re-imported when rerunning downstream
statistical analyses.

The repository contains separate scripts for the genome-wide and
RRBS-restricted histone randomization analyses.

---


# 📧 Contact

**Betina Gonzalez**

Laboratory of Epigenetics and Functional Genomics  
Instituto Tecnológico de Buenos Aires (ITBA)  
Argentina
