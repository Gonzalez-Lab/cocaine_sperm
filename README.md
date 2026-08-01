Cocaine-induced epigenetic remodeling of mouse sperm

This repository contains the analysis scripts used in the manuscript:

"Cocaine reshapes the sperm epigenome through coordinated DNA methylation and RNA remodeling in mice"

The repository includes all custom R scripts used for the analyses presented in the manuscript, including:

RRBS hotspot annotation
Histone ChIP-seq integration
Transcription factor motif analysis
RNA-seq differential expression analysis
Figure generation
Supplementary analyses
Data availability

The complete sequencing datasets have been deposited in the NCBI Gene Expression Omnibus (GEO):

RRBS: GSE341998
RNA-seq: GSE341997
Required external files

Some analyses require files that are not included in this repository because of their large size.

RRBS CX reports

The scripts:

TF_analysis.R
TE_enrichment.R

require the original Bismark *_CX_report.txt.gz files generated during the RRBS analysis.

These files are available as part of the RRBS GEO submission (GSE341998).

Download the corresponding CX_report files from GEO and place them in the working directory before running these scripts.

Histone ChIP-seq data

The script:

Histone_analysis.R

requires the processed BigWig files containing sperm histone ChIP-seq signal.

These files should be downloaded from the original GEO study:

GSE79227

The required BigWig files should be placed inside the folder:

histone data/

before running the analysis.

Software

Analyses were performed using:

R 4.5.1
Bioconductor
GenomicRanges
rtracklayer
motifmatchr
JASPAR2024
DESeq2
pheatmap
EnhancedVolcano
tidyverse

Additional package dependencies are listed within each script.

Reproducibility

All scripts were developed using the GRCm39/mm39 mouse genome assembly.

The repository corresponds to the revised version of the manuscript following peer review, including:

complete RRBS reanalysis using a fully GRCm39-consistent workflow,
transcription factor motif analysis,
updated histone integration analysis,
revised RNA-seq analyses,
generation of all manuscript figures.

<img width="672" height="480" alt="fractal-large" src="https://github.com/user-attachments/assets/fc21c83d-6d40-4465-851a-c80312ebd654" />
