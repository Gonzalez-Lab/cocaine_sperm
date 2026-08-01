###############################################################################
# Distance between cocaine-associated DMR hotspots and paternally imprinted genes
# Genome assembly: GRCm39
###############################################################################

library(GenomicRanges)
library(IRanges)
library(biomaRt)
library(GenomeInfoDb)
library(dplyr)

#----------------------------------------------------------
# 1. Load cocaine-associated DMR hotspots
#----------------------------------------------------------

hotspots <- read.csv2(
  "hotspots.csv",
  stringsAsFactors = FALSE
)

# Check chromosome format before adding "chr"
head(hotspots$chr)

# Add "chr" only if chromosome names do not already contain it
hotspots$chr <- ifelse(
  grepl("^chr", hotspots$chr),
  hotspots$chr,
  paste0("chr", hotspots$chr)
)

# Convert hotspots to GRanges
hotspots_gr <- GRanges(
  seqnames = hotspots$chr,
  ranges = IRanges(
    start = hotspots$start,
    end = hotspots$end
  ),
  diffMethy = hotspots$diff.Methy,
  dirMethy = hotspots$dir.Methy,
  gene = hotspots$gene,
  hotspot = hotspots$hotspot
)

# Explicitly indicate genome assembly
genome(hotspots_gr) <- "GRCm39"


#----------------------------------------------------------
# 2. Load curated paternally imprinted genes
#----------------------------------------------------------

imprinting <- read.csv2(
  "imprinting.csv",
  fileEncoding = "latin1",
  stringsAsFactors = FALSE
)

# Remove empty or duplicated gene symbols
imprinting <- imprinting %>%
  filter(
    !is.na(Gene),
    Gene != ""
  ) %>%
  distinct(Gene, .keep_all = TRUE)


#----------------------------------------------------------
# 3. Retrieve GRCm39 gene coordinates from Ensembl
#----------------------------------------------------------

options(timeout = 300)

mart <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl"
)

coords <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand",
    "gene_biotype"
  ),
  filters = "external_gene_name",
  values = imprinting$Gene,
  mart = mart
)

# Retain standard chromosomes and one genomic interval per Ensembl gene
coords <- coords %>%
  filter(
    chromosome_name %in% c(as.character(1:19), "X", "Y"),
    external_gene_name != ""
  ) %>%
  distinct(
    ensembl_gene_id,
    .keep_all = TRUE
  )

# Report genes that were not found in Ensembl
genes_not_found <- setdiff(
  imprinting$Gene,
  coords$external_gene_name
)

if (length(genes_not_found) > 0) {
  message(
    "Genes not found in Ensembl: ",
    paste(genes_not_found, collapse = ", ")
  )
}

# Inspect retrieved coordinates
head(coords)


#----------------------------------------------------------
# 4. Convert imprinted genes to GRanges
#----------------------------------------------------------

imprinting_gr <- GRanges(
  seqnames = paste0("chr", coords$chromosome_name),
  ranges = IRanges(
    start = coords$start_position,
    end = coords$end_position
  ),
  strand = ifelse(
    coords$strand == 1,
    "+",
    "-"
  ),
  ensembl_gene_id = coords$ensembl_gene_id,
  gene_symbol = coords$external_gene_name,
  gene_biotype = coords$gene_biotype
)

genome(imprinting_gr) <- "GRCm39"


#----------------------------------------------------------
# 5. Recreate hotspots_gr to preserve all 24 hotspots
#----------------------------------------------------------

hotspots_gr <- GRanges(
  seqnames = hotspots$chr,
  ranges = IRanges(
    start = hotspots$start,
    end = hotspots$end
  ),
  diffMethy = hotspots$diff.Methy,
  dirMethy = hotspots$dir.Methy,
  gene = hotspots$gene,
  hotspot = hotspots$hotspot
)

genome(hotspots_gr) <- "GRCm39"

# Check chromosome naming
seqlevels(hotspots_gr)
seqlevels(imprinting_gr)

# Both objects should use the same format:
# chr1, chr2, ..., chr19, chrX, chrY


#----------------------------------------------------------
# 6. Find nearest paternally imprinted gene
#----------------------------------------------------------

nearest_hits <- distanceToNearest(
  hotspots_gr,
  imprinting_gr,
  ignore.strand = TRUE
)

query_idx <- queryHits(nearest_hits)
subject_idx <- subjectHits(nearest_hits)

# Initialize one result for every hotspot
nearest_gene <- rep(NA_character_, length(hotspots_gr))
nearest_ensembl <- rep(NA_character_, length(hotspots_gr))
distance_bp <- rep(NA_integer_, length(hotspots_gr))

# Add results only for hotspots with a valid nearest gene
nearest_gene[query_idx] <-
  imprinting_gr$gene_symbol[subject_idx]

nearest_ensembl[query_idx] <-
  imprinting_gr$ensembl_gene_id[subject_idx]

distance_bp[query_idx] <-
  mcols(nearest_hits)$distance


#----------------------------------------------------------
# 7. Return results to the original hotspots table
#----------------------------------------------------------

hotspots$nearest_imprinted_gene <- nearest_gene

hotspots$nearest_imprinted_ensembl <- nearest_ensembl

hotspots$dist_imprinting_bp <- distance_bp

hotspots$dist_imprinting_kb <-
  hotspots$dist_imprinting_bp / 1000


#----------------------------------------------------------
# 8. Inspect results
#----------------------------------------------------------

hotspots %>%
  dplyr::select(
    hotspot,
    chr,
    start,
    end,
    nearest_imprinted_gene,
    dist_imprinting_bp,
    dist_imprinting_kb
  ) %>%
  arrange(dist_imprinting_bp)