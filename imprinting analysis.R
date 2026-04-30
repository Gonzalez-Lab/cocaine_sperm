###############################################################################
#paternal imprinting
###############################################################################

library(GenomicRanges)
library(data.table)
library(biomaRt)
library(dplyr)

#upload cocaine DMRs hotspots
hotspots <- read.csv2("hotspots.csv")

#convert to granges
hotspots_gr <- GRanges(
  seqnames = paste0("chr", hotspots$chr),
  ranges = IRanges(hotspots$start, hotspots$end),
  diffMethy = hotspots$diff.Methy,
  dirMethy = hotspots$Dir.Methy,
  gene = hotspots$gene,
  hotspot = hotspots$hotspot
)

#upload mouse paternal imprinting genes curated from Geneimprint database
imprinting <- read.csv2("imprinting.csv", fileEncoding = "latin1")

# Conectar a Ensembl (mouse)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# obtain mm10 coordinates
options(timeout = 300)  # 5 minutos
coords <- getBM(
  attributes = c("external_gene_name",
                 "chromosome_name",
                 "start_position",
                 "end_position",
                 "strand",
                 "gene_biotype"),
  filters = "external_gene_name",
  values = imprinting$Gene,
  mart = mart
)

# Limpiar cromosomas raros
coords <- coords %>%
  filter(chromosome_name %in% c(1:19, "X", "Y")) %>%
  distinct()

head(coords)


imprinting_gr <- GRanges(
  seqnames = paste0("chr",coords$chromosome_name),
  ranges = IRanges(coords$start_position, coords$end_position),
  locus = coords$gene_biotype
)

dist_imprinting <- distanceToNearest(hotspots_gr, imprinting_gr)

hotspots_gr$dist_imprinting <- NA

hits <- queryHits(dist_imprinting)

hotspots_gr$dist_imprinting[hits] <- mcols(dist_imprinting)$distance


hotspots$dist_imprinting_kb <- hotspots_gr$dist_imprinting/1000
