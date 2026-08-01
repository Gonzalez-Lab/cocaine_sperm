################################################################################
# HISTONE ChIP SIGNAL IN COCAINE-ASSOCIATED METHYLATION HOTSPOTS
#
# Hotspots: mm39
# Histone WIG files: mm9
################################################################################


# ==============================================================================
# 1. LIBRARIES
# ==============================================================================

library(rtracklayer)
library(GenomicRanges)
library(IRanges)
library(pheatmap)
library(R.utils)


# ==============================================================================
# 2. PARAMETERS
# ==============================================================================

hotspot_file <- "hotspots.csv"
data_directory <- "histone data"

chain_gz <- "mm39ToMm9.over.chain.gz"
chain_file <- "mm39ToMm9.over.chain"

standard_chr <- paste0("chr", c(1:19, "X", "Y"))


# ==============================================================================
# 3. READ HOTSPOTS
# ==============================================================================

hotspots <- read.csv2(
  hotspot_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "chr",
  "start",
  "end",
  "hotspot",
  "gene_name",
  "dir.Methy"
)

missing_columns <- setdiff(
  required_columns,
  colnames(hotspots)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing columns in hotspots.csv: ",
    paste(missing_columns, collapse = ", ")
  )
}


# Chromosome names
hotspots$chr <- as.character(hotspots$chr)

hotspots$chr <- ifelse(
  grepl("^chr", hotspots$chr),
  hotspots$chr,
  paste0("chr", hotspots$chr)
)


# Coordinates
hotspots$start <- as.integer(hotspots$start)
hotspots$end <- as.integer(hotspots$end)

if (anyNA(hotspots$start) || anyNA(hotspots$end)) {
  stop("Some hotspot coordinates are missing or non-numeric.")
}

if (any(hotspots$start < 1)) {
  stop("Some hotspot start coordinates are smaller than 1.")
}

if (any(hotspots$end < hotspots$start)) {
  stop("Some hotspot end coordinates are smaller than their start coordinates.")
}

# Hotspot IDs
hotspots$hotspot <- as.character(hotspots$hotspot)
hotspots$gene_name <- as.character(hotspots$gene_name)

if (anyDuplicated(hotspots$hotspot)) {
  stop("Duplicated hotspot identifiers were found.")
}

hotspots_gr <- GRanges(
  seqnames = hotspots$chr,
  ranges = IRanges(
    start = hotspots$start,
    end = hotspots$end
  ),
  hotspot_id = hotspots$hotspot,
  gene_name = hotspots$gene_name
)


hotspot_labels <- paste(
  hotspots$chr,
  hotspots$hotspot,
  hotspots$gene_name,
  sep = "_"
)

if (anyDuplicated(hotspot_labels)) {
  stop("Duplicated combined hotspot labels were found.")
}

cat("Hotspots loaded:", length(hotspots_gr), "\n")

# ==============================================================================
# 4. HISTONE FILE INFORMATION
# ==============================================================================

files_chip <- list(
  
  H3K9ac = list(
    gsm = c("GSM2088389", "GSM2401437"),
    files = c(
      "GSM2088389_SPERM_H3K9ac.wig.gz",
      "GSM2401437_SPERM_H3K9ac_replicate2.wig.gz"
    )
  ),
  
  H3K27ac = list(
    gsm = c("GSM2088387", "GSM2401435"),
    files = c(
      "GSM2088387_SPERM_H3K27AC.wig.gz",
      "GSM2401435_SPERM_H3K27AC_replicate2.wig.gz"
    )
  ),
  
  H3K4me1 = list(
    gsm = c("GSM2088390", "GSM2401438"),
    files = c(
      "GSM2088390_SPERM_H3K4me1.wig.gz",
      "GSM2401438_SPERM_H3K4me1_replicate2.wig.gz"
    )
  ),
  
  H3K4me3 = list(
    gsm = c("GSM2088391", "GSM2401439"),
    files = c(
      "GSM2088391_Sperm_H3K4me3.wig.gz",
      "GSM2401439_SPERM_H3K4me3_replicate2.wig.gz"
    )
  ),
  
  H3K36me3 = list(
    gsm = c("GSM2088385", "GSM2401433"),
    files = c(
      "GSM2088385_SPERM_H3K36me3.wig.gz",
      "GSM2401433_SPERM_H3K36me3_replicate2.wig.gz"
    )
  ),
  
  H3K27me3 = list(
    gsm = c("GSM2088386", "GSM2401434"),
    files = c(
      "GSM2088386_SPERM_H3K27me3.wig.gz",
      "GSM2401434_SPERM_H3K27me3_replicate2.wig.gz"
    )
  ),
  
  H3K9me3 = list(
    gsm = c("GSM2088388", "GSM2401436"),
    files = c(
      "GSM2088388_SPERM_H3K9me3.wig.gz",
      "GSM2401436_SPERM_H3K9me3_replicate2.wig.gz"
    )
  )
)


# ==============================================================================
# 5. LOCATE WIG FILES IN LOCAL HISTONE DATA FOLDER
# ==============================================================================

if (!dir.exists(data_directory)) {
  stop(
    "The folder '", data_directory,
    "' was not found. Current working directory: ",
    getwd()
  )
}

paths_chip <- lapply(
  files_chip,
  function(x) {
    
    paths <- file.path(
      data_directory,
      x$files
    )
    
    missing_files <- paths[!file.exists(paths)]
    
    if (length(missing_files) > 0) {
      stop(
        "Missing histone files:\n",
        paste(missing_files, collapse = "\n")
      )
    }
    
    paths
  }
)

cat("All histone WIG files were found in:", data_directory, "\n")


# ==============================================================================
# 6. DOWNLOAD AND IMPORT mm39 -> mm9 CHAIN
# ==============================================================================

if (!file.exists(chain_file)) {
  
  if (!file.exists(chain_gz)) {
    
    download.file(
      "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/liftOver/mm39ToMm9.over.chain.gz",
      chain_gz,
      mode = "wb"
    )
  }
  
  gunzip(
    chain_gz,
    destname = chain_file,
    remove = FALSE,
    overwrite = TRUE
  )
}

chain_mm39_to_mm9 <- import.chain(chain_file)


# ==============================================================================
# 7. LIFT HOTSPOTS mm39 -> mm9
# ==============================================================================

hotspots_mm9_list <- liftOver(
  hotspots_gr,
  chain_mm39_to_mm9
)

# Retain only hotspots with exactly one lifted interval
unique_lift <- elementNROWS(hotspots_mm9_list) == 1

hotspots_mm9 <- unlist(
  hotspots_mm9_list[unique_lift],
  use.names = FALSE
)

# Retain standard chromosomes
hotspots_mm9 <- hotspots_mm9[
  as.character(seqnames(hotspots_mm9)) %in% standard_chr
]

if (!"hotspot_id" %in% colnames(mcols(hotspots_mm9))) {
  stop("hotspot_id metadata were lost during liftOver.")
}

if (anyDuplicated(hotspots_mm9$hotspot_id)) {
  stop("Duplicated hotspot IDs were generated by liftOver.")
}


not_lifted_ids <- setdiff(
  hotspots_gr$hotspot_id,
  hotspots_mm9$hotspot_id
)

cat(
  "Hotspots retained after liftOver:",
  length(hotspots_mm9),
  "of",
  length(hotspots_gr),
  "\n"
)

if (length(not_lifted_ids) > 0) {
  warning(
    "Hotspots not evaluable in mm9: ",
    paste(not_lifted_ids, collapse = ", ")
  )
}


# ==============================================================================
# 8. CLEAN WIG TRACK
# ==============================================================================

clean_track <- function(track) {
  
  chromosome_names <- seqlevels(track)
  
  # Examples:
  # chr1      -> chr1
  # mm9.chr1  -> chr1
  cleaned_names <- sub(
    "^.*\\.(chr[0-9]+|chrX|chrY|chrM)$",
    "\\1",
    chromosome_names
  )
  
  if (!anyDuplicated(cleaned_names)) {
    seqlevels(track) <- cleaned_names
  }
  
  track <- track[
    as.character(seqnames(track)) %in% standard_chr
  ]
  
  if (!"score" %in% colnames(mcols(track))) {
    stop("The imported WIG file does not contain a score column.")
  }
  
  track$score <- as.numeric(track$score)
  
  track <- track[
    !is.na(track$score)
  ]
  
  track
}


# ==============================================================================
# 9. WEIGHTED MEAN SIGNAL WITHIN HOTSPOTS
# ==============================================================================

signal_over_regions <- function(regions, track) {
  
  signal <- numeric(length(regions))
  
  hits <- findOverlaps(
    regions,
    track,
    ignore.strand = TRUE
  )
  
  if (length(hits) == 0) {
    return(signal)
  }
  
  overlap_ranges <- pintersect(
    regions[queryHits(hits)],
    track[subjectHits(hits)],
    ignore.strand = TRUE
  )
  
  weighted_signal <- width(overlap_ranges) *
    track$score[subjectHits(hits)]
  
  signal_sum <- rowsum(
    weighted_signal,
    group = queryHits(hits),
    reorder = FALSE
  )
  
  indices <- as.integer(rownames(signal_sum))
  
  signal[indices] <- signal_sum[, 1] / width(regions)[indices]
  
  signal
}


# ==============================================================================
# 10. ANALYZE ONE HISTONE MARK
# ==============================================================================

analyze_histone_mark <- function(
    mark_name,
    files,
    hotspots_original,
    hotspots_lifted
) {
  
  message("Processing ", mark_name)
  
  track_1 <- clean_track(import(files[1]))
  track_2 <- clean_track(import(files[2]))
  
  lifted_signal_1 <- signal_over_regions(
    hotspots_lifted,
    track_1
  )
  
  lifted_signal_2 <- signal_over_regions(
    hotspots_lifted,
    track_2
  )
  
  
  # NA means that the hotspot could not be evaluated after liftOver
  signal_1 <- rep(NA_real_, length(hotspots_original))
  signal_2 <- rep(NA_real_, length(hotspots_original))
  
  original_index <- match(
    hotspots_lifted$hotspot_id,
    hotspots_original$hotspot_id
  )
  
  signal_1[original_index] <- lifted_signal_1
  signal_2[original_index] <- lifted_signal_2
  
  
  signal_mean <- rowMeans(
    cbind(signal_1, signal_2)
  )
  
  
  positive_signal <- signal_mean[
    !is.na(signal_mean) &
      signal_mean > 0
  ]
  
  
  threshold <- if (length(positive_signal) > 0) {
    as.numeric(
      quantile(
        positive_signal,
        probs = 0.75,
        names = FALSE
      )
    )
  } else {
    Inf
  }
  
  
  signal_high <- rep(NA_integer_, length(signal_mean))
  
  evaluable <- !is.na(signal_mean)
  
  signal_high[evaluable] <- as.integer(
    signal_mean[evaluable] > 0 &
      signal_mean[evaluable] >= threshold
  )
  
  
  data.frame(
    hotspot_id = hotspots_original$hotspot_id,
    rep1 = signal_1,
    rep2 = signal_2,
    signal_mean = signal_mean,
    signal_high = signal_high,
    threshold_q75 = threshold
  )
}


# ==============================================================================
# 11. ANALYZE ALL HISTONE MARKS
# ==============================================================================

results <- lapply(
  names(paths_chip),
  function(mark_name) {
    analyze_histone_mark(
      mark_name = mark_name,
      files = paths_chip[[mark_name]],
      hotspots_original = hotspots_gr,
      hotspots_lifted = hotspots_mm9
    )
  }
)

names(results) <- names(paths_chip)


# ==============================================================================
# 12. CONTINUOUS SIGNAL MATRIX
# ==============================================================================

histone_signal <- do.call(
  cbind,
  lapply(results, `[[`, "signal_mean")
)

colnames(histone_signal) <- names(results)
rownames(histone_signal) <- hotspot_labels


# ==============================================================================
# 13. CONTINUOUS HEATMAP
# ==============================================================================

histone_log <- log2(histone_signal + 1)

pheatmap(
  t(histone_log),
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Relative histone ChIP signal across cocaine hotspots",
  border_color = "grey80",
  fontsize_row = 10,
  fontsize_col = 8,
  na_col = "grey80"
)


# ==============================================================================
# 14. BINARY HIGH-SIGNAL MATRIX
# ==============================================================================

histones_high <- do.call(
  cbind,
  lapply(results, `[[`, "signal_high")
)

colnames(histones_high) <- names(results)
rownames(histones_high) <- hotspot_labels


#export file for further analysis
write.csv2(
  histones_high,
  "histone high signal binary.csv",
  na = "NA"
)


# ==============================================================================
# 15. THRESHOLDS
# ==============================================================================

histone_thresholds <- data.frame(
  Histone_mark = names(results),
  threshold_q75 = vapply(
    results,
    function(x) unique(x$threshold_q75),
    numeric(1)
  )
)


# ==============================================================================
# 16. HOTSPOT ANNOTATION
# ==============================================================================

hotspots$dir.Methy <- tolower(
  trimws(hotspots$dir.Methy)
)

annotation_col <- data.frame(
  Direction = hotspots$dir.Methy,
  row.names = hotspot_labels
)

direction_values <- unique(hotspots$dir.Methy)

annotation_colors <- list(
  Direction = setNames(
    rep(c("deeppink3", "dodgerblue4"), length.out = length(direction_values)),
    direction_values
  )
)


# ==============================================================================
# 17. BINARY HEATMAP
# ==============================================================================

pheatmap(
  t(histones_high),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = c("floralwhite", "cadetblue1"),
  breaks = c(-0.5, 0.5, 1.5),
  na_col = "grey80",
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  main = "High relative histone ChIP signal across cocaine hotspots",
  border_color = "grey80",
  fontsize_row = 10,
  fontsize_col = 8,
  legend = FALSE
)


# ==============================================================================
# 18. SUMMARY
# ==============================================================================

cat("\nHistone analysis completed\n")

cat(
  "Original hotspots:",
  length(hotspots_gr),
  "\n"
)

cat(
  "Hotspots evaluable after liftOver:",
  length(hotspots_mm9),
  "\n"
)

cat(
  "Hotspots not evaluable:",
  length(not_lifted_ids),
  "\n"
)

cat("\nHigh-signal hotspots per histone mark:\n")

print(
  colSums(
    histones_high,
    na.rm = TRUE
  )
)

cat("\nNon-evaluable hotspots per histone mark:\n")

print(
  colSums(
    is.na(histones_high)
  )
)