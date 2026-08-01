# ==============================================================================
# TF motif enrichment RRBS-AWARE
# Mus musculus mm39 / GRCm39
#
# This code construct a DMR universe based on the RRBS cx files from the experiment
# selecting 200 random DMRs for each cocaine hotspot, matching for chromosome,
# length and CpG density

# ==============================================================================
# 1. libraries
# ==============================================================================

library(data.table)
library(GenomicRanges)
library(IRanges)
library(Biostrings)
library(motifmatchr)
library(JASPAR2024)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm39)
library(circlize)
library(RSQLite)
library(grid)
library(pheatmap)

# ==============================================================================
# 2. parameters
# ==============================================================================

hotspot_file <- "hotspots.csv" #file with cocaine DMRs or methylation hotspots

cx_files <- c(                      # files with all CpGs sampled in RRBS experiment
  "control_1_CX_report.txt.gz",
  "control_2_CX_report.txt.gz",
  "cocaine_1_CX_report.txt.gz",
  "cocaine_2_CX_report.txt.gz"
)

# regions are scanned in a 200 bp window
window_width <- 200L

# Matching RRBS
length_tolerance <- 0.15
density_tolerance <- 0.25
cpg_tolerance <- 2L

# best candidates for each DMR
top_candidates <- 200L

# Treshold for scanned motifs
motif_p_cutoff <- 1e-5

# Permutations
n_perm <- 10000L
random_seed <- 111L

# Chromosomes shared in Ensembl GRCm39 y UCSC mm39
standard_chr <- paste0("chr", c(1:19, "X", "Y"))

# define genome mm39
genome_mm39 <- BSgenome.Mmusculus.UCSC.mm39

# ==============================================================================
# 3. upload cx reports
# ==============================================================================

cx_list <- lapply(cx_files, function(file_name) {
    fread(file_name)
  }
)

cx <- rbindlist(cx_list, use.names = TRUE, fill = TRUE)

rm(cx_list)

gc()

setnames(
  cx,
  old = names(cx)[1:6],
  new = c(
    "chr",
    "pos",
    "strand",
    "meth_count",
    "unmeth_count",
    "context"
  )
)

cx[, chr := as.character(chr)]
cx[, pos := as.integer(pos)]
cx[, meth_count := as.numeric(meth_count)]
cx[, unmeth_count := as.numeric(unmeth_count)]
cx[, context := as.character(context)]

cx[
  !grepl("^chr", chr),
  chr := paste0("chr", chr)
]

cx[
  chr == "chrMT",
  chr := "chrM"
]

cx[
  ,
  coverage := meth_count + unmeth_count
]


# ==============================================================================
# 4. obtain covered CpGs
# ==============================================================================

cx_cpg <- cx[
  context == "CG" &
    coverage > 0 &
    chr %in% standard_chr &
    strand %in% c("+", "-")
]

# Collapse the two strand-specific cytosines of each CpG dyad to the
# coordinate of the cytosine on the forward strand.
cx_cpg[
  ,
  cpg_pos := fifelse(
    strand == "-",
    pos - 1L,
    pos
  )
]

cpgs <- unique(
  cx_cpg[
    cpg_pos >= 1L,
    .(
      chr,
      pos = cpg_pos
    )
  ]
)

setorder(
  cpgs,
  chr,
  pos
)

cat(
  "\nUnique covered CpGs:",
  nrow(cpgs),
  "\n"
)

rm(cx, cx_cpg)
gc()


# ==============================================================================
# 5. RRBS universe construction
# ==============================================================================

# We use same DSS parameters:
#
# minCG = 3
# minlen = 50
# dis.merge = 100

cpgs[
  ,
  gap := c(
    NA_integer_,
    diff(pos)
  ),
  by = chr
]

cpgs[
  ,
  new_region := is.na(gap) | gap > 100L
]

cpgs[
  ,
  region_number := cumsum(new_region),
  by = chr
]

rrbs_universe <- cpgs[
  ,
  .(
    start = min(pos),
    end = max(pos),
    n_cpg = .N
  ),
  by = .(
    chr,
    region_number
  )
]

rrbs_universe[
  ,
  length_bp := end - start + 1L
]

rrbs_universe[
  ,
  cpg_density := n_cpg / length_bp
]

rrbs_universe <- rrbs_universe[
  n_cpg >= 3L &
    length_bp >= 50L
]

rrbs_universe[
  ,
  region_id := paste0(
    chr,
    "_RRBS_",
    region_number
  )
]

setcolorder(
  rrbs_universe,
  c(
    "region_id",
    "chr",
    "start",
    "end",
    "length_bp",
    "n_cpg",
    "cpg_density",
    "region_number"
  )
)

setorder(
  rrbs_universe,
  chr,
  start
)

cat(
  "Regions in RRBS universe:",
  nrow(rrbs_universe),
  "\n"
)

# ==============================================================================
# 6. upload hotspots
# ==============================================================================

hotspots <- fread(hotspot_file, sep = "auto", header = TRUE)

hotspots[, chr := as.character(chr)]

hotspots[
  !grepl("^chr", chr),
  chr := paste0("chr", chr)
]

if (!"length_bp" %in% names(hotspots)) {
  hotspots[, length_bp := end - start + 1L]
}

# we use coordinate-derived length.
hotspots[, length_bp := end - start + 1L]
hotspots[, cpg_density := n_cpg / length_bp]

# Use the hotspot nomenclature already defined in hotspots.csv
hotspots[, hotspot := as.character(hotspot)]

#convert hotspots to granges object
hotspots_gr <- GRanges(
  seqnames = hotspots$chr,
  ranges = IRanges(
    start = hotspots$start,
    end = hotspots$end
  )
)

mcols(hotspots_gr)$hotspot <- hotspots$hotspot

# ==============================================================================
# 7. exclude hotspots DMRs from backgroud
# ==============================================================================

rrbs_gr_all <- GRanges(
  seqnames = rrbs_universe$chr,
  ranges = IRanges(
    start = rrbs_universe$start,
    end = rrbs_universe$end
  )
)

overlap_real <- overlapsAny(
  rrbs_gr_all,
  hotspots_gr,
  ignore.strand = TRUE
)

rrbs_background <- rrbs_universe[!overlap_real]

cat(
  "Available RRBS background regions:",
  nrow(rrbs_background),
  "\n"
)

# ==============================================================================
# 8. construct matched RRBS pools
# ==============================================================================

candidate_pools <- vector(mode = "list",length = nrow(hotspots))

for (i in seq_len(nrow(hotspots))) {
  
  h_chr <- hotspots$chr[i]
  h_length <- hotspots$length_bp[i]
  h_n_cpg <- hotspots$n_cpg[i]
  h_density <- hotspots$cpg_density[i]
  
  candidates <- rrbs_background[
    chr == h_chr &
      length_bp >= h_length * (1 - length_tolerance) &
      length_bp <= h_length * (1 + length_tolerance) &
      abs(n_cpg - h_n_cpg) <= cpg_tolerance &
      cpg_density >= h_density * (1 - density_tolerance) &
      cpg_density <= h_density * (1 + density_tolerance)
  ]
  
  if (nrow(candidates) == 0L) {
    stop(
      "No RRBS candidates were found for ",
      hotspots$hotspot[i]
    )
  }
  
  candidates[
    ,
    match_score :=
      abs(length_bp - h_length) / h_length +
      abs(n_cpg - h_n_cpg) / max(h_n_cpg, 1L) +
      abs(cpg_density - h_density) / max(h_density, 1e-10)
  ]
  
  setorder(
    candidates,
    match_score
  )
  
  n_keep <- min(
    top_candidates,
    nrow(candidates)
  )
  
  candidate_pools[[i]] <- copy(
    candidates[
      seq_len(n_keep)
    ]
  )
  
  cat(
    hotspots$hotspot[i],
    ":",
    nrow(candidates),
    "compatible candidates conserved",
    n_keep,
    "\n"
  )
}


# ==============================================================================
# 9. join unique candidates
# ==============================================================================

candidate_universe <- unique(
  rbindlist(
    candidate_pools,
    use.names = TRUE
  ),
  by = "region_id"
)

candidate_gr <- GRanges(
  seqnames = candidate_universe$chr,
  ranges = IRanges(
    start = candidate_universe$start,
    end = candidate_universe$end
  )
)

mcols(candidate_gr)$region_id <- candidate_universe$region_id

cat(
  "\nUnique background regions:",
  length(candidate_gr),
  "\n"
)


# ==============================================================================
# 10. function that creates 200 bp window
# ==============================================================================

make_centered_windows <- function(
    gr,
    width_bp,
    chromosome_lengths
) {
  
  chromosome <- as.character(seqnames(gr)  )
  
  midpoint <- floor(
    ( start(gr) + end(gr) ) / 2
  )
  
  new_start <- midpoint - floor((width_bp - 1L) / 2L)
  
  new_end <- new_start +  width_bp -  1L
  
  chr_lengths <- unname(chromosome_lengths[chromosome])
  
  invalid_windows <-
    is.na(chr_lengths) |
    new_start < 1L |
    new_end > chr_lengths
  
  if (any(invalid_windows)) {
    stop("some 200 bp windows exceed mm39 limits.")
  }
  
  new_gr <- GRanges(
    seqnames = chromosome,
    ranges = IRanges(
      start = new_start,
      end = new_end
    )
  )
  
  mcols(new_gr) <- mcols(gr)
  
  new_gr
}

#define mm39 chromosomes length
mm39_lengths <- seqlengths(BSgenome.Mmusculus.UCSC.mm39)

hotspot_windows <- make_centered_windows(
  gr = hotspots_gr,
  width_bp = window_width,
  chromosome_lengths = mm39_lengths
)

candidate_windows <- make_centered_windows(
  gr = candidate_gr,
  width_bp = window_width,
  chromosome_lengths = mm39_lengths
)

# Exclude background windows that overlap any scanned hotspot window.
background_overlaps_hotspots <- overlapsAny(
  candidate_windows,
  hotspot_windows,
  ignore.strand = TRUE
)

candidate_windows <- candidate_windows[
  !background_overlaps_hotspots
]

valid_background_ids <- as.character(
  mcols(candidate_windows)$region_id
)

# Prune every hotspot-specific candidate pool to the windows that remain valid.
candidate_pool_ids <- lapply(
  candidate_pools,
  function(pool) {
    intersect(
      pool$region_id,
      valid_background_ids
    )
  }
)

if (any(lengths(candidate_pool_ids) == 0L)) {
  stop(
    "At least one hotspot has no valid candidates after excluding ",
    "background windows that overlap hotspot windows."
  )
}

hotspot_window_ids <- as.character(
  mcols(hotspot_windows)$hotspot
)

candidate_window_ids <- as.character(
  mcols(candidate_windows)$region_id
)

cat("DMR windows:", length(hotspot_windows),"\n")

cat("Background windows:",length(candidate_windows),"\n")

# ==============================================================================
# 11. upload mouse JASPAR2024 motifs
# ==============================================================================

#connect to jaspar database
jaspar_connection <- dbConnect(
  SQLite(),
  db(JASPAR2024())
)

pfm_list <- getMatrixSet(
  x = jaspar_connection,
  opts = list(
    species = 10090,
    collection = "CORE",
    all_versions = FALSE
  )
)

dbDisconnect(jaspar_connection) #close conection to jaspar database

cat("\nMouse motifs JASPAR2024:",length(pfm_list),"\n")

motif_id <- ID(pfm_list)
motif_label <- name(pfm_list)

if (anyDuplicated(motif_id)) {
  stop("JASPAR returned duplicated motif IDs.")
}

# Use the JASPAR matrix ID as the unique identifier throughout the analysis.
names(pfm_list) <- motif_id


# ==============================================================================
# 12. scan motifs
# ==============================================================================

#scan motifs in hotspots
hotspot_motifs <- matchMotifs(
  pwms = pfm_list,
  subject = hotspot_windows,
  genome = genome_mm39,
  bg = "genome",
  p.cutoff = motif_p_cutoff
)

#scan motifs in background
background_motifs <- matchMotifs(
  pwms = pfm_list,
  subject = candidate_windows,
  genome = genome_mm39,
  bg = "genome",
  p.cutoff = motif_p_cutoff
)

hotspot_match_matrix <- as.matrix(motifMatches(hotspot_motifs))

background_match_matrix <- as.matrix(motifMatches(background_motifs))

colnames(hotspot_match_matrix) <- motif_id
colnames(background_match_matrix) <- motif_id

rownames(hotspot_match_matrix) <- hotspot_window_ids
rownames(background_match_matrix) <- candidate_window_ids


# ==============================================================================
# 13. prepare permutation pools
# ==============================================================================

candidate_pool_ids <- lapply(
  candidate_pool_ids,
  function(ids) {
    intersect(
      ids,
      rownames(background_match_matrix)
    )
  }
)

if (any(lengths(candidate_pool_ids) == 0L)) {
  stop("At least one hotspot has no candidates in background_match_matrix.")
}

# Sample without replacement, randomizing assignment priority in each permutation.
sample_background_set <- function(
    pool_ids,
    max_attempts = 100L
) {

  for (attempt in seq_len(max_attempts)) {

    selected_ids <- rep(NA_character_, length(pool_ids))
    used_ids <- character()
    random_order <- sample(seq_along(pool_ids))
    success <- TRUE

    for (i in random_order) {

      available <- setdiff(pool_ids[[i]], used_ids)

      if (length(available) == 0L) {
        success <- FALSE
        break
      }

      chosen <- sample(available, size = 1L)
      selected_ids[i] <- chosen
      used_ids <- c(used_ids, chosen)
    }

    if (success && !anyNA(selected_ids)) {
      return(selected_ids)
    }
  }

  stop(
    "Could not assign a unique background region to every hotspot after ",
    max_attempts,
    " attempts."
  )
}

# ==============================================================================
# 14. execute 10.000 permutations
# ==============================================================================

observed_counts <- colSums(hotspot_match_matrix)

set.seed(random_seed)

null_counts <- matrix(
  0L,
  nrow = n_perm,
  ncol = ncol(background_match_matrix)
)

colnames(null_counts) <- colnames(background_match_matrix)

for (b in seq_len(n_perm)) {
  
  selected_ids <- sample_background_set(
    candidate_pool_ids
  )
  
  null_counts[b, ] <- colSums(
    background_match_matrix[
      selected_ids,
      ,
      drop = FALSE
    ]
  )
  
  if (b %% 500L == 0L) {
    
    message(
      "permutations completed: ", b," / ",n_perm)
  }
}

# ==============================================================================
# 15. calculate enrichment
# ==============================================================================

expected_counts <- colMeans(null_counts)

empirical_p <- vapply(seq_along(observed_counts), function(j) { (
      1 +
        sum(
          null_counts[, j] >=
            observed_counts[j]
        )
    ) / (
      n_perm + 1
    )
  },
  numeric(1)
)

fold_enrichment <- (observed_counts + 0.5) / (expected_counts + 0.5)

results <- data.table(
                      motif_id = colnames(hotspot_match_matrix),
                      tf_name = motif_label[
                        match(
                          colnames(hotspot_match_matrix),
                          motif_id
                        )
                      ],
                      observed_DMRs = as.integer(observed_counts),
                      observed_fraction = as.numeric(observed_counts/nrow(hotspots)),
                      expected_DMRs = as.numeric(expected_counts),
                      expected_fraction = as.numeric(expected_counts/nrow(hotspots)),
                      fold_enrichment = as.numeric(fold_enrichment),
                      empirical_p = as.numeric(empirical_p))

results[, adjusted_p_BH := p.adjust(empirical_p,method = "BH")
]

setorder(
  results,
  adjusted_p_BH,
  empirical_p,
  -fold_enrichment,
  -observed_DMRs
)

# Exploratory candidates. No FDR-significant claim is made.
TFresults_report <- results[
    observed_DMRs >= 2L &
    fold_enrichment > 1 &
    empirical_p < 0.05
]

##############################################################
# Save complete and exploratory results.
write.csv2(results,"TF_motif_enrichment_all_results.csv")


cat(
  "\nMotifs with adjusted_p_BH < 0.05:",
  results[adjusted_p_BH < 0.05, .N],
  "\n"
)

cat(
  "Exploratory candidates:",
  nrow(TFresults_report),
  "\n"
)

###############################################################################
# heatmap

# ==============================================================================
# HEATMAP TF × DMR
# ==============================================================================

if (nrow(TFresults_report) == 0L) {
  stop("No exploratory TF candidates passed the heatmap criteria.")
}

# Locate selected JASPAR motifs in the hotspot motif matrix
tf_columns <- match(
  TFresults_report$motif_id,
  colnames(hotspot_match_matrix)
)

if (anyNA(tf_columns)) {
  stop("Some selected motifs are absent from hotspot_match_matrix.")
}

# Rows = motifs/TFs; columns = hotspots
heatmap_matrix_tf <- t(
  hotspot_match_matrix[
    ,
    tf_columns,
    drop = FALSE
  ] * 1L
)

# Use TF names for display; make.unique distinguishes duplicated TF names
rownames(heatmap_matrix_tf) <- make.unique(
  TFresults_report$tf_name
)

# Order hotspots by methylation direction and genomic position
column_order <- order(
  hotspots$dir.Methy,
  hotspots$region,
  hotspots$chr,
  hotspots$start
)

heatmap_matrix_tf <- heatmap_matrix_tf[
  ,
  column_order,
  drop = FALSE
]

# Labels: chr_HS_gene
colnames(heatmap_matrix_tf) <- paste0(
  hotspots$chr[column_order], "_",
  hotspots$hotspot[column_order], "_",
  hotspots$gene_name[column_order]
)

write.csv2(
  heatmap_matrix_tf,
  "TF binding at hotspots matrix.csv"
)

# Column annotation
annotation_col <- data.frame(
  Direction = hotspots$dir.Methy[column_order],
  row.names = colnames(heatmap_matrix_tf)
)

# Row annotation
tf_family <- c(
  "Ar"         = "Nuclear receptor",
  "Zfx"        = "C2H2 ZFX",
  "Neurod2"    = "bHLH",
  "Olig2"      = "bHLH",
  "Tcf12"      = "bHLH",
  "Twist2"     = "bHLH",
  "Zic1::Zic2" = "C2H2 ZIC",
  "Zic3"       = "C2H2 ZIC"
)

heatmap_tf_names <- TFresults_report$tf_name

heatmap_families <- unname(
  tf_family[heatmap_tf_names]
)

heatmap_families[
  is.na(heatmap_families)
] <- "Other"

annotation_row <- data.frame(
  Family = heatmap_families,
  row.names = rownames(heatmap_matrix_tf),
  stringsAsFactors = FALSE
)

annotation_colors <- list(
  Direction = c(
    hyper = "deeppink3",
    hypo = "dodgerblue4"
  ),
  Family = structure(
    hcl.colors(
      length(unique(annotation_row$Family)),
      palette = "Set 3"
    ),
    names = unique(annotation_row$Family)
  )
)

pheatmap(
  heatmap_matrix_tf,
  color = c("floralwhite", "gold"),
  breaks = c(-0.5, 0.5, 1.5),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  border_color = "grey80",
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 8,
  main = "Exploratory TF motif candidates across cocaine hotspots",
  legend = FALSE
)
