################################################################################
# HISTONE ChIP ENRICHMENT IN COCAINE-ASSOCIATED METHYLATION HOTSPOTS
#
# Questions:
#   1) Are cocaine-associated RRBS hotspots enriched for each histone mark?
#   2) Are hotspots enriched for histone-marked chromatin overall?
#
# Design:
#   - Hotspots and RRBS universe are defined in mm39.
#   - Background regions are RRBS-aware and matched to each hotspot by:
#       chromosome, length, CpG number, and CpG density.
#   - Histone ChIP WIG tracks are in mm9.
#   - For EACH replicate of EACH histone mark, "high ChIP signal" is defined
#     BEFORE looking at hotspots, using the genome-wide positive-signal Q75.
#   - A region is positive for a mark only if BOTH biological replicates show
#     mean regional ChIP signal >= that replicate's genome-wide Q75.
#   - Enrichment is tested with 10,000 matched permutations.
#
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(R.utils)
})

# ==============================================================================
# 1. PARAMETERS
# ==============================================================================

hotspot_file   <- "hotspots.csv"
data_directory <- "histone data"

cx_files <- c(
  "control_1_CX_report.txt.gz",
  "control_2_CX_report.txt.gz",
  "cocaine_1_CX_report.txt.gz",
  "cocaine_2_CX_report.txt.gz"
)

chain_gz   <- "mm39ToMm9.over.chain.gz"
chain_file <- "mm39ToMm9.over.chain"

standard_chr <- paste0("chr", c(1:19, "X", "Y"))

# Same RRBS matching used in the TF-motif analysis
length_tolerance  <- 0.15
density_tolerance <- 0.25
cpg_tolerance     <- 2L
top_candidates    <- 200L

n_perm      <- 10000L
random_seed <- 111L

# Histone tracks: two replicates per mark
files_chip <- list(
  H3K9ac = c(
    "GSM2088389_SPERM_H3K9ac.wig.gz",
    "GSM2401437_SPERM_H3K9ac_replicate2.wig.gz"
  ),
  H3K27ac = c(
    "GSM2088387_SPERM_H3K27AC.wig.gz",
    "GSM2401435_SPERM_H3K27AC_replicate2.wig.gz"
  ),
  H3K4me1 = c(
    "GSM2088390_SPERM_H3K4me1.wig.gz",
    "GSM2401438_SPERM_H3K4me1_replicate2.wig.gz"
  ),
  H3K4me3 = c(
    "GSM2088391_Sperm_H3K4me3.wig.gz",
    "GSM2401439_SPERM_H3K4me3_replicate2.wig.gz"
  ),
  H3K36me3 = c(
    "GSM2088385_SPERM_H3K36me3.wig.gz",
    "GSM2401433_SPERM_H3K36me3_replicate2.wig.gz"
  ),
  H3K27me3 = c(
    "GSM2088386_SPERM_H3K27me3.wig.gz",
    "GSM2401434_SPERM_H3K27me3_replicate2.wig.gz"
  ),
  H3K9me3 = c(
    "GSM2088388_SPERM_H3K9me3.wig.gz",
    "GSM2401436_SPERM_H3K9me3_replicate2.wig.gz"
  )
)

# ==============================================================================
# 2. HELPERS
# ==============================================================================

clean_track <- function(track) {
  old <- seqlevels(track)
  new <- sub("^.*\\.(chr[0-9]+|chrX|chrY|chrM)$", "\\1", old)
  if (!anyDuplicated(new)) seqlevels(track) <- new

  track <- track[as.character(seqnames(track)) %in% standard_chr]

  if (!"score" %in% colnames(mcols(track))) {
    stop("Imported WIG has no 'score' column.")
  }

  track$score <- as.numeric(track$score)
  track <- track[!is.na(track$score) & width(track) > 0]
  track
}

# Weighted quantile: WIG intervals can have different widths.
# This estimates the quantile across covered genomic bases, not across WIG rows.
weighted_quantile <- function(x, w, prob = 0.75) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]

  if (!length(x)) return(NA_real_)

  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)

  x[which(cw >= prob)[1]]
}

genome_q75_positive <- function(track) {
  positive <- track[track$score > 0]
  if (!length(positive)) return(Inf)

  weighted_quantile(
    x = positive$score,
    w = width(positive),
    prob = 0.75
  )
}

# Width-weighted mean ChIP signal across each query region.
# Bases not represented by a WIG interval contribute 0, as in the previous script.
signal_over_regions <- function(regions, track) {
  out <- numeric(length(regions))

  hits <- findOverlaps(regions, track, ignore.strand = TRUE)
  if (!length(hits)) return(out)

  ov <- pintersect(
    regions[queryHits(hits)],
    track[subjectHits(hits)],
    ignore.strand = TRUE
  )

  weighted <- width(ov) * track$score[subjectHits(hits)]
  sums <- rowsum(weighted, group = queryHits(hits), reorder = FALSE)
  idx <- as.integer(rownames(sums))
  out[idx] <- sums[, 1] / width(regions)[idx]
  out
}

# Lift mm39 regions to mm9; keep only unique, standard-chromosome mappings.
lift_unique_mm39_to_mm9 <- function(gr, chain) {
  x <- liftOver(gr, chain)
  keep <- elementNROWS(x) == 1L

  ans <- unlist(x[keep], use.names = FALSE)
  ans <- ans[as.character(seqnames(ans)) %in% standard_chr]

  list(
    gr = ans,
    original_index = which(keep)[
      as.character(seqnames(unlist(x[keep], use.names = FALSE))) %in% standard_chr
    ]
  )
}

sample_background_set <- function(pool_ids, max_attempts = 100L) {
  for (attempt in seq_len(max_attempts)) {
    selected <- rep(NA_character_, length(pool_ids))
    used <- character()
    random_order <- sample(seq_along(pool_ids))
    success <- TRUE

    for (i in random_order) {
      available <- setdiff(pool_ids[[i]], used)
      if (!length(available)) {
        success <- FALSE
        break
      }
      selected[i] <- sample(available, 1L)
      used <- c(used, selected[i])
    }

    if (success && !anyNA(selected)) return(selected)
  }

  stop("Could not assign a unique matched background region to every hotspot.")
}

empirical_upper_p <- function(observed, null) {
  (1 + sum(null >= observed)) / (length(null) + 1)
}

# ==============================================================================
# 3. BUILD THE RRBS UNIVERSE FROM EXPERIMENTAL CX REPORTS
# ==============================================================================

message("Building RRBS universe...")

cx <- rbindlist(lapply(cx_files, fread), use.names = TRUE, fill = TRUE)

setnames(
  cx,
  old = names(cx)[1:6],
  new = c("chr", "pos", "strand", "meth_count", "unmeth_count", "context")
)

cx[, chr := as.character(chr)]
cx[!grepl("^chr", chr), chr := paste0("chr", chr)]
cx[chr == "chrMT", chr := "chrM"]

cx[, pos := as.integer(pos)]
cx[, meth_count := as.numeric(meth_count)]
cx[, unmeth_count := as.numeric(unmeth_count)]
cx[, coverage := meth_count + unmeth_count]

cx_cpg <- cx[
  context == "CG" &
    coverage > 0 &
    chr %in% standard_chr &
    strand %in% c("+", "-")
]

# Collapse both cytosines of a CpG dyad to the forward-strand coordinate.
cx_cpg[, cpg_pos := fifelse(strand == "-", pos - 1L, pos)]

cpgs <- unique(cx_cpg[cpg_pos >= 1L, .(chr, pos = cpg_pos)])
setorder(cpgs, chr, pos)

rm(cx, cx_cpg)
gc()

# Same regional logic as the TF analysis / DSS-compatible universe.
cpgs[, gap := c(NA_integer_, diff(pos)), by = chr]
cpgs[, new_region := is.na(gap) | gap > 100L]
cpgs[, region_number := cumsum(new_region), by = chr]

rrbs_universe <- cpgs[, .(
  start = min(pos),
  end   = max(pos),
  n_cpg = .N
), by = .(chr, region_number)]

rrbs_universe[, length_bp := end - start + 1L]
rrbs_universe[, cpg_density := n_cpg / length_bp]

rrbs_universe <- rrbs_universe[n_cpg >= 3L & length_bp >= 50L]
rrbs_universe[, region_id := paste0(chr, "_RRBS_", region_number)]

message("RRBS universe regions: ", nrow(rrbs_universe))

# ==============================================================================
# 4. HOTSPOTS + MATCHED DMR-LIKE BACKGROUND POOLS
# ==============================================================================

hotspots <- fread(hotspot_file, sep = "auto", header = TRUE)
hotspots[, chr := as.character(chr)]
hotspots[!grepl("^chr", chr), chr := paste0("chr", chr)]

hotspots[, length_bp := end - start + 1L]
hotspots[, cpg_density := n_cpg / length_bp]
hotspots[, hotspot := as.character(hotspot)]

hotspots_gr <- GRanges(
  seqnames = hotspots$chr,
  ranges = IRanges(hotspots$start, hotspots$end),
  hotspot_id = hotspots$hotspot
)

rrbs_gr_all <- GRanges(
  seqnames = rrbs_universe$chr,
  ranges = IRanges(rrbs_universe$start, rrbs_universe$end),
  region_id = rrbs_universe$region_id
)

# Real hotspots cannot be used as null regions.
rrbs_background <- rrbs_universe[
  !overlapsAny(rrbs_gr_all, hotspots_gr, ignore.strand = TRUE)
]

candidate_pools <- vector("list", nrow(hotspots))

for (i in seq_len(nrow(hotspots))) {
  h_chr     <- hotspots$chr[i]
  h_length  <- hotspots$length_bp[i]
  h_n_cpg   <- hotspots$n_cpg[i]
  h_density <- hotspots$cpg_density[i]

  candidates <- rrbs_background[
    chr == h_chr &
      length_bp >= h_length * (1 - length_tolerance) &
      length_bp <= h_length * (1 + length_tolerance) &
      abs(n_cpg - h_n_cpg) <= cpg_tolerance &
      cpg_density >= h_density * (1 - density_tolerance) &
      cpg_density <= h_density * (1 + density_tolerance)
  ]

  if (!nrow(candidates)) {
    stop("No matched RRBS candidates for hotspot: ", hotspots$hotspot[i])
  }

  candidates[, match_score :=
    abs(length_bp - h_length) / h_length +
    abs(n_cpg - h_n_cpg) / max(h_n_cpg, 1L) +
    abs(cpg_density - h_density) / max(h_density, 1e-10)
  ]

  setorder(candidates, match_score)
  candidate_pools[[i]] <- copy(candidates[seq_len(min(top_candidates, .N))])
}

candidate_universe <- unique(
  rbindlist(candidate_pools, use.names = TRUE),
  by = "region_id"
)

candidate_gr <- GRanges(
  seqnames = candidate_universe$chr,
  ranges = IRanges(candidate_universe$start, candidate_universe$end),
  region_id = candidate_universe$region_id
)

# Exclude any candidate that overlaps a real hotspot.
candidate_gr <- candidate_gr[
  !overlapsAny(candidate_gr, hotspots_gr, ignore.strand = TRUE)
]

valid_ids <- as.character(candidate_gr$region_id)

candidate_pool_ids <- lapply(
  candidate_pools,
  function(x) intersect(x$region_id, valid_ids)
)

if (any(lengths(candidate_pool_ids) == 0L)) {
  stop("At least one hotspot has no valid matched background candidates.")
}

# ==============================================================================
# 5. LIFT HOTSPOTS AND BACKGROUND mm39 -> mm9
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

chain <- import.chain(chain_file)

# Keep IDs explicitly so matching survives liftOver.
hotspots_gr$source_id  <- hotspots_gr$hotspot_id
candidate_gr$source_id <- candidate_gr$region_id

lift_with_id <- function(gr, chain) {
  lo <- liftOver(gr, chain)
  keep <- elementNROWS(lo) == 1L

  out <- unlist(lo[keep], use.names = FALSE)
  original_ids <- gr$source_id[keep]

  standard <- as.character(seqnames(out)) %in% standard_chr
  out <- out[standard]
  out$source_id <- original_ids[standard]

  out
}

hotspots_mm9   <- lift_with_id(hotspots_gr, chain)
candidate_mm9  <- lift_with_id(candidate_gr, chain)

# The observed test requires every real hotspot to be evaluable.
missing_hotspots <- setdiff(hotspots$hotspot, hotspots_mm9$source_id)
if (length(missing_hotspots)) {
  stop(
    "These real hotspots do not map uniquely mm39->mm9: ",
    paste(missing_hotspots, collapse = ", "),
    ". Decide explicitly how to handle them before testing."
  )
}

# Reorder lifted hotspots to original hotspot order.
hotspots_mm9 <- hotspots_mm9[match(hotspots$hotspot, hotspots_mm9$source_id)]

# Prune matched pools to background regions that also map uniquely to mm9.
valid_mm9_background <- candidate_mm9$source_id

candidate_pool_ids <- lapply(
  candidate_pool_ids,
  function(ids) intersect(ids, valid_mm9_background)
)

if (any(lengths(candidate_pool_ids) == 0L)) {
  stop("At least one hotspot has no matched candidates after mm39->mm9 liftOver.")
}

# Reorder candidate GRanges by ID for fast lookup.
candidate_mm9 <- candidate_mm9[match(valid_mm9_background, candidate_mm9$source_id)]
names(candidate_mm9) <- candidate_mm9$source_id

message("Hotspots evaluable in mm9: ", length(hotspots_mm9))
message("Unique matched background regions evaluable in mm9: ", length(candidate_mm9))

# ==============================================================================
# 6. PRE-GENERATE THE SAME 10,000 MATCHED NULL SETS FOR ALL HISTONE MARKS
# ==============================================================================

set.seed(random_seed)

permutation_ids <- matrix(
  NA_character_,
  nrow = n_perm,
  ncol = nrow(hotspots)
)

for (b in seq_len(n_perm)) {
  permutation_ids[b, ] <- sample_background_set(candidate_pool_ids)

  if (b %% 500L == 0L) {
    message("Matched permutations generated: ", b, " / ", n_perm)
  }
}

# ==============================================================================
# 7. HISTONE MARK ANALYSIS
# ==============================================================================

# For every mark we store:
#   - genome-wide Q75 in replicate 1 and replicate 2
#   - mean regional signal in each replicate
#   - consensus binary call: BOTH replicates >= their own genome-wide Q75

mark_region_calls <- list()
threshold_table   <- list()

for (mark in names(files_chip)) {

  message("\nProcessing ", mark)

  paths <- file.path(data_directory, files_chip[[mark]])
  if (!all(file.exists(paths))) {
    stop("Missing WIG file(s) for ", mark)
  }

  track1 <- clean_track(import(paths[1]))
  track2 <- clean_track(import(paths[2]))

  # CRITICAL: thresholds come from the genome-wide tracks,
  # independently of hotspots and RRBS background.
  q75_rep1 <- genome_q75_positive(track1)
  q75_rep2 <- genome_q75_positive(track2)

  message(
    mark, " genome-wide positive-signal Q75: rep1=",
    signif(q75_rep1, 5), "; rep2=", signif(q75_rep2, 5)
  )

  hs1 <- signal_over_regions(hotspots_mm9, track1)
  hs2 <- signal_over_regions(hotspots_mm9, track2)

  bg1 <- signal_over_regions(candidate_mm9, track1)
  bg2 <- signal_over_regions(candidate_mm9, track2)

  # Conservative replicate-consensus definition.
  hs_high <- (hs1 >= q75_rep1) & (hs2 >= q75_rep2)
  bg_high <- (bg1 >= q75_rep1) & (bg2 >= q75_rep2)

  names(hs_high) <- hotspots_mm9$source_id
  names(bg_high) <- candidate_mm9$source_id

  mark_region_calls[[mark]] <- list(
    hotspot_high = hs_high,
    background_high = bg_high,
    hotspot_signal_rep1 = hs1,
    hotspot_signal_rep2 = hs2,
    background_signal_rep1 = bg1,
    background_signal_rep2 = bg2
  )

  threshold_table[[mark]] <- data.table(
    histone_mark = mark,
    q75_rep1 = q75_rep1,
    q75_rep2 = q75_rep2
  )

  rm(track1, track2)
  gc()
}

threshold_table <- rbindlist(threshold_table)

# ==============================================================================
# 8. QUESTION 1: ENRICHMENT FOR EACH HISTONE MARK
# ==============================================================================

per_mark_results <- vector("list", length(mark_region_calls))
names(per_mark_results) <- names(mark_region_calls)

for (mark in names(mark_region_calls)) {

  hs <- mark_region_calls[[mark]]$hotspot_high
  bg <- mark_region_calls[[mark]]$background_high

  observed_count <- sum(hs)
  observed_fraction <- mean(hs)

  null_counts <- integer(n_perm)

  for (b in seq_len(n_perm)) {
    null_counts[b] <- sum(bg[permutation_ids[b, ]])
  }

  expected_count <- mean(null_counts)
  expected_fraction <- expected_count / length(hs)

  empirical_p <- empirical_upper_p(observed_count, null_counts)

  fold_enrichment <- (observed_count + 0.5) / (expected_count + 0.5)

  per_mark_results[[mark]] <- data.table(
    histone_mark = mark,
    observed_high_DMRs = observed_count,
    total_DMRs = length(hs),
    observed_fraction = observed_fraction,
    expected_high_DMRs = expected_count,
    expected_fraction = expected_fraction,
    fold_enrichment = fold_enrichment,
    empirical_p = empirical_p
  )
}

per_mark_results <- rbindlist(per_mark_results)
per_mark_results[, adjusted_p_BH := p.adjust(empirical_p, method = "BH")]
setorder(per_mark_results, adjusted_p_BH, empirical_p)

# ==============================================================================
# 9. QUESTION 2: INTEGRATED HISTONE-MARKED CHROMATIN
# ==============================================================================

# Matrix: regions x histone marks
hotspot_binary <- do.call(
  cbind,
  lapply(mark_region_calls, `[[`, "hotspot_high")
)

background_binary <- do.call(
  cbind,
  lapply(mark_region_calls, `[[`, "background_high")
)

colnames(hotspot_binary) <- names(mark_region_calls)
colnames(background_binary) <- names(mark_region_calls)

# Integrated endpoint A:
# number of hotspots carrying >=1 assayed high histone mark
hotspot_any_mark <- rowSums(hotspot_binary) >= 1L
background_any_mark <- rowSums(background_binary) >= 1L

observed_any_count <- sum(hotspot_any_mark)
null_any_counts <- integer(n_perm)

# Integrated endpoint B:
# total number of high-mark calls across all hotspot x mark combinations.
# This captures multiplicity of marked states, not just presence/absence.
observed_total_mark_calls <- sum(hotspot_binary)
null_total_mark_calls <- integer(n_perm)

# Optional descriptive endpoint:
# mean number of high marks per region
observed_mean_marks <- mean(rowSums(hotspot_binary))
null_mean_marks <- numeric(n_perm)

for (b in seq_len(n_perm)) {

  ids <- permutation_ids[b, ]
  mat <- background_binary[ids, , drop = FALSE]

  null_any_counts[b] <- sum(rowSums(mat) >= 1L)
  null_total_mark_calls[b] <- sum(mat)
  null_mean_marks[b] <- mean(rowSums(mat))
}

integrated_results <- data.table(
  test = c(
    "Regions with >=1 high histone mark",
    "Total high histone-mark calls",
    "Mean number of high histone marks per region"
  ),
  observed = c(
    observed_any_count,
    observed_total_mark_calls,
    observed_mean_marks
  ),
  expected = c(
    mean(null_any_counts),
    mean(null_total_mark_calls),
    mean(null_mean_marks)
  ),
  empirical_p = c(
    empirical_upper_p(observed_any_count, null_any_counts),
    empirical_upper_p(observed_total_mark_calls, null_total_mark_calls),
    empirical_upper_p(observed_mean_marks, null_mean_marks)
  )
)

# These three integrated statistics are related descriptions of the same
# biological question; report them transparently. Do not present them as three
# independent discoveries.

# ==============================================================================
# 10. EXPORT
# ==============================================================================

write.csv(
  threshold_table,
  "histone_genomewide_Q75_thresholds.csv",
  row.names = FALSE
)

write.csv(
  per_mark_results,
  "histone_enrichment_per_mark.csv",
  row.names = FALSE
)

write.csv(
  integrated_results,
  "histone_enrichment_integrated.csv",
  row.names = FALSE
)

write.csv(
  cbind(
    hotspot_id = rownames(hotspot_binary),
    as.data.frame(hotspot_binary * 1L)
  ),
  "histone_high_signal_hotspot_matrix.csv",
  row.names = FALSE
)

# Save null distributions for reviewer/reproducibility.
saveRDS(
  list(
    permutation_ids = permutation_ids,
    integrated_null_any = null_any_counts,
    integrated_null_total_calls = null_total_mark_calls,
    integrated_null_mean_marks = null_mean_marks
  ),
  "histone_enrichment_null_distributions.rds"
)

# ==============================================================================
# 11. CONSOLE SUMMARY
# ==============================================================================

cat("\n\n============================================================\n")
cat("GENOME-WIDE Q75 THRESHOLDS\n")
cat("============================================================\n")
print(threshold_table)

cat("\n============================================================\n")
cat("PER-MARK RRBS-AWARE ENRICHMENT\n")
cat("============================================================\n")
print(per_mark_results)

cat("\n============================================================\n")
cat("INTEGRATED HISTONE-MARK ANALYSIS\n")
cat("============================================================\n")
print(integrated_results)

cat("\nAnalysis completed successfully.\n")
