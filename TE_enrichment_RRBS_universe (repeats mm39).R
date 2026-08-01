library(GenomicRanges)
library(data.table)

#upload cocaine DMRs hotspots
hotspots <- read.csv2("hotspots.csv")

hotspots_gr <- GRanges(
  seqnames = paste0("chr", hotspots$chr),
  ranges = IRanges(hotspots$start, hotspots$end),
  diffMethy = hotspots$diff.Methy
)

################################################################################
# Transposable Elements TE analysis (repeats)
################################################################################

#download repeats 
repeats <- fread(
  "http://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/rmsk.txt.gz",
  header = FALSE
)

colnames(repeats) <- c(
  "bin", "swScore", "milliDiv", "milliDel", "milliIns",
  "genoName", "genoStart", "genoEnd", "genoLeft",
  "strand",
  "repName", "repClass", "repFamily",
  "repStart", "repEnd", "repLeft",
  "id"
)


#convert to granges
repeats_gr <- GRanges(
  seqnames = repeats$genoName,
  ranges = IRanges(repeats$genoStart + 1, repeats$genoEnd),
  strand = repeats$strand,
  repName = repeats$repName,
  repClass = repeats$repClass,
  repFamily = repeats$repFamily
)

standard_chr <- paste0("chr", c(1:19, "X", "Y"))

repeats_gr <- repeats_gr[seqnames(repeats_gr) %in% standard_chr]

##
LINE <- repeats_gr[repeats_gr$repClass == "LINE"]
SINE <- repeats_gr[repeats_gr$repClass == "SINE"]
LTR  <- repeats_gr[repeats_gr$repClass == "LTR"]

##################
#Repeats analysis
#################

# direct overlaps
ov_LINE <- findOverlaps(hotspots_gr, LINE)
ov_SINE <- findOverlaps(hotspots_gr, SINE)
ov_LTR  <- findOverlaps(hotspots_gr, LTR)

# asign
hotspots_gr$LINE_overlap <- FALSE
hotspots_gr$SINE_overlap <- FALSE
hotspots_gr$LTR_overlap  <- FALSE

hotspots_gr$LINE_overlap[queryHits(ov_LINE)] <- TRUE
hotspots_gr$SINE_overlap[queryHits(ov_SINE)] <- TRUE
hotspots_gr$LTR_overlap[queryHits(ov_LTR)]   <- TRUE

# combined
hotspots_gr$repeat_overlap <-
  hotspots_gr$LINE_overlap |
  hotspots_gr$SINE_overlap |
  hotspots_gr$LTR_overlap

dist_LINE <- distanceToNearest(hotspots_gr, LINE)
dist_SINE <- distanceToNearest(hotspots_gr, SINE)
dist_LTR  <- distanceToNearest(hotspots_gr, LTR)

hotspots_gr$dist_LINE <- mcols(dist_LINE)$distance
hotspots_gr$dist_SINE <- mcols(dist_SINE)$distance
hotspots_gr$dist_LTR  <- mcols(dist_LTR)$distance

#############################################
# Define window treshold
#############################################

threshold <- 1000  # 1 kb

hotspots_gr$near_repeat <-
  hotspots_gr$dist_LINE < threshold |
  hotspots_gr$dist_SINE < threshold |
  hotspots_gr$dist_LTR  < threshold

hotspots_gr$near_LINE <-  hotspots_gr$dist_LINE < threshold
hotspots_gr$near_SINE <-  hotspots_gr$dist_SINE < threshold
hotspots_gr$near_LTR <-  hotspots_gr$dist_LTR < threshold

hotspots_gr$class <- ifelse(hotspots_gr$near_repeat, "near_repeat","no")
hotspots_gr$class[hotspots_gr$repeat_overlap] <- "repeat_direct"


near_cols <- c("near_LINE", "near_SINE", "near_LTR")
te_names  <- c("LINE", "SINE", "LTR")

hotspots_gr$near_TE <- apply(
  as.data.frame(mcols(hotspots_gr)[, near_cols]),
  1,
  function(x) {
    hits <- te_names[which(x)]
    if (length(hits) == 0) return("none")
    paste(hits, collapse = "-")
  }
)

table(hotspots_gr$near_TE)

hotspots$TEassociation <- hotspots_gr$class
hotspots$near_TE <- hotspots_gr$near_TE

##############################################################################

################################################################################
# RRBS-AWARE TRANSPOSABLE ELEMENT ANALYSIS OF COCAINE-ASSOCIATED DMR HOTSPOTS
#
# This script:
#   1. Builds an RRBS-testable genomic universe from four Bismark CX reports.
#   2. Reads the cocaine-associated DMR hotspots and calculates their CpG content.
#   3. Matches each hotspot to comparable RRBS regions from the same chromosome.
#   4. Tests association with LINE, SINE, and LTR elements using 10,000 permutations.
#
# Required input files:
#   - control_1_CX_report.txt.gz
#   - control_2_CX_report.txt.gz
#   - cocaine_1_CX_report.txt.gz
#   - cocaine_2_CX_report.txt.gz
#   - hotspots.csv
#
# Required objects already loaded in the R environment:
#   - LINE: GRanges with mm39 LINE annotations
#   - SINE: GRanges with mm39 SINE annotations
#   - LTR:  GRanges with mm39 LTR annotations
################################################################################


# ==============================================================================
# 1. PACKAGES
# ==============================================================================

library(data.table)
library(GenomicRanges)
library(IRanges)


# ==============================================================================
# 2. INPUT FILES AND ANALYSIS PARAMETERS
# ==============================================================================

cx_files <- c(
  "control_1_CX_report.txt.gz",
  "control_2_CX_report.txt.gz",
  "cocaine_1_CX_report.txt.gz",
  "cocaine_2_CX_report.txt.gz"
)

hotspot_file <- "hotspots.csv"

# Parameters used to build DSS-compatible RRBS regions
max_cpg_gap <- 100L
minimum_cpgs <- 3L
minimum_region_length <- 50L

# TE proximity and coverage parameters
threshold_bp <- 1000L
flank_bp <- 1000L

# Matching parameters
length_tolerance <- 0.15
density_tolerance <- 0.25
cpg_tolerance <- 2L
top_candidates <- 100L

# Permutation parameters
n_perm <- 10000L
random_seed <- 123L

standard_chr <- paste0("chr", c(1:19, "X", "Y"))


# ==============================================================================
# 3. CHECK INPUT FILES AND TE OBJECTS
# ==============================================================================

missing_input_files <- c(cx_files, hotspot_file)[
  !file.exists(c(cx_files, hotspot_file))
]

if (length(missing_input_files) > 0) {
  stop(
    "The following input files were not found:\n",
    paste(missing_input_files, collapse = "\n")
  )
}

if (!exists("LINE") || !is(LINE, "GRanges")) {
  stop("The object LINE does not exist or is not a GRanges object.")
}

if (!exists("SINE") || !is(SINE, "GRanges")) {
  stop("The object SINE does not exist or is not a GRanges object.")
}

if (!exists("LTR") || !is(LTR, "GRanges")) {
  stop("The object LTR does not exist or is not a GRanges object.")
}


# ==============================================================================
# 4. READ AND COMBINE THE FOUR CX REPORTS
# ==============================================================================

message("Reading CX reports...")

cx1 <- fread(cx_files[1], header = FALSE)
cx2 <- fread(cx_files[2], header = FALSE)
cx3 <- fread(cx_files[3], header = FALSE)
cx4 <- fread(cx_files[4], header = FALSE)

cx <- rbindlist(
  list(cx1, cx2, cx3, cx4),
  use.names = TRUE,
  fill = TRUE
)

# The four temporary sample-specific objects are no longer needed.
rm(cx1, cx2, cx3, cx4)
invisible(gc())

if (ncol(cx) < 6) {
  stop("The CX reports contain fewer than six columns.")
}

# Standard Bismark CX report columns:
# V1 chromosome, V2 position, V3 strand, V4 methylated reads,
# V5 unmethylated reads, and V6 sequence context.
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

cx[, coverage := meth_count + unmeth_count]

# Retain only CpG cytosines with actual sequencing coverage.
cx_cpg <- cx[
  context == "CG" & coverage > 0
]

cat("Covered CpG entries:", nrow(cx_cpg), "\n")

cx_cpg[, chr := as.character(chr)]
cx_cpg[, pos := as.integer(pos)]
cx_cpg[!grepl("^chr", chr), chr := paste0("chr", chr)]

# Keep only unique covered CpG coordinates across all four libraries.
cpgs <- unique(cx_cpg[, .(chr, pos)])
setorder(cpgs, chr, pos)

cat("Unique covered CpG coordinates:", nrow(cpgs), "\n")

# The full combined CX table is no longer required after extracting CpG positions.
rm(cx, cx_cpg)
invisible(gc())


# ==============================================================================
# 5. BUILD THE RRBS-TESTABLE GENOMIC UNIVERSE
# ==============================================================================

# Start a new region whenever consecutive CpGs are separated by more than 100 bp.
cpgs[, gap := c(NA_integer_, diff(pos)), by = chr]
cpgs[, new_region := is.na(gap) | gap > max_cpg_gap]
cpgs[, region_number := cumsum(new_region), by = chr]

rrbs_universe <- cpgs[
  ,
  .(
    start = min(pos),
    end = max(pos),
    n_cpg = .N
  ),
  by = .(chr, region_number)
]

rrbs_universe[, length_bp := end - start + 1L]
rrbs_universe[, cpg_density := n_cpg / length_bp]

# Apply criteria equivalent to the region-level DSS settings.
rrbs_universe <- rrbs_universe[
  n_cpg >= minimum_cpgs &
    length_bp >= minimum_region_length
]

rrbs_universe[
  ,
  region_id := paste0(chr, "_RRBS_", region_number)
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

setorder(rrbs_universe, chr, start)

cat("RRBS-testable regions:", nrow(rrbs_universe), "\n")


# ==============================================================================
# 6. READ AND PREPARE THE REAL DMR HOTSPOTS
# ==============================================================================

hotspots <- fread(
  hotspot_file,
  sep = ";",
  header = TRUE,
  data.table = TRUE
)

required_hotspot_columns <- c("chr", "start", "end")
missing_hotspot_columns <- setdiff(required_hotspot_columns, names(hotspots))

if (length(missing_hotspot_columns) > 0) {
  stop(
    "Missing columns in hotspots.csv: ",
    paste(missing_hotspot_columns, collapse = ", ")
  )
}

hotspots[, chr := as.character(chr)]
hotspots[!grepl("^chr", chr), chr := paste0("chr", chr)]
hotspots[, start := as.integer(start)]
hotspots[, end := as.integer(end)]

if (anyNA(hotspots[, .(chr, start, end)])) {
  stop("The hotspot coordinates contain missing or non-numeric values.")
}

if (any(hotspots$end < hotspots$start)) {
  stop("At least one hotspot has an end coordinate smaller than its start.")
}

hotspots[, hotspot_id := sprintf("HS_%02d", seq_len(.N))]
hotspots[, length_bp := end - start + 1L]

cpg_gr <- GRanges(
  seqnames = cpgs$chr,
  ranges = IRanges(start = cpgs$pos, end = cpgs$pos)
)

hotspots_gr <- GRanges(
  seqnames = hotspots$chr,
  ranges = IRanges(start = hotspots$start, end = hotspots$end),
  hotspot_id = hotspots$hotspot_id
)

hits <- findOverlaps(
  query = hotspots_gr,
  subject = cpg_gr,
  ignore.strand = TRUE
)

hotspots[, n_cpg := tabulate(queryHits(hits), nbins = .N)]
hotspots[, cpg_density := n_cpg / length_bp]

if (any(hotspots$n_cpg == 0L)) {
  warning(
    "At least one hotspot contains no covered CpGs in the combined CX reports."
  )
}

cat("Real hotspots loaded:", nrow(hotspots), "\n")


# ==============================================================================
# 7. PREPARE TE CLASSES
# ==============================================================================

te_classes <- list(
  LINE = LINE,
  SINE = SINE,
  LTR = LTR
)

te_classes <- lapply(
  te_classes,
  function(te_gr) {
    te_gr[
      as.character(seqnames(te_gr)) %in% standard_chr
    ]
  }
)

if (any(vapply(te_classes, length, integer(1)) == 0L)) {
  stop(
    "At least one TE class became empty after chromosome filtering. ",
    "Confirm that LINE, SINE, and LTR use chr1, chr2, etc."
  )
}


# ==============================================================================
# 8. FUNCTIONS FOR TE ASSOCIATION METRICS
# ==============================================================================

nearest_TE_distance <- function(query_gr, subject_gr) {
  distances <- rep(NA_integer_, length(query_gr))

  hits <- suppressWarnings(
    distanceToNearest(
      query_gr,
      subject_gr,
      ignore.strand = TRUE
    )
  )

  if (length(hits) > 0) {
    distances[queryHits(hits)] <- mcols(hits)$distance
  }

  distances
}


TE_fraction_in_window <- function(query_gr, subject_gr, flank = 1000L) {
  windows_gr <- GRanges(
    seqnames = seqnames(query_gr),
    ranges = IRanges(
      start = pmax(1L, start(query_gr) - flank),
      end = end(query_gr) + flank
    )
  )

  # Merge overlapping elements of the same class to avoid double-counting bases.
  subject_union <- reduce(subject_gr, ignore.strand = TRUE)

  hits <- suppressWarnings(
    findOverlaps(
      windows_gr,
      subject_union,
      ignore.strand = TRUE
    )
  )

  covered_bp <- integer(length(windows_gr))

  if (length(hits) > 0) {
    intersections <- pintersect(
      ranges(windows_gr)[queryHits(hits)],
      ranges(subject_union)[subjectHits(hits)]
    )

    overlap_table <- data.table(
      window_id = queryHits(hits),
      overlap_bp = width(intersections)
    )

    bp_by_window <- overlap_table[
      ,
      .(covered_bp = sum(overlap_bp)),
      by = window_id
    ]

    covered_bp[bp_by_window$window_id] <- bp_by_window$covered_bp
  }

  covered_bp / width(windows_gr)
}


calculate_TE_metrics <- function(
    query_gr,
    subject_gr,
    threshold = 1000L,
    flank = 1000L
) {
  distances <- nearest_TE_distance(query_gr, subject_gr)
  fractions <- TE_fraction_in_window(query_gr, subject_gr, flank)

  data.table(
    distance = distances,
    near_1kb = !is.na(distances) & distances < threshold,
    window_fraction = fractions
  )
}


# ==============================================================================
# 9. CALCULATE OBSERVED TE METRICS
# ==============================================================================

observed_by_hotspot <- data.table(
  hotspot_id = hotspots$hotspot_id,
  chr = hotspots$chr,
  start = hotspots$start,
  end = hotspots$end
)

observed_summary_list <- vector("list", length(te_classes))
names(observed_summary_list) <- names(te_classes)

for (te_name in names(te_classes)) {
  metrics <- calculate_TE_metrics(
    query_gr = hotspots_gr,
    subject_gr = te_classes[[te_name]],
    threshold = threshold_bp,
    flank = flank_bp
  )

  distance_column <- paste0("distance_", te_name)
  near_column <- paste0("near_", te_name)
  fraction_column <- paste0("fraction_", te_name)

  observed_by_hotspot[, (distance_column) := metrics$distance]
  observed_by_hotspot[, (near_column) := metrics$near_1kb]
  observed_by_hotspot[, (fraction_column) := metrics$window_fraction]

  observed_summary_list[[te_name]] <- data.table(
    te_class = te_name,
    n_within_1kb = sum(metrics$near_1kb),
    median_distance = median(metrics$distance, na.rm = TRUE),
    mean_window_fraction = mean(metrics$window_fraction, na.rm = TRUE)
  )
}

observed_summary <- rbindlist(observed_summary_list)

cat("\nOBSERVED RESULTS\n")
cat("================\n")
print(observed_summary)


# ==============================================================================
# 10. BUILD MATCHED RRBS CANDIDATE POOLS
# ==============================================================================

candidate_pools <- vector("list", nrow(hotspots))

for (i in seq_len(nrow(hotspots))) {
  h_chr <- hotspots$chr[i]
  h_length <- hotspots$length_bp[i]
  h_n_cpg <- hotspots$n_cpg[i]
  h_density <- hotspots$cpg_density[i]

  candidates <- rrbs_universe[
    chr == h_chr &
      length_bp >= h_length * (1 - length_tolerance) &
      length_bp <= h_length * (1 + length_tolerance) &
      abs(n_cpg - h_n_cpg) <= cpg_tolerance &
      cpg_density >= h_density * (1 - density_tolerance) &
      cpg_density <= h_density * (1 + density_tolerance)
  ]

  if (nrow(candidates) == 0L) {
    stop("No compatible RRBS candidates were found for ", hotspots$hotspot_id[i], ".")
  }

  candidates[
    ,
    match_score :=
      abs(length_bp - h_length) / h_length +
      abs(n_cpg - h_n_cpg) / h_n_cpg +
      abs(cpg_density - h_density) / h_density
  ]

  setorder(candidates, match_score)
  number_to_keep <- min(top_candidates, nrow(candidates))
  candidate_pools[[i]] <- copy(candidates[seq_len(number_to_keep)])

  cat(
    hotspots$hotspot_id[i], ":", nrow(candidates),
    "compatible candidates; retained", number_to_keep, "\n"
  )
}


# ==============================================================================
# 11. PRECALCULATE TE METRICS FOR ALL UNIQUE CANDIDATES
# ==============================================================================

candidate_universe <- unique(
  rbindlist(candidate_pools, use.names = TRUE),
  by = "region_id"
)

candidate_gr <- GRanges(
  seqnames = candidate_universe$chr,
  ranges = IRanges(
    start = candidate_universe$start,
    end = candidate_universe$end
  )
)

for (te_name in names(te_classes)) {
  metrics <- calculate_TE_metrics(
    query_gr = candidate_gr,
    subject_gr = te_classes[[te_name]],
    threshold = threshold_bp,
    flank = flank_bp
  )

  candidate_universe[, (paste0("distance_", te_name)) := metrics$distance]
  candidate_universe[, (paste0("near_", te_name)) := metrics$near_1kb]
  candidate_universe[, (paste0("fraction_", te_name)) := metrics$window_fraction]
}

metric_columns <- c(
  "region_id",
  paste0("distance_", names(te_classes)),
  paste0("near_", names(te_classes)),
  paste0("fraction_", names(te_classes))
)

candidate_metrics <- candidate_universe[, ..metric_columns]

for (i in seq_along(candidate_pools)) {
  candidate_pools[[i]] <- merge(
    candidate_pools[[i]],
    candidate_metrics,
    by = "region_id",
    all.x = TRUE,
    sort = FALSE
  )
}

cat("\nUnique candidate regions:", nrow(candidate_universe), "\n")


# ==============================================================================
# 12. SAMPLE ONE COMPLETE SET OF PSEUDO-HOTSPOTS
# ==============================================================================

sample_pseudo_set <- function(candidate_pools) {
  selected_regions <- vector("list", length(candidate_pools))
  used_region_ids <- character()

  for (i in seq_along(candidate_pools)) {
    available <- candidate_pools[[i]][!region_id %in% used_region_ids]

    if (nrow(available) == 0L) {
      stop("No unused candidates remain for hotspot ", i, ".")
    }

    chosen <- available[sample.int(nrow(available), size = 1L)]
    used_region_ids <- c(used_region_ids, chosen$region_id)
    selected_regions[[i]] <- chosen
  }

  rbindlist(selected_regions, use.names = TRUE)
}

set.seed(random_seed)
test_pseudo <- sample_pseudo_set(candidate_pools)

cat("\nEXAMPLE MATCHED PSEUDO-HOTSPOT SET\n")
cat("==================================\n")
print(
  test_pseudo[
    ,
    c(
      "region_id",
      "chr",
      "length_bp",
      "n_cpg",
      "cpg_density",
      "near_LINE",
      "near_SINE",
      "near_LTR"
    ),
    with = FALSE
  ]
)


# ==============================================================================
# 13. RUN THE PERMUTATION TEST
# ==============================================================================

set.seed(random_seed)
permutation_results <- vector("list", n_perm)

for (b in seq_len(n_perm)) {
  pseudo <- sample_pseudo_set(candidate_pools)
  one_permutation <- vector("list", length(te_classes))
  names(one_permutation) <- names(te_classes)

  for (te_name in names(te_classes)) {
    distances <- pseudo[[paste0("distance_", te_name)]]
    near_status <- pseudo[[paste0("near_", te_name)]]
    fractions <- pseudo[[paste0("fraction_", te_name)]]

    one_permutation[[te_name]] <- data.table(
      permutation = b,
      te_class = te_name,
      n_within_1kb = sum(near_status, na.rm = TRUE),
      median_distance = median(distances, na.rm = TRUE),
      mean_window_fraction = mean(fractions, na.rm = TRUE)
    )
  }

  permutation_results[[b]] <- rbindlist(one_permutation, use.names = TRUE)

  if (b %% 500L == 0L) {
    message("Permutations completed: ", b, " / ", n_perm)
  }
}

permutation_results <- rbindlist(permutation_results, use.names = TRUE)


# ==============================================================================
# 14. CALCULATE EMPIRICAL P-VALUES
# ==============================================================================

final_results_list <- vector("list", length(te_classes))
names(final_results_list) <- names(te_classes)

for (te_name in names(te_classes)) {
  observed_class <- observed_summary[te_class == te_name]
  null_class <- permutation_results[te_class == te_name]

  # Main test: more hotspots within 1 kb than expected.
  p_near <- (
    1 + sum(null_class$n_within_1kb >= observed_class$n_within_1kb)
  ) / (n_perm + 1)

  # Complementary test: smaller median distance than expected.
  p_distance <- (
    1 + sum(null_class$median_distance <= observed_class$median_distance)
  ) / (n_perm + 1)

  # Complementary test: greater repeat coverage than expected.
  p_fraction <- (
    1 + sum(
      null_class$mean_window_fraction >= observed_class$mean_window_fraction
    )
  ) / (n_perm + 1)

  final_results_list[[te_name]] <- rbindlist(
    list(
      data.table(
        te_class = te_name,
        metric = "Number within 1 kb",
        observed = as.numeric(observed_class$n_within_1kb),
        expected = mean(null_class$n_within_1kb),
        empirical_p = p_near
      ),
      data.table(
        te_class = te_name,
        metric = "Median nearest distance",
        observed = as.numeric(observed_class$median_distance),
        expected = median(null_class$median_distance),
        empirical_p = p_distance
      ),
      data.table(
        te_class = te_name,
        metric = "Mean fraction in +/-1 kb window",
        observed = as.numeric(observed_class$mean_window_fraction),
        expected = mean(null_class$mean_window_fraction),
        empirical_p = p_fraction
      )
    ),
    use.names = TRUE
  )
}

final_results <- rbindlist(final_results_list, use.names = TRUE)

# Correct across all nine tests.
final_results[
  ,
  p_adjusted_global_BH := p.adjust(empirical_p, method = "BH")
]

# Correct across LINE, SINE, and LTR within each metric.
final_results[
  ,
  p_adjusted_within_metric := p.adjust(empirical_p, method = "BH"),
  by = metric
]

principal_results <- final_results[
  metric == "Number within 1 kb",
  .(
    te_class,
    observed,
    expected,
    empirical_p,
    p_adjusted_within_metric
  )
]

cat("\nFINAL RESULTS BY TE CLASS\n")
cat("=========================\n")
print(final_results)

cat("\nMAIN TEST: TE WITHIN 1 KB\n")
cat("==========================\n")
print(principal_results)


# ==============================================================================
# 15. SAVE OUTPUT FILES
# ==============================================================================

fwrite(
  rrbs_universe,
  file = "RRBS_testable_universe_mm39.tsv",
  sep = "\t"
)

fwrite(
  hotspots,
  file = "hotspots_with_RRBS_CpG_metrics.tsv",
  sep = "\t"
)

fwrite(
  test_pseudo,
  file = "pseudo_hotspots_example.tsv",
  sep = "\t"
)

fwrite(
  observed_by_hotspot,
  file = "hotspots_TE_class_metrics_mm39.tsv",
  sep = "\t"
)

fwrite(
  observed_summary,
  file = "hotspots_TE_class_observed_summary.tsv",
  sep = "\t"
)

fwrite(
  permutation_results,
  file = "RRBS_TE_class_permutations.tsv",
  sep = "\t"
)

fwrite(
  final_results,
  file = "RRBS_TE_class_test_results.tsv",
  sep = "\t"
)

fwrite(
  principal_results,
  file = "RRBS_TE_class_main_results.tsv",
  sep = "\t"
)

message("TE analysis completed successfully.")
