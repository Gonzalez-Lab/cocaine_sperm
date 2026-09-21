################################################################################
# HISTONE ChIP-seq SIGNAL AT COCAINE-ASSOCIATED DMRs
#
# RRBS-RESTRICTED RANDOMIZATION ANALYSIS
#
# Biological question:
# Do the 24 cocaine-associated DMRs show greater sperm histone-modification
# ChIP-seq signal than expected from RRBS-accessible regions of similar length?
#
# Design:
#   - 24 real DMRs defined in mm39
#   - 10,000 requested random RRBS sets, 24 regions per set
#   - RRBS-accessible regions matched to each DMR by length (+/-20%)
#   - real DMRs excluded; no matching by chromosome, CpG count or CpG density
#   - only complete sets with unique mm39 -> mm9 mappings retained
#   - histone ChIP-seq WIG tracks interrogated in mm9
#
# Regional ChIP signal:
#   1) base-weighted mean signal over each region, separately for rep1 and rep2
#   2) signal_mean = mean(rep1, rep2)
#
# Primary per-mark statistic:
#   mean(log2(signal_mean + 1)) across the 24 regions
#
# Global statistic:
#   - standardize each histone mark relative to its random distribution
#   - calculate mean Z across the 7 histone modifications
#   - compare observed global Z with the retained RRBS random-set values
#
#-----------------------------------------------------
# IMPORTANT:
# The slow WIG-processing step is cached.
# Once "histone_RRBS_random_signal_cache.rds" exists, use:
#
#   RUN_WIG_PROCESSING <- FALSE
#-----------------------------------------------------
################################################################################

# ==============================================================================
# 1. LIBRARIES
# ==============================================================================

library(rtracklayer)
library(GenomicRanges)
library(IRanges)
library(R.utils)
library(GenomeInfoDb)
library(data.table)

# ==============================================================================
# 2. PARAMETERS
# ==============================================================================

hotspot_file <- "hotspots.csv"

data_directory <- "histone data"

chain_gz <- "mm39ToMm9.over.chain.gz"
chain_file <- "mm39ToMm9.over.chain"

standard_chr <- paste0("chr", c(1:19, "X", "Y"))

n_random_sets_requested <- 10000

random_seed <- 12345

# RRBS universe: regions are matched to each DMR only by length (+/-20%).
# Chromosome, CpG number and CpG density are NOT used for matching.
length_tolerance <- 0.20

cx_files <- c(
  "control_1_CX_report.txt.gz",
  "control_2_CX_report.txt.gz",
  "cocaine_1_CX_report.txt.gz",
  "cocaine_2_CX_report.txt.gz"
)

# Set TRUE only if the random regions and WIG signals need to be generated again.
# that part takes a lot of processing time
RUN_WIG_PROCESSING <- TRUE

random_regions_file <- "RRBS_random_regions_10000_requested_sets_mm39_mm9.rds"

signal_cache_file <- "histone_RRBS_random_signal_cache.rds"

# ==============================================================================
# 3. HISTONE FILES
# ==============================================================================

files_chip <- list(
  
  H3K9ac = c("GSM2088389_SPERM_H3K9ac.wig.gz",
             "GSM2401437_SPERM_H3K9ac_replicate2.wig.gz"),
  
  H3K27ac = c("GSM2088387_SPERM_H3K27AC.wig.gz",
              "GSM2401435_SPERM_H3K27AC_replicate2.wig.gz"),
  
  H3K4me1 = c("GSM2088390_SPERM_H3K4me1.wig.gz",
              "GSM2401438_SPERM_H3K4me1_replicate2.wig.gz"),
  
  H3K4me3 = c("GSM2088391_Sperm_H3K4me3.wig.gz",
              "GSM2401439_SPERM_H3K4me3_replicate2.wig.gz"),
  
  H3K36me3 = c("GSM2088385_SPERM_H3K36me3.wig.gz",
               "GSM2401433_SPERM_H3K36me3_replicate2.wig.gz"),
  
  H3K27me3 = c("GSM2088386_SPERM_H3K27me3.wig.gz",
               "GSM2401434_SPERM_H3K27me3_replicate2.wig.gz"),
  
  H3K9me3 = c("GSM2088388_SPERM_H3K9me3.wig.gz",
              "GSM2401436_SPERM_H3K9me3_replicate2.wig.gz")
  )

paths_chip <- lapply(files_chip, function(x) file.path(data_directory,x))

#check
if (!all( vapply( paths_chip, function(x) all(file.exists(x)), logical(1)))) {
  stop("One or more histone WIG files are missing.")}

# ==============================================================================
# 4. READ THE 24 REAL DMRs
# ==============================================================================

hotspots <- read.csv2(hotspot_file, stringsAsFactors = FALSE)

hotspots$chr <- as.character(hotspots$chr)

hotspots$chr <- ifelse(grepl("^chr", hotspots$chr), hotspots$chr,
                       paste0("chr", hotspots$chr))

hotspots$start <- as.integer(hotspots$start)
hotspots$end <- as.integer(hotspots$end)

hotspots$hotspot <- as.character(hotspots$hotspot)

hotspots_gr <- GRanges(seqnames = hotspots$chr,
                       ranges = IRanges(start=hotspots$start, end=hotspots$end),
                       hotspot_id = hotspots$hotspot)

hotspot_lengths <- width(hotspots_gr)

cat("Number of real DMRs:", length(hotspots_gr),"\n")
cat("DMR lengths:\n")
print(hotspot_lengths)

# ==============================================================================
# 5. mm39 -> mm9 CHAIN
# ==============================================================================

if (!file.exists(chain_file)) { if (!file.exists(chain_gz)) {
      download.file(
      "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/liftOver/mm39ToMm9.over.chain.gz",
      chain_gz, mode = "wb")}
  
  gunzip(chain_gz, destname = chain_file, remove = FALSE, overwrite = TRUE)
}

chain_mm39_to_mm9 <- import.chain(chain_file)

# ==============================================================================
# 6. LIFT COCAINE DMRs TO mm9
# ==============================================================================

hotspots_mm9_list <- liftOver(hotspots_gr,  chain_mm39_to_mm9)

unique_lift <- elementNROWS(hotspots_mm9_list) == 1

if (!all(unique_lift)) {
  stop("Not all 24 real DMRs have a unique mm39 -> mm9 lift.") }

hotspots_mm9 <- unlist( hotspots_mm9_list,  use.names = FALSE)

if (!all( as.character(seqnames(hotspots_mm9)) %in% standard_chr)) {
    stop("At least one real DMR lifted outside the standard mm9 chromosomes.") }

hotspots_mm9$hotspot_id <- hotspots_gr$hotspot_id

cat("All real DMRs have a unique mm39 -> mm9 lift.\n")

# ==============================================================================
# 7. BUILD RRBS UNIVERSE AND GENERATE LENGTH-MATCHED RANDOM SETS
#
# This is the ONLY conceptual difference from the genome-wide analysis.
#
# The null regions are restricted to regions accessible to the RRBS experiment.
# RRBS regions are reconstructed from the four experimental CX reports by:
#   - retaining CpGs on standard chromosomes
#   - joining consecutive CpGs separated by <=100 bp
#   - retaining regions with >=3 CpGs and length >=50 bp
#
# For each real DMR, eligible RRBS regions are matched ONLY by length (+/-20%).
# Chromosome, CpG number and CpG density are NOT used for matching.
# Real DMRs are excluded from the background.
#
# 10,000 sets are requested. After mm39 -> mm9 liftOver, only complete sets are
# retained: all 24 regions must map uniquely, to standard mm9 chromosomes, and
# preserve their original RRBS-region length.
# ==============================================================================

if (!file.exists(random_regions_file)) {

  set.seed(random_seed)

  cat("Building RRBS-accessible universe from CX reports...\n")

  if (!all(file.exists(cx_files))) {
    stop("One or more CX report files are missing.")
  }

  cx <- rbindlist(
    lapply(cx_files, fread),
    use.names = TRUE,
    fill = TRUE
  )

  # Bismark CX reports: first six columns = chr, pos, strand, methylated,
  # unmethylated, context. Only chromosome and CpG position are needed here.
  setnames(
    cx,
    old = names(cx)[1:6],
    new = c("chr", "pos", "strand", "meth_count", "unmeth_count", "context")
  )

  cx[, chr := as.character(chr)]
  cx[!grepl("^chr", chr), chr := paste0("chr", chr)]
  cx[, pos := as.integer(pos)]

  cx_cpg <- unique(
    cx[
      chr %in% standard_chr &
      context == "CG" &
      is.finite(pos),
      .(chr, pos)
    ]
  )

  setorder(cx_cpg, chr, pos)

  cx_cpg[, gap := c(NA_integer_, diff(pos)), by = chr]
  cx_cpg[, new_region := is.na(gap) | gap > 100L]
  cx_cpg[, region_number := cumsum(new_region), by = chr]

  rrbs_universe <- cx_cpg[, .(
    start = min(pos),
    end = max(pos),
    n_cpg = .N
  ), by = .(chr, region_number)]

  rrbs_universe[, length_bp := end - start + 1L]
  rrbs_universe <- rrbs_universe[n_cpg >= 3L & length_bp >= 50L]
  rrbs_universe[, region_id := paste0(chr, "_RRBS_", region_number)]

  cat("RRBS universe regions:", nrow(rrbs_universe), "\n")

  rrbs_gr <- GRanges(
    seqnames = rrbs_universe$chr,
    ranges = IRanges(rrbs_universe$start, rrbs_universe$end),
    region_id = rrbs_universe$region_id
  )

  # Exclude all RRBS regions overlapping a real cocaine-associated DMR.
  keep_background <- !overlapsAny(
    rrbs_gr,
    hotspots_gr,
    ignore.strand = TRUE
  )

  rrbs_background <- rrbs_universe[keep_background]

  cat("RRBS background after excluding real DMRs:",
      nrow(rrbs_background), "\n")

  # Hotspot-specific candidate pools matched ONLY by region length.
  candidate_pools <- vector("list", length(hotspot_lengths))

  for (i in seq_along(hotspot_lengths)) {

    target_length <- hotspot_lengths[i]

    eligible <- rrbs_background[
      length_bp >= target_length * (1 - length_tolerance) &
      length_bp <= target_length * (1 + length_tolerance)
    ]

    if (nrow(eligible) == 0L) {
      stop("No length-matched RRBS candidates for hotspot ",
           hotspots$hotspot[i])
    }

    candidate_pools[[i]] <- eligible
  }

  cat("Candidate pool sizes per DMR:\n")
  print(vapply(candidate_pools, nrow, integer(1)))

  # Generate the requested RRBS random sets in mm39.
  # One region is sampled independently for each DMR-specific length pool.
  proposed_sets <- vector("list", n_random_sets_requested)

  for (set_id in seq_len(n_random_sets_requested)) {

    if (set_id %% 500 == 0) {
      cat("Sampling RRBS set", set_id, "of", n_random_sets_requested, "\n")
    }

    one_set <- vector("list", length(hotspot_lengths))

    for (dmr_index in seq_along(hotspot_lengths)) {

      pool <- candidate_pools[[dmr_index]]
      selected <- pool[sample.int(nrow(pool), 1L)]

      gr <- GRanges(
        seqnames = selected$chr,
        ranges = IRanges(selected$start, selected$end)
      )

      gr$set_id <- set_id
      gr$dmr_index <- dmr_index
      gr$region_id <- selected$region_id

      one_set[[dmr_index]] <- gr
    }

    proposed_sets[[set_id]] <- do.call(c, one_set)
  }

  proposed_mm39 <- do.call(c, proposed_sets)

  # Lift all proposed RRBS regions.
  lifted_list <- liftOver(proposed_mm39, chain_mm39_to_mm9)

  unique_lift <- elementNROWS(lifted_list) == 1L

  lifted_unique <- unlist(
    lifted_list[unique_lift],
    use.names = FALSE
  )

  original_unique <- proposed_mm39[unique_lift]

  valid_unique <-
    as.character(seqnames(lifted_unique)) %in% standard_chr &
    width(lifted_unique) == width(original_unique)

  valid_original_indices <- which(unique_lift)[valid_unique]

  region_valid <- rep(FALSE, length(proposed_mm39))
  region_valid[valid_original_indices] <- TRUE

  # Retain only complete 24-region sets.
  set_id_vector <- proposed_mm39$set_id

  complete_set_ids <- as.integer(names(which(
    tapply(region_valid, set_id_vector, all)
  )))

  if (!length(complete_set_ids)) {
    stop("No complete RRBS random sets survived liftOver QC.")
  }

  keep_mm39 <- proposed_mm39$set_id %in% complete_set_ids
  random_mm39 <- proposed_mm39[keep_mm39]

  # Re-lift the retained complete sets so mm39 and mm9 objects have exactly
  # the same ordering: set 1 region 1...24, set 2 region 1...24, etc.
  retained_lift <- liftOver(random_mm39, chain_mm39_to_mm9)
  stopifnot(all(elementNROWS(retained_lift) == 1L))

  random_mm9 <- unlist(retained_lift, use.names = FALSE)

  stopifnot(
    all(as.character(seqnames(random_mm9)) %in% standard_chr),
    all(width(random_mm9) == width(random_mm39))
  )

  # Renumber retained sets consecutively.
  old_set_ids <- unique(random_mm39$set_id)
  new_set_ids <- match(random_mm39$set_id, old_set_ids)

  random_mm39$set_id <- new_set_ids
  random_mm9$set_id <- new_set_ids
  random_mm9$dmr_index <- random_mm39$dmr_index
  random_mm9$region_id <- random_mm39$region_id

  n_random_sets <- length(old_set_ids)

  cat("Requested RRBS sets:", n_random_sets_requested, "\n")
  cat("Complete RRBS sets retained after liftOver QC:", n_random_sets, "\n")

  saveRDS(
    list(
      random_mm39 = random_mm39,
      random_mm9 = random_mm9,
      n_random_sets_requested = n_random_sets_requested,
      n_random_sets = n_random_sets,
      random_seed = random_seed,
      length_tolerance = length_tolerance
    ),
    random_regions_file
  )

  rm(cx, cx_cpg, proposed_sets, proposed_mm39, lifted_list)
  gc()

} else {

  random_regions <- readRDS(random_regions_file)

  random_mm39 <- random_regions$random_mm39
  random_mm9 <- random_regions$random_mm9
  n_random_sets <- random_regions$n_random_sets

  cat("Existing RRBS random regions loaded from disk.\n")
}

# ==============================================================================
# 8. QC OF RRBS RANDOM SETS
# ==============================================================================

stopifnot(
  length(random_mm9) == n_random_sets * length(hotspots_gr),
  length(random_mm39) == length(random_mm9),
  all(random_mm9$dmr_index == rep(seq_along(hotspot_lengths), n_random_sets))
)

# Each random RRBS region must be within +/-20% of the corresponding DMR length.
random_width_matrix <- matrix(
  width(random_mm39),
  nrow = length(hotspots_gr),
  ncol = n_random_sets
)

min_width <- hotspot_lengths * (1 - length_tolerance)
max_width <- hotspot_lengths * (1 + length_tolerance)

stopifnot(
  all(random_width_matrix >= matrix(min_width,
                                    nrow = length(hotspots_gr),
                                    ncol = n_random_sets)),
  all(random_width_matrix <= matrix(max_width,
                                    nrow = length(hotspots_gr),
                                    ncol = n_random_sets))
)

cat("RRBS random control contains", n_random_sets, "complete sets x",
    length(hotspots_gr), "regions.\n")
cat("Every random region is RRBS-accessible and length-matched (+/-20%).\n")

# ==============================================================================
# 9. FUNCTIONS FOR REGIONAL ChIP SIGNAL
# ==============================================================================

clean_track <- function(track) {
  
  chromosome_names <- seqlevels(track)
  
  cleaned_names <- sub("^.*\\.(chr[0-9]+|chrX|chrY|chrM)$", "\\1",
                       chromosome_names)
  
  if (!anyDuplicated(cleaned_names) ) { seqlevels(track) <- cleaned_names }
  
  track <- track[ as.character(seqnames(track)) %in% standard_chr ]
  
  if (!"score" %in% colnames(mcols(track))) {
        stop("The imported WIG file does not contain a score column." )  }
  
  track$score <- as.numeric(track$score)
  
  track <- track[ !is.na(track$score)]
  
  track
}

signal_over_regions <- function( regions, track) {
  
  signal <- numeric(length(regions))
  
  hits <- findOverlaps(regions, track, ignore.strand = TRUE )
  
  if (length(hits) == 0 ) { return(signal) }
  
  overlap_ranges <- pintersect( regions[queryHits(hits)], track[subjectHits(hits)],
                                ignore.strand = TRUE)
  
  weighted_signal <- width(overlap_ranges) * track$score[subjectHits(hits)]
  
  signal_sum <- rowsum(weighted_signal, group=queryHits(hits), reorder = FALSE)
  
  indices <- as.integer(rownames(signal_sum))
  
  signal[indices] <- signal_sum[, 1] / width(regions)[indices]
  
  signal
}

# ==============================================================================
# 10. PROCESS HISTONE WIG FILES
#
# SLOW STEP.
#
# For each histone mark:
#   - calculate regional signal in the 24 DMRs
#   - calculate regional signal in all 240,000 random regions
#   - do this independently for the two biological replicates
#   - average the two replicate regional signals
#
# Results are cached after every histone mark.
# ==============================================================================

if (RUN_WIG_PROCESSING) { signal_cache <- list()
  
  for ( mark in names(paths_chip) ) {
    
    cat("\n========================================\n")
    cat("Processing", mark,"\n")
    cat("========================================\n")
    
    # --------------------------------------------------------------------------
    # Replicate 1
    # --------------------------------------------------------------------------
    
    cat("Importing replicate 1...\n")
    
    track_1 <- clean_track(  import( paths_chip[[mark]][1] )  )
    
    hotspot_rep1 <- signal_over_regions( hotspots_mm9, track_1  )
    
    random_rep1 <- signal_over_regions( random_mm9, track_1 )
    
    rm(track_1)
    
    gc()
    
    # --------------------------------------------------------------------------
    # Replicate 2
    # --------------------------------------------------------------------------
    
    cat( "Importing replicate 2...\n" )
    
    track_2 <- clean_track( import( paths_chip[[mark]][2] ) )
    
    hotspot_rep2 <- signal_over_regions( hotspots_mm9 , track_2 )
    
    random_rep2 <- signal_over_regions( random_mm9 , track_2)
    
    rm(track_2)
    
    gc()
    
    # --------------------------------------------------------------------------
    # Mean signal across biological replicates
    # --------------------------------------------------------------------------
    
    hotspot_signal_mean <- rowMeans(  cbind(hotspot_rep1, hotspot_rep2  ) )
    
    random_signal_mean <- rowMeans( cbind( random_rep1, random_rep2  ) )
    
    signal_cache[[mark]] <- list(
      
      hotspot_rep1 =  hotspot_rep1,
      hotspot_rep2 =  hotspot_rep2,
      hotspot_signal_mean = hotspot_signal_mean,
      random_rep1 = random_rep1,
      random_rep2 = random_rep2,
      random_signal_mean = random_signal_mean
    )
    
    # Save after every mark.
    saveRDS(signal_cache, signal_cache_file )
    
    cat("Saved", mark, "to", signal_cache_file, "\n")
  }
  
  
} else {
  
  if (!file.exists(signal_cache_file) ) {
    stop(signal_cache_file,
      " does not exist. Run once with RUN_WIG_PROCESSING <- TRUE.")
  }
  
  signal_cache <- readRDS(signal_cache_file )
  
  cat("Signal cache loaded. WIG files were NOT imported.\n")
}

# ==============================================================================
# 11. CHECK CACHE
# ==============================================================================

signal_cache <- readRDS(signal_cache_file)

marks <- names(signal_cache)

n_dmr <- length( signal_cache[[1]]$hotspot_signal_mean )

n_random_regions <- length( signal_cache[[1]]$random_signal_mean )

n_random_sets_from_cache <-  n_random_regions / n_dmr

stopifnot(  n_dmr == 24 )

stopifnot( n_random_sets_from_cache == n_random_sets)

stopifnot( length(marks) == 7 )

cat( "\nDMRs:", n_dmr,  "\n" )
cat( "Random sets:",  n_random_sets_from_cache, "\n")
cat( "Histone marks:", paste(marks, collapse = ", " ), "\n")


# ==============================================================================
# 12. DESCRIPTIVE CONTINUOUS SIGNAL
# ==============================================================================

descriptive_results <-  data.frame()

for ( mark in marks) {
  
  dmr <- signal_cache[[mark]]$hotspot_signal_mean
  
  random <-  signal_cache[[mark]]$random_signal_mean
  
  temp <- data.frame(
    
    histone_mark = mark,
    
    DMR_zero_prop =  mean(dmr == 0),
    
    random_zero_prop = mean(random == 0),
    
    DMR_positive_prop = mean(dmr > 0),
    
    random_positive_prop =  mean(random > 0),
    
    DMR_mean =  mean(dmr),
    
    random_mean =  mean(random),
    
    DMR_median =  median(dmr),
    
    random_median = median(random),
    
    DMR_Q75 =  as.numeric( quantile(dmr,  0.75  )  ),
    
    random_Q75 =  as.numeric( quantile(random, 0.75 )  ),
    
    DMR_Q90 =  as.numeric(quantile(dmr,0.90 ) ),
    
    random_Q90 = as.numeric(quantile(random, 0.90 ) ),
    
    DMR_Q95 = as.numeric(quantile(dmr,0.95 ) ),
    
    random_Q95 = as.numeric(quantile(random,0.95 ) ),
    
    DMR_mean_log = mean(  log2(dmr + 1)   ),
    
    random_mean_log =  mean( log2(random + 1)  )
  )
  
  
  descriptive_results <- rbind( descriptive_results, temp )
}


cat(  "\n============================================================\n")
cat("DESCRIPTIVE CONTINUOUS ChIP SIGNAL\n")
cat("============================================================\n\n")

print(descriptive_results)

write.csv2(descriptive_results, "histone RRBS continuous signal stats.csv")


# ==============================================================================
# 13. PRIMARY PER-MARK RANDOMIZATION TEST
#
# Statistic:
#
#   mean(log2(signal_mean + 1))
#
# calculated across the 24 DMRs.
#
# The exact same statistic is calculated for every one of the 10,000
# RRBS random sets.
#
# One-sided hypothesis:
#
#   H1 = ChIP signal is greater at the DMRs than expected from random
#        length-matched RRBS-accessible regions.
# ==============================================================================

test_results <- data.frame()

random_distributions <-  list()

for (mark in marks) {
  
  cat("\nTesting",  mark, "\n")
  
  dmr <- signal_cache[[mark]]$hotspot_signal_mean
  
  random <- signal_cache[[mark]]$random_signal_mean
  
  # ---------------------------------------------------------------------------
  # Observed statistic
  # ---------------------------------------------------------------------------
  
  observed_stat <- mean(log2(dmr + 1) )
  
  # ---------------------------------------------------------------------------
  # Reconstruct the 10,000 random sets
  # ---------------------------------------------------------------------------
  
  random_matrix <- matrix(random, nrow = n_dmr, ncol = n_random_sets)
  
  # ---------------------------------------------------------------------------
  # Calculate statistic for each random set
  # ---------------------------------------------------------------------------
  
  random_stats <- colMeans( log2(random_matrix + 1) )
  
  random_distributions[[mark]] <- random_stats
  
  # ---------------------------------------------------------------------------
  # One-sided empirical P
  # ---------------------------------------------------------------------------
  
  empirical_p <- (1 + sum(random_stats >= observed_stat))/(n_random_sets + 1)
  
  # ---------------------------------------------------------------------------
  # Difference from random expectation
  # ---------------------------------------------------------------------------
  
  random_mean_stat <-  mean(random_stats)
  
  difference <- observed_stat - random_mean_stat
  
  temp <- data.frame(
    
    histone_mark =  mark,
    
    observed_mean_log_signal = observed_stat,
    
    random_mean_log_signal = random_mean_stat,
    
    difference = difference,
    
    random_median = median(random_stats),
    
    random_Q025 =  as.numeric(quantile( random_stats, 0.025 ) ),
    
    random_Q975 = as.numeric( quantile( random_stats, 0.975 ) ),
    
    empirical_p = empirical_p )
  
  test_results <- rbind( test_results,  temp )
}

# ==============================================================================
# 14. MULTIPLE-TEST CORRECTION ACROSS THE 7 HISTONE MARKS
# ==============================================================================

test_results$BH_adjusted_p <- p.adjust(test_results$empirical_p, method = "BH")

cat("\n============================================================\n")
cat("PER-MARK CONTINUOUS RANDOMIZATION TESTS\n")
cat("============================================================\n\n")

print(test_results)

write.csv2(test_results,"histone continuous signal randomization test.csv")

# ==============================================================================
# 15. PER-MARK NULL DISTRIBUTION PLOTS
# ==============================================================================

for (mark in marks) {
  
  random_stats <- random_distributions[[mark]]
  
  observed_stat <- test_results[test_results$histone_mark == mark,
    "observed_mean_log_signal"]
  
  hist(random_stats,
       breaks = 50,
       main = paste(mark, "- random genomic background"),
       xlab = "Mean log2(ChIP signal + 1) in 24 regions" )
  
  abline(v = observed_stat, col = "red", lwd = 3)
}

################################################################################
################################################################################
#
# GLOBAL SEVEN-HISTONE RANDOMIZATION TEST
#
################################################################################
################################################################################

# ==============================================================================
# 16. BUILD THE 10,000 x 7 RANDOM STATISTIC MATRIX
# ==============================================================================

n_marks <- length(marks)

random_stat_matrix <- matrix( NA_real_, nrow = n_random_sets, ncol = n_marks )

colnames( random_stat_matrix ) <- marks

observed_stats <- numeric( n_marks )

names( observed_stats ) <- marks

for ( i in seq_along(marks) ) {
  
  mark <- marks[i]
  
  dmr <- signal_cache[[mark]]$hotspot_signal_mean
  
  random <- signal_cache[[mark]]$random_signal_mean
  
  # Observed statistic.
  
  observed_stats[i] <- mean(log2( dmr + 1 ) )
  
  # Reconstruct random sets.
  
  random_matrix <- matrix(random, nrow = n_dmr, ncol = n_random_sets )
  
  # Statistic for each random set.
  
  random_stat_matrix[, i] <- colMeans( log2( random_matrix + 1 )  )
}

# ==============================================================================
# 17. STANDARDIZE EACH HISTONE RELATIVE TO ITS OWN RANDOM DISTRIBUTION
#
# Means and SDs are calculated ONLY from the retained RRBS random sets.
# ==============================================================================

random_means <- colMeans( random_stat_matrix )

random_sds <- apply( random_stat_matrix,  2,  sd )

stopifnot( all(random_sds > 0) )

# Observed Z for each histone mark.

observed_Z <- ( observed_stats - random_means) / random_sds

# Z for every random set and every histone mark.

random_Z_matrix <- sweep( random_stat_matrix, 2, random_means, "-" )

random_Z_matrix <- sweep( random_Z_matrix,  2, random_sds,  "/" )

# ==============================================================================
# 18. GLOBAL STATISTIC
#-----------------------------------------------------------------
# Mean standardized signal across the seven histone modifications.
#-----------------------------------------------------------------
# Importantly, the SAME random genomic set is retained across the seven
# histone marks. Therefore the dependence among histone marks is preserved
# in the empirical global null distribution.
# ==============================================================================

observed_global_Z <- mean( observed_Z )

random_global_Z <- rowMeans( random_Z_matrix )

# ==============================================================================
# 19. GLOBAL EMPIRICAL RANDOMIZATION TEST
#
# One-sided:
#
# H1 = DMRs have greater overall histone-modification ChIP-seq signal
#      than expected from RRBS-accessible random regions.
# ==============================================================================

global_empirical_p <- ( 1 + sum(random_global_Z >= observed_global_Z)) / 
                            ( n_random_sets + 1)

observed_percentile <- mean( random_global_Z <= observed_global_Z) * 100

# ==============================================================================
# 20. GLOBAL RESULTS
# ==============================================================================

per_mark_global <- data.frame(histone_mark =  marks,
                              observed_mean_log_signal =  observed_stats,
                              random_mean_log_signal = random_means,
                              random_SD = random_sds,
                              observed_Z = observed_Z)


global_result <- data.frame(observed_global_Z =  observed_global_Z,
                            random_mean =  mean(random_global_Z),
                            random_median = median(random_global_Z),
                            random_Q025 = as.numeric(quantile(random_global_Z, 0.025)),
                            random_Q975 = as.numeric(quantile(random_global_Z,0.975)),
                            observed_percentile = observed_percentile,
                            empirical_p = global_empirical_p)


cat("\nGLOBAL HISTONE ChIP SIGNAL TEST\n\n")

cat("Observed Z scores by histone:\n")
print(round(observed_Z, 3))

cat("\nObserved global Z =", round(observed_global_Z, 3),
    "\nEmpirical P =", round(global_empirical_p, 4),
    "\nObserved percentile =", round(observed_percentile, 2), "%\n")

# ==============================================================================
# 21. GLOBAL NULL DISTRIBUTION PLOT
# ==============================================================================

hist( random_global_Z,
      breaks = 50,
      main = "Global histone ChIP signal",
      xlab = "Mean standardized ChIP signal across 7 histone marks")

abline( v = observed_global_Z, col = "red", lwd = 3)

# ==============================================================================
# 22. EXPORT FINAL RESULTS
# ==============================================================================

write.csv2( per_mark_global, "global histone signal per mark Zscores.csv",
  row.names = FALSE
)


write.csv2(global_result, "global_histone_signal_test.csv")

saveRDS(list(descriptive_results = descriptive_results,
             per_mark_tests = test_results,
             random_distributions = random_distributions,
             per_mark_global = per_mark_global,
             observed_global_Z = observed_global_Z,
             random_global_Z = random_global_Z,
             observed_percentile = observed_percentile,
             global_empirical_p = global_empirical_p),
    "FINAL_histone_genome_randomization_results.rds")

# ==============================================================================
# 23. BINARY HEATMAP FOR HISTONES SIGNAL
#
# "High relative histone ChIP signal" is defined relative to the
# length-matched RANDOM GENOMIC BACKGROUND.
#
# For each histone mark:
#
#   cutoff = Q75 of regional signal across all RRBS-accessible random regions
#
# A DMR is classified as:
#
#   1 = signal >= genomic-background Q75
#   0 = signal <  genomic-background Q75
#
# IMPORTANT:
# This binary classification is DESCRIPTIVE ONLY.
# Statistical inference is based on the continuous randomization analysis.
################################################################################

library(pheatmap)

# ==============================================================================
# CALCULATE RRBS-BACKGROUND Q75 FOR EACH HISTONE
# ==============================================================================

RRBS_Q75 <- sapply(marks, function(mark) {
    
    random_signal <-  signal_cache[[mark]]$random_signal_mean
    as.numeric(quantile(random_signal, probs = 0.75, na.rm = TRUE, names = FALSE))
  }  )

RRBS_Q75

# ==============================================================================
# DMR SIGNAL MATRIX
# ==============================================================================

dmr_signal_matrix <- sapply( marks, function(mark) {
    
  signal_cache[[mark]]$hotspot_signal_mean
  } )

rownames(dmr_signal_matrix) <- hotspots$hotspot

colnames(dmr_signal_matrix) <- marks

# ==============================================================================
# BINARIZE DMR SIGNAL USING GENOMIC Q75
# ==============================================================================

dmr_high_matrix <- matrix(0L, nrow = nrow(dmr_signal_matrix),
                              ncol = ncol(dmr_signal_matrix),
                              dimnames = dimnames(dmr_signal_matrix) )

for (mark in marks) {
  
  dmr_high_matrix[, mark] <- as.integer(dmr_signal_matrix[, mark] >=
      RRBS_Q75[mark]  )
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\nNumber of DMRs classified as high signal per histone:\n")

print(colSums(dmr_high_matrix))

cat("\nNumber of histone marks per DMR:\n")

print(rowSums(dmr_high_matrix))

table(rowSums(dmr_high_matrix))

cat("\nDMRs with at least one high-signal histone mark:\n")

# ==============================================================================
# HEATMAP
# ==============================================================================

heatmap_binary <- t(dmr_high_matrix)

hotspot_labels <- paste(hotspots$chr, hotspots$hotspot, hotspots$gene_name, sep = "_")

# Extraer HS de hotspot_labels
hs_labels <- sub( ".*_(HS[0-9]+)_.*","\\1", hotspot_labels)

# Match con las columnas del heatmap
idx <- match(colnames(heatmap_binary),hs_labels)

# Comprobar que todos encontraron match
stopifnot(!any(is.na(idx)))

# Asignar labels correctos
colnames(heatmap_binary) <- hotspot_labels[idx]

# ==============================================================================
#  HOTSPOT ANNOTATION
# ==============================================================================

hotspots$dir.Methy <- tolower(trimws(hotspots$dir.Methy))

annotation_col <- data.frame(Direction = hotspots$dir.Methy,
                             row.names = hotspot_labels)

direction_values <- unique(hotspots$dir.Methy)

annotation_colors <- list(
  Direction = setNames(rep(c("deeppink3", "dodgerblue4"), 
  length.out = length(direction_values)), direction_values)  )


# ==============================================================================
#BINARY HEATMAP
# ==============================================================================

pheatmap(
  heatmap_binary,
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

