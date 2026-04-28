library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(GenomeInfoDb)
library(pheatmap)

#upload cocaine DMRs hotspots
hotspots <- read.csv2("hotspots.csv")

hotspots_gr <- GRanges(
  seqnames = paste0("chr", hotspots$chr),
  ranges = IRanges(hotspots$start, hotspots$end),
  diffMethy = hotspots$diff.Methy,
  dirMethy = hotspots$Dir.Methy,
  gene = hotspots$gene,
  hotspot = hotspots$hotspot
)

################################################################################
#DOWNLOAD ATAC AND ChIPseq DATA
###############################################################################

# function to download GSM (WIG histones)
download_geo_gsm <- function(gsm, filename, outdir = "data"){
  
  if(!dir.exists(outdir)) dir.create(outdir)
  
  prefix <- paste0(substr(gsm, 1, nchar(gsm)-3), "nnn")
  
  url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/",
    prefix, "/", gsm, "/suppl/", filename
  )
  
  dest <- file.path(outdir, filename)
  
  if(!file.exists(dest)){
    download.file(url, dest, mode = "wb")
  }
  
  return(dest)
}


# function to download GSE (ATAC nucleosomes)
download_geo_gse <- function(gse, filename, outdir = "data"){
  
  if(!dir.exists(outdir)) dir.create(outdir)
  
  prefix <- paste0(substr(gse, 1, nchar(gse)-3), "nnn")
  
  url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/",
    prefix, "/", gse, "/suppl/", filename
  )
  
  dest <- file.path(outdir, filename)
  
  if(!file.exists(dest)){
    download.file(url, dest, mode = "wb")
  }
  
  return(dest)
}

######################################################
#list ChIP files
files_chip <- list(
  H3K36me3 = list(
    gsm = c("GSM2088385", "GSM2401433"),
    files = c("GSM2088385_SPERM_H3K36me3.wig.gz",
              "GSM2401433_SPERM_H3K36me3_replicate2.wig.gz")
  ),
  
  H3K27me3 = list(
    gsm = c("GSM2088386", "GSM2401434"),
    files = c("GSM2088386_SPERM_H3K27me3.wig.gz",
              "GSM2401434_SPERM_H3K27me3_replicate2.wig.gz")
  ),
  
  H3K27ac = list(
    gsm = c("GSM2088387", "GSM2401435"),
    files = c("GSM2088387_SPERM_H3K27AC.wig.gz",
              "GSM2401435_SPERM_H3K27AC_replicate2.wig.gz")
  ),
  
  H3K9me3 = list(
    gsm = c("GSM2088388", "GSM2401436"),
    files = c("GSM2088388_SPERM_H3K9me3.wig.gz",
              "GSM2401436_SPERM_H3K9me3_replicate2.wig.gz")
  ),
  
  H3K9ac = list(
    gsm = c("GSM2088389", "GSM2401437"),
    files = c("GSM2088389_SPERM_H3K9ac.wig.gz",
              "GSM2401437_SPERM_H3K9ac_replicate2.wig.gz")
  ),
  
  H3K4me1 = list(
    gsm = c("GSM2088390", "GSM2401438"),
    files = c("GSM2088390_SPERM_H3K4me1.wig.gz",
              "GSM2401438_SPERM_H3K4me1_replicate2.wig.gz")
  ),
  
  H3K4me3 = list(
    gsm = c("GSM2088391", "GSM2401439"),
    files = c("GSM2088391_Sperm_H3K4me3.wig.gz",
              "GSM2401439_SPERM_H3K4me3_replicate2.wig.gz")
  )
)

#download chip files
paths_chip <- lapply(files_chip, function(x){
  mapply(download_geo_gsm,
         gsm = x$gsm,
         filename = x$files,
         SIMPLIFY = TRUE)
})

paths_chip <- as.list(paths_chip)

################################################################################
################################################################################
#Hotspots overlap analysis with sperm ATAC signal 
################################################################################
################################################################################

#define ATAC
gse_atac <- "GSE134825"

files_atac <- c(
  "GSE134825_SV129.Sp.Nucleosomes.bed.gz",
  "GSE134825_SV129.Sp.Subnucleosome.bed.gz"
)

paths_atac <- sapply(files_atac, function(f){
  download_geo_gse(gse_atac, f)
})

paths_atac <- as.list(paths_atac)

nuc <- import(paths_atac[[1]])
subnuc <- import(paths_atac[[2]])

#########################################
df_nuc <- as.data.frame(nuc)
df_subnuc <- as.data.frame(subnuc)

#convert to GRanges object
gr.nuc <- GRanges(
  seqnames = df_nuc$seqnames,
  ranges = IRanges(df_nuc$start, df_nuc$end),
  score = df_nuc$score
)

gr.subnuc <- GRanges(
  seqnames = df_subnuc$seqnames,
  ranges = IRanges(df_subnuc$start, df_subnuc$end),
  score = df_subnuc$score
)

list_nuc <- list(gr.nuc)
list_subnuc <- list(gr.subnuc)

standard_chr <- paste0("chr", c(1:19, "X", "Y"))

######
df_nuc <- do.call(rbind, lapply(list_nuc, function(gr) {
  df <- as.data.frame(gr)
  df[, c("seqnames", "start", "end")]
}))

df_nuc <- df_nuc %>%
  filter(seqnames %in% standard_chr)

df_nuc <- df_nuc %>%
  arrange(seqnames, start, end)

# Group overlapped regions
df_nuc$group <- cumsum(
  c(1, 
    with(df_nuc[-1,], 
         seqnames[-1] != seqnames[-nrow(df_nuc)] |
           start[-1] > end[-nrow(df_nuc)])
  )
)

######
df_subnuc <- do.call(rbind, lapply(list_subnuc, function(gr) {
  df <- as.data.frame(gr)
  df[, c("seqnames", "start", "end")]
}))

df_subnuc <- df_subnuc %>%
  filter(seqnames %in% standard_chr)

df_subnuc <- df_subnuc %>%
  arrange(seqnames, start, end)

# group overlapped regions manually (avoid reduce crashing)
df_subnuc$group <- cumsum(
  c(1, 
    with(df_subnuc[-1,], 
         seqnames[-1] != seqnames[-nrow(df_subnuc)] |
           start[-1] > end[-nrow(df_subnuc)])
  )
)

consensus_nuc <- df_nuc %>%
  group_by(seqnames, group) %>%
  summarise(
    start = min(start),
    end = max(end),
    n = n(),
    .groups = "drop"
  )

consensus_subnuc <- df_subnuc %>%
  group_by(seqnames, group) %>%
  summarise(
    start = min(start),
    end = max(end),
    n = n(),
    .groups = "drop"
  )

################
ATACnuc.gr <- GRanges(
  seqnames = consensus_nuc$seqnames,
  ranges = IRanges(consensus_nuc$start, consensus_nuc$end),
  score = consensus_nuc$n
)

ATACsubnuc.gr <- GRanges(
  seqnames = consensus_subnuc$seqnames,
  ranges = IRanges(consensus_subnuc$start, consensus_subnuc$end),
  score = consensus_subnuc$n
)

hits.nuc <- findOverlaps(hotspots_gr, ATACnuc.gr)
hits.subnuc <- findOverlaps(hotspots_gr, ATACsubnuc.gr)

ATAC_nuc <- tapply(
  ATACnuc.gr$score[subjectHits(hits.nuc)],
  queryHits(hits.nuc),
  mean
)

ATAC_subnuc <- tapply(
  ATACsubnuc.gr$score[subjectHits(hits.subnuc)],
  queryHits(hits.subnuc),
  mean
)

ATAC_nuc_full <- rep(0, length(hotspots_gr))
ATAC_nuc_full[as.numeric(names(ATAC_nuc))] <- ATAC_nuc

ATAC_subnuc_full <- rep(0, length(hotspots_gr))
ATAC_subnuc_full[as.numeric(names(ATAC_subnuc))] <- ATAC_subnuc

hotspots_gr$chromatin_state <- case_when(
  ATAC_subnuc_full > 0 & ATAC_nuc_full == 0 ~ "open",
  ATAC_subnuc_full > 0 & ATAC_nuc_full > 0 ~ "poised",
  ATAC_subnuc_full == 0 & ATAC_nuc_full > 0 ~ "nucleosomal",
  TRUE ~ "closed"
)

ATAC <- ifelse(hotspots_gr$chromatin_state=="closed",0,1)

################################################################################
################################################################################
#Hotspots overlap with modified histones analysis
################################################################################
################################################################################

#Upload wig files from Jung et al 2017 GSE79227 (this takes a while to process, come tomorrow)

marks <- list(
  H3k9ac   = paths_chip$H3K9ac,
  H3k27ac  = paths_chip$H3K27ac,
  H3k4me1  = paths_chip$H3K4me1,
  H3k4me3  = paths_chip$H3K4me3,
  H3k36me3 = paths_chip$H3K36me3,
  H3k27me3 = paths_chip$H3K27me3,
  H3k9me3  = paths_chip$H3K9me3
)

#create function to process wig files
analyze_histone_mark <- function(hotspots_gr, 
                                 file_rep1, 
                                 file_rep2, 
                                 mark_name,
                                 quantile_thr = 0.75, 
                                 plot = TRUE){
  
  library(GenomicRanges)
  library(rtracklayer)
  
  rep1 <- import(file_rep1)
  rep2 <- import(file_rep2)
  
  seqlevels(rep1) <- sub("\\.chr.*", "", seqlevels(rep1))
  seqlevels(rep2) <- sub("\\.chr.*", "", seqlevels(rep2))
  
  hits1 <- findOverlaps(hotspots_gr, rep1)
  hits2 <- findOverlaps(hotspots_gr, rep2)
  
  signal1 <- tapply(rep1$score[subjectHits(hits1)], queryHits(hits1), mean)
  signal2 <- tapply(rep2$score[subjectHits(hits2)], queryHits(hits2), mean)
  
  #convert to clean vector
  vec1 <- rep(NA_real_, length(hotspots_gr))
  vec2 <- rep(NA_real_, length(hotspots_gr))
  
  if(length(signal1) > 0){
    vec1[as.numeric(names(signal1))] <- as.numeric(signal1)
  }
  if(length(signal2) > 0){
    vec2[as.numeric(names(signal2))] <- as.numeric(signal2)
  }
  
  df <- data.frame(
    rep1 = vec1,
    rep2 = vec2
  )
  
  df$signalWT <- rowMeans(df, na.rm = TRUE)
  
  if(plot){
    pheatmap::pheatmap(df)
  }
  
  all_signal <- c(mcols(rep1)$score, mcols(rep2)$score)
  all_signal_nonzero <- all_signal[all_signal > 0]
  
  thr <- quantile(all_signal_nonzero, probs = quantile_thr, na.rm = TRUE)
  
  df$signal_binary <- ifelse(df$signalWT > thr, 1, 0)
  
  return(df)
}

#apply the function
results <- lapply(names(marks), function(m){
  
  files <- marks[[m]]
  
  analyze_histone_mark(
    hotspots_gr,
    files[1],
    files[2],
    mark_name = m,
    plot = FALSE
  )
})

names(results) <- names(marks)

###############################################################################
###############################################################################
#create data frame with binary signals to plot heatmap

histones <- data.frame(
  ATAC = ATAC,
  H3k9ac = results$H3k9ac$signal_binary,
  H3k27ac = results$H3k27ac$signal_binary,
  H3k4me1 = results$H3k4me1$signal_binary,
  H3k4me3 = results$H3k4me3$signal_binary,
  H3k36me3 = results$H3k36me3$signal_binary,
  H3k27me3 = results$H3k27me3$signal_binary,
  H3k9me3 = results$H3k9me3$signal_binary
)

pheatmap(t(histones), cluster_rows = FALSE, cluster_cols = FALSE)
