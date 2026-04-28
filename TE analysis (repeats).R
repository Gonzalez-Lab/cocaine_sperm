library(GenomicRanges)
library(data.table)

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
# Transposable Elements TE analysis (repeats)
################################################################################

#download repeats 
repeats <- fread(
  "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz",
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

library(tidyr)
library(dplyr)
library(ggplot2)

#plot percentages of TE class in hybrid hotspots

df_plot_TEhybrid <- hotspots %>%
  separate_rows(TE.class, sep = "-") %>%
  filter(TE.class %in% c("LINE","SINE","LTR")) %>% 
  count(TE.class) %>%
  mutate(prop = n / sum(n),
         label = paste0(TE.class, "\n", round(prop*100, 1), "%"))

ggplot(df_plot_TEhybrid, aes(x = 2, y = n, fill = TE.class)) +
  geom_col(color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 4) +
  theme_void() +
  labs(fill = "TE type in hotspot") +
  scale_fill_manual(values = c(
    "#56B4E9", "#009E73", "#CC79A7"
  ))


#Plot percentages of TE class near hotspots (<1kb) 

df_plot_nearTE <- hotspots %>%
  separate_rows(near_TE, sep = "-") %>%
  count(near_TE) %>%
  mutate(prop = n / sum(n),
         label = paste0(near_TE, "\n", round(prop*100, 1), "%"))

ggplot(df_plot_nearTE, aes(x = 2, y = n, fill = near_TE)) +
  geom_col(color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 4) +
  theme_void() +
  labs(fill = "TE type near hotspot (< 1 kb)") +
  scale_fill_manual(values = c(
    "#56B4E9", "#009E73", "#CC79A7"
  ))


