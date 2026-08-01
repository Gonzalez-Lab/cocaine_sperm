###############################################################################
###############################################################################
# hotspots length barplot
###############################################################################
###############################################################################
library(ggplot2)

hotspots <- read.csv2("hotspots.csv")

# labels: chr_HS_gene
hotspots$label <- paste0(
  "crh",
  hotspots$chr, "_",
  hotspots$hotspot, "_",
  hotspots$gene_name
)

# Ordenar por longitud
hotspots$label <- factor(
  hotspots$label,
  levels = hotspots$label[order(hotspots$length_bp)]
)

ggplot(
  hotspots,
  aes(
    x = length_bp,
    y = label,
    fill = dir.Methy
  )
) +
  geom_col(width = 0.75) +
  geom_text(
    aes(label = length_bp),
    hjust = -0.15,
    size = 3.5
  ) +
  scale_fill_manual(
    values = c(
      hyper = "deeppink3",
      hypo = "dodgerblue4"
    ),
    name = "Direction"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Length of cocaine-associated methylation hotspots",
    x = "Length (bp)",
    y = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.text.y = element_text(
      color = "black",
      size = 9
    )
  )


################################################################################
# OVERLAP BETWEEN TF MOTIFS AND HIGH HISTONE SIGNAL
################################################################################
#to run this you need to run first the TF motif enrichment and histone analysis scripts

histone_matrix <- read.csv2("histone high signal binary.csv")
rownames(histone_matrix) <- histone_matrix$X
histone_matrix <- histone_matrix[,-1]

tf_matrix <- read.csv2("TF binding at hotspots matrix.csv")
rownames(tf_matrix) <- tf_matrix$X
tf_matrix <- tf_matrix[,-1]

# Convert explicitly to matrices
histone_matrix <- as.matrix(histone_matrix)
tf_matrix <- as.matrix(heatmap_matrix_tf)

# --------------------------------------------------------------------------
# 1. Validate hotspot identifiers
# --------------------------------------------------------------------------

histone_ids <- rownames(histone_matrix)
tf_ids <- colnames(tf_matrix)

if (is.null(histone_ids) || is.null(tf_ids)) {
  stop("One of the matrices does not contain hotspot names.")
}

if (anyDuplicated(histone_ids)) {
  stop("Duplicated hotspot names were found in histones_high.")
}

if (anyDuplicated(tf_ids)) {
  stop("Duplicated hotspot names were found in heatmap_matrix_tf.")
}

# Do not use union(): both matrices must contain the same hotspot universe
if (!setequal(histone_ids, tf_ids)) {
  
  missing_in_tf <- setdiff(histone_ids, tf_ids)
  missing_in_histones <- setdiff(tf_ids, histone_ids)
  
  stop(
    paste0(
      "The hotspot identifiers in both matrices do not match.\n",
      "Present in histones but absent from TF matrix: ",
      paste(missing_in_tf, collapse = ", "),
      "\nPresent in TF matrix but absent from histone matrix: ",
      paste(missing_in_histones, collapse = ", ")
    )
  )
}

# --------------------------------------------------------------------------
# 2. Put both matrices in exactly the same hotspot order
# --------------------------------------------------------------------------

all_hotspots <- histone_ids

tf_matrix <- tf_matrix[
  ,
  all_hotspots,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(histone_matrix),
    colnames(tf_matrix)
  )
)

cat(
  "Common hotspot universe:",
  length(all_hotspots),
  "\n"
)

# --------------------------------------------------------------------------
# 3. Validate matrix values
# --------------------------------------------------------------------------

if (!all(histone_matrix %in% c(0, 1) | is.na(histone_matrix))) {
  stop("histones_high contains values other than 0, 1 or NA.")
}

if (!all(tf_matrix %in% c(0, 1) | is.na(tf_matrix))) {
  stop("heatmap_matrix_tf contains values other than 0, 1 or NA.")
}

# --------------------------------------------------------------------------
# 4. Determine which hotspots are evaluable for histones
# --------------------------------------------------------------------------

# A hotspot is non-evaluable only if every histone value is NA
histone_evaluable <- rowSums(
  !is.na(histone_matrix)
) > 0

# Histone-positive means at least one high-signal histone mark
histone_positive <- rowSums(
  histone_matrix,
  na.rm = TRUE
) > 0

# TF-positive means at least one selected TF motif
tf_positive <- colSums(
  tf_matrix,
  na.rm = TRUE
) > 0

# Check correspondence after alignment
stopifnot(
  length(histone_positive) == length(all_hotspots),
  length(tf_positive) == length(all_hotspots)
)

names(histone_evaluable) <- all_hotspots
names(histone_positive) <- all_hotspots
names(tf_positive) <- all_hotspots

# --------------------------------------------------------------------------
# 5. Assign one mutually exclusive category to every hotspot
# --------------------------------------------------------------------------

category <- rep(
  NA_character_,
  length(all_hotspots)
)

names(category) <- all_hotspots

category[!histone_evaluable] <-
  "Histones not evaluable"

category[
  histone_evaluable &
    !histone_positive &
    !tf_positive
] <- "Neither"

category[
  histone_evaluable &
    !histone_positive &
    tf_positive
] <- "TF only"

category[
  histone_evaluable &
    histone_positive &
    !tf_positive
] <- "Histones only"

category[
  histone_evaluable &
    histone_positive &
    tf_positive
] <- "Histones + TF"

if (anyNA(category)) {
  stop("At least one hotspot could not be assigned to a category.")
}

# --------------------------------------------------------------------------
# 6. Create plotting table
# --------------------------------------------------------------------------

category_levels <- c(
  "Neither",
  "TF only",
  "Histones only",
  "Histones + TF",
  "Histones not evaluable"
)

counts <- data.frame(
  hotspot = all_hotspots,
  category = factor(
    category,
    levels = category_levels
  ),
  stringsAsFactors = FALSE
)

counts <- aggregate(
  hotspot ~ category,
  data = counts,
  FUN = function(x) paste(x, collapse = "\n"),
  drop = FALSE
)

counts$n <- vapply(
  counts$hotspot,
  function(x) {
    
    if (is.na(x) || x == "") {
      return(0L)
    }
    
    length(strsplit(x, "\n", fixed = TRUE)[[1]])
  },
  integer(1)
)

# Remove empty categories from the plot
counts <- counts[
  counts$n > 0,
  ,
  drop = FALSE
]

counts$label <- paste0(
  counts$category,
  " (",
  counts$n,
  ")\n\n",
  counts$hotspot
)

# Final validation: every hotspot must be counted exactly once
stopifnot(
  sum(counts$n) == length(all_hotspots)
)

print(
  counts[, c("category", "n")]
)

cat(
  "Total hotspots represented:",
  sum(counts$n),
  "\n"
)

# --------------------------------------------------------------------------
# 7. Plot
# --------------------------------------------------------------------------

ggplot(
  counts,
  aes(
    x = "",
    y = n,
    fill = category
  )
) +
  geom_col(
    width = 0.7,
    color = "white",
    position = position_stack(reverse = TRUE)
  ) +
  geom_text(
    aes(label = label),
    position = position_stack(
      vjust = 0.5,
      reverse = TRUE
    ),
    size = 4,
    lineheight = 0.9
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Neither" = "grey75",
      "TF only" = "gold",
      "Histones only" = "cadetblue1",
      "Histones + TF" = "lightpink",
      "Histones not evaluable" = "grey40"
    ),
    drop = FALSE
  ) +
  labs(
    title = paste0(
      "Regulatory features across cocaine-associated ",
      "methylation hotspots"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )
