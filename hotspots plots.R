################################################################################
# FIGURE 1E
# OVERLAP BETWEEN TF MOTIFS AND HIGH RELATIVE HISTONE SIGNAL
################################################################################

library(ggplot2)

# ------------------------------------------------------------------------------
# 1. Read binary matrices
# ------------------------------------------------------------------------------

histone_matrix <- read.csv2(
  "histone high signal binary.csv"
)

rownames(histone_matrix) <- histone_matrix$X
histone_matrix <- histone_matrix[, -1, drop = FALSE]
histone_matrix <- as.matrix(histone_matrix)


tf_matrix <- read.csv2(
  "TF binding at hotspots matrix.csv",check.names = FALSE)

rownames(tf_matrix) <- tf_matrix[,1]
tf_matrix <- tf_matrix[, -1, drop = FALSE]
tf_matrix <- as.matrix(tf_matrix)


# ------------------------------------------------------------------------------
# 2. Check hotspot identifiers
# ------------------------------------------------------------------------------

histone_ids <- rownames(histone_matrix)
tf_ids <- colnames(tf_matrix)

if (!setequal(histone_ids, tf_ids)) {
  
  missing_in_tf <- setdiff(histone_ids, tf_ids)
  missing_in_histones <- setdiff(tf_ids, histone_ids)
  
  stop(
    paste0(
      "Hotspot identifiers do not match.\n",
      "Histones but not TF: ",
      paste(missing_in_tf, collapse = ", "),
      "\nTF but not histones: ",
      paste(missing_in_histones, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------------------------
# 3. Put matrices in the same hotspot order
# ------------------------------------------------------------------------------

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


# ------------------------------------------------------------------------------
# 4. Define histone-positive and TF-positive hotspots
# ------------------------------------------------------------------------------

histone_evaluable <- rowSums(
  !is.na(histone_matrix)
) > 0

histone_positive <- rowSums(
  histone_matrix,
  na.rm = TRUE
) > 0

tf_positive <- colSums(
  tf_matrix,
  na.rm = TRUE
) > 0

names(histone_evaluable) <- all_hotspots
names(histone_positive) <- all_hotspots
names(tf_positive) <- all_hotspots


# ------------------------------------------------------------------------------
# 5. Assign categories
# ------------------------------------------------------------------------------

category <- rep(
  NA_character_,
  length(all_hotspots)
)

names(category) <- all_hotspots

category[
  !histone_evaluable
] <- "Histones not evaluable"

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

stopifnot(!anyNA(category))


# ------------------------------------------------------------------------------
# 6. Table hotspot by hotspot
# ------------------------------------------------------------------------------

classification <- data.frame(
  hotspot = all_hotspots,
  histone_positive = histone_positive,
  TF_positive = tf_positive,
  category = category,
  stringsAsFactors = FALSE
)

classification$HS_number <- as.numeric(
  sub("HS", "", classification$hotspot)
)

classification <- classification[
  order(classification$HS_number),
]

classification$HS_number <- NULL

print(classification, row.names = FALSE)


# ------------------------------------------------------------------------------
# 7. Sanity checks
# ------------------------------------------------------------------------------

cat(
  "\nHistone-positive:",
  sum(histone_positive),
  "/",
  length(histone_positive),
  "\n"
)

cat(
  "TF-positive:",
  sum(tf_positive),
  "/",
  length(tf_positive),
  "\n\n"
)

print(
  table(classification$category)
)


# ------------------------------------------------------------------------------
# 8. Prepare Figure 1E
# ------------------------------------------------------------------------------

category_levels <- c(
  "Neither",
  "TF only",
  "Histones only",
  "Histones + TF",
  "Histones not evaluable"
)

classification$category <- factor(
  classification$category,
  levels = category_levels
)

counts <- aggregate(
  hotspot ~ category,
  data = classification,
  FUN = function(x) paste(x, collapse = "\n"),
  drop = FALSE
)

counts$n <- vapply(
  counts$hotspot,
  function(x) {
    
    if (is.na(x) || x == "") {
      return(0L)
    }
    
    length(
      strsplit(x, "\n", fixed = TRUE)[[1]]
    )
  },
  integer(1)
)

counts <- counts[
  counts$n > 0,
  ,
  drop = FALSE
]

counts$percent <- 100 * counts$n / length(all_hotspots)

counts$label <- paste0(
  counts$category,
  " (",
  counts$n,
  ")\n(",
  round(counts$percent, 1),
  "%)\n\n",
  counts$hotspot
)

stopifnot(
  sum(counts$n) == length(all_hotspots)
)

print(
  counts[, c(
    "category",
    "n",
    "percent"
  )]
)


# ------------------------------------------------------------------------------
# 9. Plot
# ------------------------------------------------------------------------------

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
    title = "Regulatory features across cocaine hotspots",
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
