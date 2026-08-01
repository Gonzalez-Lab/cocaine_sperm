#----------------------------------------------------------
# Volcano plot with EnhancedVolcano
#----------------------------------------------------------

library(EnhancedVolcano)
library(dplyr)

# Convert DESeq2 results to data frame
volcano.data <- as.data.frame(res)

# Use only gene symbols as labels
volcano.data$label <- ifelse(
  is.na(volcano.data$symbol) | volcano.data$symbol == "",
  "",
  volcano.data$symbol
)

# Genes to label (only significant DEGs)
select.genes <- volcano.data %>%
  filter(
    !is.na(symbol),
    symbol != "",
    !is.na(padj),
    padj < 0.1,
    abs(log2FoldChange) > 1
  ) %>%
  pull(symbol)

# Volcano plot
volcano.plot <- EnhancedVolcano(
  volcano.data,
  
  lab = volcano.data$label,
  selectLab = select.genes,
  
  x = "log2FoldChange",
  y = "padj",
  
  pCutoff = 0.1,
  FCcutoff = 1,
  
  xlim = c(-5, 5),
  ylim = c(0, 3),
  
  title = "Differentially expressed genes in sperm after cocaine",
  subtitle = NULL,
  caption = NULL,
  
  xlab = expression(Log[2]~fold~change),
  ylab = expression(-Log[10]~adjusted~italic(P)),
  
  pointSize = 3,
  labSize = 4,
  
  col = c(
    "grey80",
    "grey80",
    "grey80",
    "#C44E52"
  ),
  
  colAlpha = 0.9,
  
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = "grey40",
  
  boxedLabels = FALSE,
  
  cutoffLineType = "dashed",
  cutoffLineWidth = 0.6,
  cutoffLineCol = "black",
  
  legendPosition = "top",
  legendLabSize = 12,
  legendIconSize = 4,
  
  legendLabels = c(
    "Not significant",
    "Log2 FC",
    "Adjusted P",
    "padj < 0.1 & |Log2FC| > 1"
  ),
  
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  
  border = "full"
)

# Show plot
volcano.plot

###########################################################################
#----------------------------------------------------------
# PCA plot from rlog-transformed RNA-seq data
#----------------------------------------------------------

library(ggplot2)
library(ggrepel)

# Extract PCA coordinates from the rlog-transformed object
pca.data <- plotPCA(
  rld,
  intgroup = "group",
  returnData = TRUE
)

# Percentage of variance explained by each component
percent.var <- round(
  100 * attr(pca.data, "percentVar"),
  digits = 1
)

# Improve group labels
pca.data$group <- factor(
  pca.data$group,
  levels = c("veh", "coc"),
  labels = c("Vehicle", "Cocaine")
)

# Add sample names
pca.data$sample <- rownames(pca.data)

# Create PCA plot
pca.plot <- ggplot(
  pca.data,
  aes(
    x = PC1,
    y = PC2,
    color = group,
    shape = group
  )
) +
  geom_point(
    size = 5,
    alpha = 0.9
  ) +
  geom_text_repel(
    aes(label = sample),
    size = 4,
    color = "black",
    box.padding = 0.5,
    point.padding = 0.4,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Vehicle" = "#4F81BD",
      "Cocaine" = "#C44E52"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Vehicle" = 16,
      "Cocaine" = 17
    )
  ) +
  labs(
    title = "Principal component analysis of sperm RNA-seq samples",
    x = paste0("PC1 (", percent.var[1], "%)"),
    y = paste0("PC2 (", percent.var[2], "%)"),
    color = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    legend.position = "top",
    legend.text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    )
  )

# Show plot
pca.plot
