#----------------------------------------------------------
# Diverging bar plot: direction and transcript biotype
# Aesthetic matched to the volcano plot
#----------------------------------------------------------

library(ggplot2)
library(dplyr)

# Verify these counts before final use
transcript.summary <- data.frame(
  biotype = c(
    "Coding",
    "Coding",
    "Non-coding",
    "Non-coding"
  ),
  direction = c(
    "DOWN",
    "UP",
    "DOWN",
    "UP"
  ),
  count = c(
    70,
    1,
    3,
    1
  )
)

# Negative values are used to place downregulated transcripts
# on the left side of the plot
transcript.summary <- transcript.summary %>%
  mutate(
    plot_count = ifelse(
      direction == "DOWN",
      -count,
      count
    )
  )

# Order of categories
transcript.summary$biotype <- factor(
  transcript.summary$biotype,
  levels = c("Non-coding", "Coding")
)

transcript.summary$direction <- factor(
  transcript.summary$direction,
  levels = c("DOWN", "UP")
)

# Create plot
transcript.summary.plot <- ggplot(
  transcript.summary,
  aes(
    x = plot_count,
    y = biotype,
    fill = direction
  )
) +
  geom_col(
    width = 0.55,
    color = "black",
    linewidth = 0.35
  ) +
  
  # Number of transcripts
  geom_text(
    aes(
      label = count,
      hjust = ifelse(plot_count < 0, 1.15, -0.35)
    ),
    size = 4.5,
    fontface = "bold"
  ) +
  
  # Central line separating down- and upregulated transcripts
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    color = "black"
  ) +
  
  scale_fill_manual(
    values = c(
      "DOWN" = "dodgerblue2",
      "UP" = "deeppink2"
    )
  ) +
  
  scale_x_continuous(
    breaks = c(-70, -50, -30, -10, 0, 10),
    labels = abs,
    limits = c(-78, 10),
    expand = expansion(mult = c(0, 0))
  ) +
  
  labs(
    title = "Direction and biotype of differentially enriched transcripts",
    x = "Number of transcripts",
    y = NULL,
    fill = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16,
      margin = margin(b = 12)
    ),
    
    legend.position = "top",
    legend.text = element_text(
      size = 11,
      color = "black"
    ),
    legend.key.size = unit(0.7, "cm"),
    
    axis.title.x = element_text(
      face = "bold",
      size = 13,
      margin = margin(t = 8)
    ),
    
    axis.text.x = element_text(
      size = 11,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 12,
      face = "bold",
      color = "black"
    ),
    
    axis.ticks.y = element_blank(),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    panel.grid = element_blank(),
    
    plot.margin = margin(
      t = 10,
      r = 25,
      b = 10,
      l = 10
    )
  )

# Show plot
transcript.summary.plot

