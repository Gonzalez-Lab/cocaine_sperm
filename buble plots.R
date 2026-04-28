library(ggplot2)
library(ggrepel)

func.enrich.meth <- read.csv2("functional enrichment methyl.csv")

func.enrich.meth$Term <- factor(func.enrich.meth$Term, 
                           levels = rev(func.enrich$Term)) 

ggplot(func.enrich.meth, aes(x = as.numeric(Odds.Ratio), y = Term)) +
  # bubles: size by gene count, color by p-adj
  geom_point(aes(size = Ngenes, color = Adjusted.p.value)) +
  
  # Add labels
  geom_text_repel(
    aes(label = genes), 
    size = 3,                  
    box.padding = 0.5,         
    max.overlaps = Inf,        
    segment.color = "grey50"  
  ) +
  
  # labels and title
  labs(
    x = "Odds Ratio",
    y = "GO moleculr function",
    title = "Enrichment Analysis",
    color = "p-adj",
    size = "Gene Count"
  ) +
  
  # visual style
  theme_minimal() +
  theme(
    axis.text.y 
    = element_text(size = 10),
    panel.grid.major = element_line(linetype = "dotted", color = "grey80")
  )

################################################################################
#lncRNAseq data
################################################################################

func.enrich.cluster1 <- read.csv2("func enrich cluster1.csv")

func.enrich.cluster1$Term <- factor(func.enrich.cluster1$Term,
                           levels = rev(func.enrich.cluster1$Term))

ggplot(func.enrich.cluster1, aes(x = as.numeric(strength), y = Term)) +
  # Las burbujas: tamaño por conteo, color por p-valor
  geom_point(aes(size = Ngene, color = p.adj)) +
  
  # Escala de colores (de rojo a azul o similar)
  #  scale_color_gradient(low = "red", high = "blue") +
  # Añadir etiquetas (usamos ggrepel para que se muevan si hay choque)
  geom_text_repel(
    aes(label = genes),
    size = 3,                  # Tamaño de letra pequeño
    box.padding = 0.5,         # Espacio alrededor del punto
    max.overlaps = Inf,        # Forzar a que muestre todos (si no son demasiados)
    segment.color = "grey50"   # Línea que une el texto con el punto
  ) +
  
  # Etiquetas y títulos
  labs(
    x = "Strength",
    y = "",
    title = "Enrichment Analysis cluster 1",
    color = "p-adj",
    size = "Gene Count"
  ) +
  
  # Estilo visual limpio
  theme_minimal() +
  theme(
    axis.text.y
    = element_text(size = 10),
    panel.grid.major = element_line(linetype = "dotted", color = "grey80")
  )

