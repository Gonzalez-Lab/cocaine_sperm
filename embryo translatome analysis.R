################################################################################
#analysis for active translation in embryo
library(dplyr)
library(ggwordcloud)

DEGs <- read.csv2("DEGs RNA sperm.csv")

embryo.translatome = read.csv2("translatome embryo data.csv")

embryo.translatome <- filter(embryo.translatome,!is.nan(TE)) %>%  filter(TE > 0)

DEGs.in.embryo = DEGs[which(DEGs$symbol %in% embryo.translatome$symbol),]

#wordcloud plot

traslatome.DEGs <- embryo.translatome[which(embryo.translatome$symbol %in% DEGs.in.embryo$symbol),]

cluster1 <- c("Atg9a","Ulk2","Cib1","Stk35","Tbc1d5","Epn1")
cluster2 <- c("Tpgs2","Lrrc8b")

traslatome.DEGs$cluster1 <- ifelse(traslatome.DEGs$symbol %in% cluster1, "1",
                                   ifelse(traslatome.DEGs$symbol %in% cluster2,"2","0"))

ggplot(data = traslatome.DEGs,
       aes(label = symbol, size = TE, col = cluster1)) +
  geom_text_wordcloud(rm_outside = TRUE, max_steps = 2,
                      grid_size = 5, eccentricity =0.9)+
  theme_void()+
  scale_color_manual(values = c("black", "#d52954","#515c9a"))


