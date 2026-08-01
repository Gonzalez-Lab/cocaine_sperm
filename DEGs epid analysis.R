#############################################################################
#Epidydimosomes analysis

library(pheatmap)

DEGs <- read.csv2("DEGs RNA sperm.csv")

sharma <- read.csv2("table sharma2016.csv")

epid <- sharma[which(sharma$Name %in% DEGs$symbol),]

#mir6326 appear as mir2 and gene categories
epid <- epid[which(epid$Class == "genes"),]

rownames(epid) <- epid$Name

epid <- epid[,-c(1:3,6:9,15)]

epid <- log2(epid+1)

pheatmap(epid, 
         #color=turbo(100),
         cutree_cols = 4, 
         cutree_rows = 5,
         fontsize = 7.5
)

