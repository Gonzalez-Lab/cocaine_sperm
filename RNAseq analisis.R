library(dplyr)

#RNAseq analysis cocaine vs vehicle sperm samples
data = read.delim('counts_matrix.txt', header=T)

#matrix preparation
count.matrix = data[,c(7:10)]
colnames(count.matrix) = c("veh1","veh2","coc1","coc2")
rownames(count.matrix) = data$Geneid
count.matrix = na.omit(count.matrix)
count.matrix <- as.matrix(count.matrix)

#Statistical analysis
library(DESeq2)

#sample table preparation
sample.table <- as.data.frame(matrix(c("veh","veh","coc","coc"), ncol = 1))
rownames(sample.table) <- colnames(count.matrix)
colnames(sample.table) <- c("group")
sample.table$group <- factor(sample.table$group, levels = c("veh","coc"))

#Check tables are OK
all(rownames(sample.table) == colnames(count.matrix))

#We make a *DESeqDataSet* (dds) from a count matrix and column data
dds <- DESeqDataSetFromMatrix(countData=count.matrix, 
                              colData=sample.table, 
                              design=~group)

#filter low count genes
keep <- rowSums(counts(dds)) >= 30
dds <- dds[keep, ]
print(paste("Genes retained after filtering:", nrow(dds)))

# add analysis to dds object
dds <- DESeq(dds)

#results
res <- results(dds)

summary(res)

table(res$padj < 0.1)
table(res$padj < 0.1 & res$log2FoldChange< -1) # downregulated
table(res$padj < 0.1 & res$log2FoldChange> 1) # upregulated

# edit ensembl IDs for mapID
genes <- rownames(res)
genes <- substr(genes,1,18)

#convert ensembl to symbol and add it to results
library(org.Mm.eg.db)

res$symbol <- mapIds(org.Mm.eg.db, keys = genes, keytype = "ENSEMBL", column="SYMBOL")

head(res)

UP = res[which(res$padj < 0.1 & res$log2FoldChange > 1),] %>% na.omit()
DOWN = res[which(res$padj < 0.1 & res$log2FoldChange< -1),] %>% na.omit()

DEGs = rbind(UP,DOWN)

#search genes of interest
res[which(res$symbol=="Sgms1"),]

#heatmap de DEGs
library(pheatmap)
library(viridis)

indexDEGs <- which(res$symbol %in% DEGs$symbol)

#covert rlog
rld <- rlog(dds)

#A PCA plot and a heatmap of the top genes: we need to use the rld object
plotPCA(rld, intgroup="group")

#heatmap for DEGs
mat <- assay(rld)[indexDEGs,]
mat <- mat - rowMeans(mat) # hago el z score al restar cada valor menos la media

mat.symbol <- mapIds(org.Mm.eg.db, keys = substr(rownames(mat),1,18), 
                     keytype = "ENSEMBL", column="SYMBOL")

rownames(mat) <- mat.symbol

head(mat)

pheatmap(mat, 
         fontsize = 6, 
         color=magma(100),
         cutree_rows = 2
)



################################################################################
#analysis for active translation in embryo
embryo.translatome = read.csv2("translatome embryo data.csv")

embryo.translatome <- filter(embryo.translatome, TE>0) %>% filter(TE!= NaN)

DEGs.in.embryo = DEGs[which(DEGs$symbol %in% embryo.translatome$symbol),]

#wordcloud plot

traslatome.DEGs <- embryo.translatome[which(embryo.translatome$symbol %in% DEGs.in.embryo$symbol),]

cluster1 <- c("Atg9a","Ulk2","Cib1","Stk35","Tbc1d5","Epn1")
cluster2 <- c("Tpgs2","Lrrc8b")

traslatome.DEGs$cluster1 <- ifelse(traslatome.DEGs$symbol %in% cluster1, "1",
                                   ifelse(traslatome.DEGs$symbol %in% cluster2,"2","0"))

library(ggwordcloud)
ggplot(data = traslatome.DEGs,
       aes(label = symbol, size = TE, col = cluster1)) +
  geom_text_wordcloud(rm_outside = TRUE, max_steps = 2,
                      grid_size = 5, eccentricity =0.9)+
  theme_void()+
  scale_color_manual(values = c("black", "#d52954","#515c9a"))


#############################################################################
#Epidydimosomes analysis

sharma <- read.csv2("table sharma2016.csv")

epid <- sharma[which(sharma$Name %in% DEGs$symbol),]

#mir6326 appear as mir2 and gene categories
epid <- epid[which(epid$Class == "genes"),]

rownames(epid) <- epid$Name

epid <- epid[,-c(1:3,7:9,15)]

epid <- log2(epid+1)

pheatmap(epid, 
         color=magma(100),
         cutree_cols = 4, 
         cutree_rows = 5,
         fontsize = 6
)

