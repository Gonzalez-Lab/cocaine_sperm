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

#export results
write.csv2(DEGs, "DEGs RNA sperm.csv")

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
         fontsize = 7.5, 
         #color=turbo(100),
         cutree_rows = 2)
