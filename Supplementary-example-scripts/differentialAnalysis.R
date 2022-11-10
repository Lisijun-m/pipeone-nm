#!/usr/bin/env Rscript
## use grass carp RNA-seq as test data, the gene_count_matrix is circRNA host gene count matrix
## grass carp GEO: GSE185170

library("DESeq2","ggplot2","pheatmap","RColorBrewer")

countData<-read.csv("gene_count_matrix.csv")
rownames(countData)<-countData[,1]
countData<-countData[,-1]

condition <- factor(c("healthy","6persalineInjured","6persalineInjured",
"6persalineInjured","6per","6per","6per","3persalineInjured" ,"3persalineInjured","3persalineInjured","3per","3per","3per","healthy","healthy"),levels=c("healthy","6persalineInjured","6per","3per","3persalineInjured"))
colData<-data.frame(row.names=colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData,design = ~ condition )

dds <- dds [ rowSums(counts(dds)) >10, ]
dds <- DESeq(dds)
rld <- rlog(dds, blind = FALSE)
#Plot PCA Plot
plotPCA(rld,intgroup=c("condition"))
#Plot Correlation Heatmap
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors)\

results(dds)
res <- results(dds,contrast = c("condition","3persalineInjured","3per"))
diff<-res
diff <- na.omit(diff)
foldChange = 2
padj = 0.05
#differential expressed genes
diffsig <- diff[(diff$pvalue < padj & abs(diff$log2FoldChange) > foldChange),]
write.csv(diffsig, "All_diffsig.csv")




