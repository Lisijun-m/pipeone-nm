## Differential analysis of circRNA, the gene_count_matrix is circRNA host gene count matrix
## gene_count_matrix can be found in output directory ciri/ciriquant/diff

## replace condition in the next line according to sample treatment
condition <- factor(c("condition1","condition1","condition1","condition2","condition2","condition2","condition3","condition3" ,"condition3"),levels=c("condition1","condition2","condition3"))


library("DESeq2","ggplot2","pheatmap","RColorBrewer")

countData<-read.csv("gene_count_matrix.csv")
rownames(countData)<-countData[,1]
countData<-countData[,-1]

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

## make comparision between different samples
## replace the control/treatment in the next line
res <- results(dds,contrast = c("condition","treatment","controlr"))

diff<-res
diff <- na.omit(diff)
foldChange = 2
padj = 0.05
#differential expressed genes
diffsig <- diff[(diff$pvalue < padj & abs(diff$log2FoldChange) > foldChange),]
write.csv(diffsig, "All_diffsig.csv")




