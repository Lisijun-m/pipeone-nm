## This script is for circRNA differential analysis without biological replications
## gene_count_matrix.csv can be found in directory ciri/ciriquant/diff/project

library(edgeR)

rawdata <- read.csv("gene_count_matrix.csv",header=T,sep=',',row.names='gene_id')

## ## replace condition in the next line according to sample treatment
group <- factor(c('condition1','condition2','condition3'))

y <- DGEList(counts = rawdata, genes = rownames(rawdata), group = group)

keep <- rowSums(cpm(y)>1) >= 1
y <- y[keep,,keep.lib.sizes=FALSE]


y <- calcNormFactors(y)
bcv <- 0.2

## replace control/treatment in the next line
et <- exactTest(y, dispersion=bcv^2,pair=c('control','treatment'))

topTags(et)
summary(de <- decideTestsDGE(et))
detags <- rownames(y)[as.logical(de)]

plotSmear(et, de.tags=detags)
abline(h=c(-2, 2), col="blue")
result = topTags(et, n = nrow(et$table))$table
foldChange = 2
padj = 0.05
diffsig <- result[(result$PValue < padj & abs(result$logFC) > foldChange),]
write.csv(diffsig,'diff_treatment_control.csv',row.names=F)

library(clusterProfiler)
go_1<-read.table("go_forStats.txt",header=FALSE)
go_2<-read.csv("go_anno_3.csv",header=F)

go_term2gene=data.frame(go_1$V2,go_1$V1)
go_term2name=data.frame(go_2$V1,go_2$V3)
names(go_term2gene)=c("go_term","gene")
names(go_term2name)=c("go_term","name")

geneName<-diffsig$genes
go_enrich <- enricher(gene=geneName,pvalueCutoff = 0.05,pAdjustMethod = "BH",TERM2GENE = go_term2gene,TERM2NAME = go_term2name)
## Second Plot
barplot(go_enrich,showCategory=15,font.size=15)

head(as.data.frame(go_enrich))
dim(go_enrich)
diff_go=as.data.frame(go_enrich)
write.csv(diff_go,"GO.csv",row.names = F)
