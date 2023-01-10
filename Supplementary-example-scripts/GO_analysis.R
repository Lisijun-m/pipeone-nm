## GO analysis script for differential expressed genes
## go_forStats.txt can be found in output directory 
## All_diffsig.csv is created in differential analysis part

library(”ggplot2”,”clusterProfiler”)

go_1<-read.table("go_forStats.txt",header=FALSE)
go_2<-read.csv("go_anno_3.csv",header=F)
go_term2gene=data.frame(go_1$V2,go_1$V1)
go_term2name=data.frame(go_2$V1,go_2$V3)
names(go_term2gene)=c("go_term","gene")
names(go_term2name)=c("go_term","name")

diff<-read.csv("All_diffsig.csv",header=T)
geneName<-diff$X
go_enrich <- enricher(gene=geneName,pvalueCutoff = 0.05,pAdjustMethod = "BH",TERM2GENE = go_term2gene,TERM2NAME = go_term2name)
diff_go=as.data.frame(go_enrich)

pdf("results.pdf")
barplot(go_enrich,showCategory=30,font.size=10)
dev.off()
