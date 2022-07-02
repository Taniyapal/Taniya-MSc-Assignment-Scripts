setwd("/Users/taniyapal/Downloads")
library(edgeR)
counts=read.delim("HTSeq_counts.tabular", header=T, stringsAsFactors=F, col.names=c("Gene_ID","w0_R1", "w0_R2", "w0_R3", "w1_R1", "w1_R2", "w1_R3" ))
counts
rownames(counts)=counts[,1]
counts=counts[,-1]
counts
group=factor(c(1,1,1,2,2,2))
y=DGEList(counts=counts, group=group)
PCAMODEL=prcomp(y, retx=T, center=T, scale=F)
library(factoextra)
fviz_pca_var(PCAMODEL, col.var="blue",geom="text")
