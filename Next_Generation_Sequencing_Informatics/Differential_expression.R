setwd("/Users/taniyapal/Downloads")

BiocManager::install("edgeR")
library(edgeR)
#importing the counts table
counts=read.delim("HTSeq_counts.tabular", header=T, stringsAsFactors=F)
rownames(counts)=counts[,1]
counts=counts[,-1]
counts
group=factor(c(1,1,1,2,2,2))
#storing the data in asimple list based data object
y=DGEList(counts=counts, group=group)

#filtering the lowly expressed genes
keep=filterByExpr(y)
y=y[keep,,keep.lib.sizes=F]#lib.sizes refers to library size or sequencing depth for each sample
#Normalizing the library size
y=calcNormFactors(y)
#estimating the common Dispersion
dge=estimateCommonDisp(y)
#estimating the Tagwise Dispersion
dge=estimateTagwiseDisp(dge)
#Testing for differentially expressed genes with Fisher's Test
dgeTest=exactTest(dge)
resNoFilt=topTags(dgeTest, n=nrow(dgeTest$table))
#Filtering out the down regulated genes
sigDownReg=resNoFilt$table[resNoFilt$table$FDR<0.05,]
sigDownReg=sigDownReg[order(sigDownReg$logFC),]
#Filtering out the up regulated genes
sigUpReg=sigDownReg[order(sigDownReg$logFC, decreasing=T),]
#Filtering out the top 100 up regulated differentially expressed genes
Top100UpRegGene=sigUpReg[1:100, ]``
write.csv(sigDownReg, file="Downregulated_genes.csv")
write.csv(sigUpReg, file="UpRegulated_genes.csv")
#Filtering out the Fold change and negative p values
volcanoData=cbind(resNoFilt$table$logFC, -log10(resNoFilt$table$FDR))
colnames(volcanoData)=c("logFC", "negLogPval")
#plotting the volcano plot
plot(volcanoData, pch=19, main="Volcano plot for differentially expressed genes")

#annotation

BiocManager::install("biomaRt")
library(biomaRt)
ensembl=useMart(biomart="plants_mart", dataset="stuberosum_eg_gene", host="plants.ensembl.org")
View(listAttributes(ensembl, page="feature_page"))
View(listFilters(ensembl))
r=rownames(Top100UpRegGene)
s=as.vector(r)
annot.biomart=getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id", "description","chromosome_name", "start_position", "end_position", "transcript_start", "transcript_end", "strand", "kegg_enzyme"), mart=ensembl, values=s)
annot.biomart
annot.table=merge(x=Top100UpRegGene, y=annot.biomart)
annot.table
write.csv(annot.table, file="annot_biomart_potato.csv", row.names=T)
