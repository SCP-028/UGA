# GSE78011
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

# phenotype data for the samples
SRX <- c("SRX1589783" "SRX1589784" "SRX1589785" "SRX1589789" "SRX1589790")
cell_type <- c(rep("MCF-7", 3), rep("MDA-MB-231", 2))
pheno_data <- data.frame(SRX, cell_type)

# read expression data calculated by StringTie
bg_chrX <- ballgown(dataDir = "ballgown", samplePattern = "SRX", pData=pheno_data)

# filter out low abundance genes
bg_chrX_filt <- subset(bg_chrX,"rowVars(texpr(bg_chrX)) >1",genomesubset=TRUE)

# identify differentially expressed transcripts and genes
results_transcripts <- stattest(bg_chrX_filt, feature="transcript",covariate="cell_type", getFC=TRUE, meas="FPKM")
results_transcripts =data.frame(geneNames=ballgown::geneNames(bg_chrX_filt), geneIDs=ballgown::geneIDs(bg_chrX_filt), results_transcripts)
results_genes <- stattest(bg_chrX_filt, feature="gene",covariate="cell_type", getFC=TRUE, meas="FPKM")

results_transcripts <- arrange(results_transcripts,pval)
results_genes <- arrange(results_genes,pval)

write.csv(results_transcripts, "TNBC_transcript.csv", row.names=FALSE)
write.csv(results_genes, "TNBC_gene.csv", row.names=FALSE)

subset(results_transcripts,results_transcripts$qval<0.05)
subset(results_genes,results_genes$qval<0.05)

write.csv(results_transcripts, "TNBC_transcript_q_0.05.csv", row.names=FALSE)
write.csv(results_genes, "TNBC_gene_q_0.05.csv", row.names=FALSE)
