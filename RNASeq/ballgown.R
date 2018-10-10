# GSE78011
source("https://bioconductor.org/biocLite.R")
biocLite(c("ballgown", "genefilter", "devtools", "tidyverse"))
devtools::install_github('alyssafrazee/RSkittleBrewer')
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)

setwd("~/data/ballgown")
# phenotype data for the samples
project <- "GSE75168"
SRX <- c("SRX1438068", "SRX1438069", "SRX1438070", "SRX1438074", "SRX1438075", "SRX1438076")
cell_type <- c(rep("MCF10A", 3), rep("MDA-MB-231", 3))
pheno_data <- data.frame(SRX, cell_type)

# read expression data calculated by StringTie
bg <- ballgown(dataDir = project, samplePattern = "SRX", pData=pheno_data)

# filter out low abundance genes
bg_filt <- subset(bg,"rowVars(texpr(bg)) >1",genomesubset=TRUE)

# identify differentially expressed transcripts and genes
results_transcripts <- stattest(bg_filt, feature="transcript",covariate="cell_type", getFC=TRUE, meas="FPKM")
results_transcripts <- data.frame(geneNames=ballgown::geneNames(bg_filt), geneIDs=ballgown::geneIDs(bg_filt), results_transcripts)

results_genes <- stattest(bg_filt, feature="gene",covariate="cell_type", getFC=TRUE, meas="FPKM")
indices <- match(results_genes$id, texpr(bg, 'all')$gene_id)
gene_names_for_result <- texpr(bg, 'all')$gene_name[indices]
results_genes <- data.frame(geneNames=gene_names_for_result, results_genes)

results_transcripts <- arrange(results_transcripts,pval)
results_genes <- arrange(results_genes,pval)

write.csv(results_transcripts, paste0(project, "/TNBC_transcript.csv"), row.names=FALSE)
write.csv(results_genes, paste0(project, "/TNBC_gene.csv"), row.names=FALSE)

results_transcripts <- subset(results_transcripts,results_transcripts$qval<0.05)
results_genes <- subset(results_genes,results_genes$qval<0.05)

write.csv(results_transcripts, paste0(project, "/TNBC_transcript_q_0.05.csv"), row.names=FALSE)
write.csv(results_genes, paste0(project, "/TNBC_gene_q_0.05.csv"), row.names=FALSE)
