setwd("/lustre1/yz73026/Bisulfite/")
library(TCGAbiolinks)
library(foreach)
library(doMC)
registerDoMC(9)
query <- list()
a <- c("TCGA-BLCA","TCGA-GBM","TCGA-BRCA","TCGA-LUAD","TCGA-UCEC","TCGA-LUSC","TCGA-STAD","TCGA-COAD","TCGA-READ")
foreach(i=1:length(a)) %dopar%
{
	query[[i]] <- GDCquery(project = a[i],data.category = "DNA methylation",data.type = "Methylation percentage",legacy = T)
	GDCdownload(query[[i]])
}