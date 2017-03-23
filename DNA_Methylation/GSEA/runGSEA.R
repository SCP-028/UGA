setwd("/lustre1/yz73026/array/")
library(foreach)
library(doMC)
registerDoMC(12)

source("./GSEA_func.R")
load("./pathways.list.GOMaco.CORUM.RCT.HGNC.Msigdb.gsc.RData")
# Get the pathways' names right
pathways <- lapply(gsc[[1]],
                   function(x) sub("(.*)\\|.*", "\\1", x))
pathways <- pathways[!duplicated(names(pathways))]
rm(gsc)

tag.R2 <- "./glmnet_TF/results/moreTF/summary/"
files.R2 <- list.files(tag.R2, pattern = "\\.RData")
projects <- sub("(.*)\\.RData", "\\1", files.R2)
weighted.score.type <- 1  # Step length consider R^2 value
valley.r2.value <- 0
gene.list <- getGenelist("./mergeData/")

foreach (i = seq_along(projects)) %dopar% {
    load(paste0(tag.R2, projects[i], ".RData"))
    for (j in 2:5){
        temp <- extractR2(results, j, valley.r2.value)
        gene.set <- as.character(temp$gene_symbol)
        t.stat <- as.numeric(temp[ ,"r_squared"])
        #t.stat <- t.stat - mean(t.stat)
        names(t.stat) <- gene.set
        temp <- calc.gsea.1(t.stat, weighted.score.type)
        write.table(temp, file = paste0("./GSEA_TF/", projects[i],
                                        ".coef", j, ".txt"))
        print(paste0("Finished GSEA for project: ", projects[i],
                     "; Group: coef", j))
    }
}
