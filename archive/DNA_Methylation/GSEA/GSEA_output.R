# Code to use GSEA_with_p_value.R results and filter genes (p <= 0.05),
# and then perform GSEA, getting pathways' ES (and p-values)

## model.normal (length 4469): 0 (numeric) and lm (list)
## R2.p.tumor (data frame nrow 4170): Gene_Symbol   R_Square    p_value

setwd("/lustre1/yz73026/array")
## Source functions and load data
source("./gsea.R")  # GSEA.EnrichmentScore2, calc.gsea.1
load("./pathways.list.GOMaco.CORUM.RCT.HGNC.Msigdb.gsc.RData")
library(piano)
library(foreach)
library(doMC)
registerDoMC(24)  # 12 projects, normal & tumor
pathways <- lapply(gsc[[1]],
                   function(x) sub("(.*)\\|.*", "\\1", x))
pathways <- pathways[unique(names(pathways))]
weighted.score.type <- 1  # Step length consider R^2 value
valley <- r2.value <- 0.5  # R^2 cutoff
coef.num <- 2

files.p <- list.files("./GSEA161123/",
                      pattern = paste0("coef", coef.num))
projects <- sub("(\\w*)\\.nor.*",  
                "\\1", files.p)  # BLCA.normalGLM_coef2_lm.RData

## Separate normal and tumor files
dir.create("./GSEA161123/sep/", showWarnings = F)
for (i in seq_along(projects)) {
    load(paste0("./GSEA161123/", files.p[i]))
    save(R2.p.normal, model.normal,
         file = paste0("./GSEA161123/sep/", projects[i],
                       "_normal_coef", coef.num, "_lm.RData"))
    save(R2.p.tumor, model.tumor,
         file = paste0("./GSEA161123/sep/", projects[i],
                       "_tumor_coef", coef.num, "_lm.RData"))
}
files.p <- list.files("./GSEA161123/sep/",
                      pattern = paste0("coef", coef.num))
file.names <- sub("(.*)_lm\\.RData", "\\1", files.p)

## Get gene list
load(paste0("./mergeData/", projects[1], ".RData"))
gene.list <- as.character(c(rownames(datan), rownames(datat)))
gene.list <- gene.list[unique(gene.list)]

## Run GSEA
foreach (i = seq_along(files.p)) %dopar% {
    load(paste0("./GSEA161123/sep/", files.p[i]))
    
}