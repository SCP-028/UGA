# Uses glmnet results and perform a Gene Set Enrichment Analysis
# with Enrichment Scores and p-values calculated.

## Resample using permutation test (shuffle order of gene expression)

## Run GSEA with the original model and the new 1000 models, 
## getting a p-value and a Enrichment Score

## Path and libraries
try(setwd("/lustre1/yz73026/array"))
rm(list = ls())
library(foreach)
library(doMC)
registerDoMC(12)  # 12 projects, 4 coefs

## Source codes and load data
source("./get_R2_of_glmnet.R")  # Code for messing with glmnet results
merge.data <- list.files("./mergeData/")
merge.data <- merge.data[order(merge.data)]
glmnet.normal <- list.files("./glmnet/", pattern = "\\.normal")
glmnet.normal <- glmnet.normal[order(glmnet.normal)]
glmnet.tumor <- list.files("./glmnet/", pattern = "\\.tumor")
glmnet.tumor <- glmnet.tumor[order(glmnet.tumor)]
projects <- sub("(.*)\\.RData", "\\1", glmnet.normal)
merge.data <- paste0("./mergeData", merge.data)
glmnet.normal <- paste0("./glmnet", glmnet.normal)
glmnet.tumor <- paste0("./glmnet", glmnet.tumor)

## Decide coefficient number to take out (2:5)
coef.num <- 2

dir.create("./GSEA161123", showWarnings = FALSE)
foreach (i = seq_along(projects)) %dopar% {
    load(merge.data[i])  # Load methy and exp. data
    load(glmnet.normal[i])
    load(glmnet.tumor[i])  # Load glmnet results

    ## Get matrices of CpG islands
    CpG.normal <- getCpGs(normalGLM, coef.num)
    CpG.tumor <- getCpGs(tumorGLM, coef.num)
    
    ## fit linear model, get R^2 and p-value
    lm.normal <- buildLnr(CpG.normal, methyn, datan, p.value = T)
    model.normal <- lm.normal[1]
    R2.p.normal <- lm.normal[[2]]
    lm.tumor <- buildLnr(CpG.tumor, methyt, datat, p.value = T)
    model.tumor <- lm.tumor[1]
    R2.p.tumor <- lm.tumor[[2]]
    save(model.normal, R2.p.normal, model.tumor, R2.p.tumor, file = paste0("./GSEA161123/", projects[i], "_coef", coef.num, "_lm.RData"))
    print(paste0("Finished: ", projects[i]))
}
