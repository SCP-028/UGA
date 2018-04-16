# Uses glmnet results and calculate p-values and R^2-values for
# each gene, saving them in different .RData files

## Path and libraries
try(setwd("/lustre1/yz73026/array"))
rm(list = ls())
library(foreach)
library(doMC)

## Source codes and load data
source("./get_R2_of_glmnet.R")  # Code for using glmnet results
coef.num <- 2  # 2:5
glmnet.normal <- what.project(input = "./glmnet/",
                              output = "./GSEA161123/",
                              coef.num = coef.num,
                              type = "normal")
glmnet.normal <- glmnet.normal[order(glmnet.normal)]
glmnet.tumor <- sub("normalGLM", "tumorGLM", glmnet.normal)
projects <- sub("(.*)\\..*\\.RData", "\\1", glmnet.normal)
merge.data <- paste0("./mergeData/", projects, ".RData")
merge.data <- merge.data[order(merge.data)]
glmnet.normal <- paste0("./glmnet/", glmnet.normal)
glmnet.tumor <- paste0("./glmnet/", glmnet.tumor)

registerDoMC(length(projects))  # 12 projects, 4 coefs

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
    save(model.normal, R2.p.normal, model.tumor, R2.p.tumor,
         file = paste0("./GSEA161123/", projects[i],
                       "_coef", coef.num, "_lm.RData"))
    print(paste0("Finished: ", projects[i]))
}