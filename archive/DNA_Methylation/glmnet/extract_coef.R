# Extract GLM results of different coefs from glmnet models
try(setwd("C:/Users/jzhou/Desktop"))
try(setwd("/lustre1/yz73026/array/glmnet_TF"))
source("./GLM_func.R")
# normalGLM
filename <- list.files("./results", pattern = "normalGLM")
projects <- sub("(.*)\\.RData", "\\1", filename)
filename <- paste0("./results/", filename)
OUTPUT <- paste0("./results/summary/", projects)
rlist[annot, trans] <- readAnnot(methyPath = "./annot_removed.csv",
                              tfPath = "./trrust_rawdata.txt")

library(foreach)
library(DoMC)
registerDoMC(12)
foreach (i = seq_along(filename)) %dopar% {
    load(filename[i])
    results <- getSummary(normalGLM)
    coef2 <- detailSummary(normalGLM, 2, results, annot, trans)
    coef3 <- detailSummary(normalGLM, 3, results, annot, trans)
    coef4 <- detailSummary(normalGLM, 4, results, annot, trans)
    coef5 <- detailSummary(normalGLM, 5, results, annot, trans)
    save(results, coef2, coef3, coef4, coef5, file = OUTPUT[i])
}