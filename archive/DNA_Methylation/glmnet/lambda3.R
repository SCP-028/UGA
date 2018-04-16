# Processing of cv.glmnet data.
# Find lambdas that make temp == 3 and see how these CpG islands predict gene expression.

## Load required data
setwd("~/data/DNA_Methylation/Methy_array/RData/Glmnet")
rm(list = ls())
load("normalGLM.RData")
load("tumorGLM.RData")
load("Raw_Data/GLM.RData")

source("./get_R2_of_glmnet.R")
coefNum <- 3
## Normal and Tumor lists of Xs
CpG_nlist <- getCpGs(normalGLM, coefNum)
CpG_tlist <- getCpGs(tumorGLM, coefNum)

normal <- buildLnr(CpG_nlist, methyn, datan)
#normalModel <- normal[1]
normalR2 <- normal[[2]]
# time1 <- proc.time()
tumor <- buildLnr(CpG_tlist, methyt, datat)
# time2 <- proc.time() - time1
#tumorModel <- tumor[[1]]
tumorR2 <- tumor[[2]]

save(normalModel, normalR2, tumorModel, tumorR2, file = "./lmModel.RData")
