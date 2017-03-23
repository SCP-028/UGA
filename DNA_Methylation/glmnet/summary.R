setwd("/lustre1/yz73026/array/glmnet_TF")
source("./GLM_func.R")
load("./data/annotation.RData")
tag.GLM <- "./results/moreTF/"
tag.summary <- "./results/moreTF/summary/"
minDev <- 0
GLM.normal <- list.files(tag.GLM, pattern = "normal")
summary.normal <- list.files(tag.summary,
                            pattern = "normal")
projects <- sub("(^\\w{4}).*", "\\1", GLM.normal)
normal.df <- vector("list", length(projects) * 4)

library(foreach)
library(doMC)
registerDoMC(12)
foreach (i = seq_along(projects)) %dopar% {
    load(paste0(tag.GLM, GLM.normal[i]))
    load(paste0(tag.summary, summary.normal[i]))
    df.list <- list(coef2, coef3, coef4, coef5)
    df.listname <- c("coef2", "coef3", "coef4", "coef5")
    names(df.list) <- df.listname
    for (j in seq_along(df.list)) {
        coef.df <- df.list[[df.listname[j]]]
        df <- retriveCoef(GLM=normalGLM, df=coef.df, minDev)
        df <- coefGroup(coef.list=df, annot=annot, coef.num=(j+1))
        df$cancer_type <- projects[i]
        df$sample_type <- "normal"
        df$coef_num <- df.listname[j]
        normal.df[[4*(i-1)+j]] <- df
    }
}
save(normal.df,
     file = paste0(tag.summary, "all_normal.RData"))
