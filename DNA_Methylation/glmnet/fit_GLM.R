# Pre-process of methylation and gene expression data,
# and fit them into a GLM model.

setwd("~/data/DNA_Methylation/Methy_array/RData/")
## Methy data colnames -> GX colnames
load("./COAD.TN_noNA.RData")
methyn <- datan_final
methyt <- datat_final
rm(datan_final,datat_final)
colnames(methyn) <- sub("(TCGA)-(..)-(\\d{4}).*","\\1_\\2_\\3",colnames(methyn))
colnames(methyt) <- sub("(TCGA)-(..)-(\\d{4}).*","\\1_\\2_\\3",colnames(methyt))

## Annotation and GX data
setwd("./Glmnet/")
annot <- read.csv("annot_removed",sep = "\t")
load("~/data/DNA_Methylation/Methy_array/RData/COAD.update1.genes.normalized.RData")
## Pre-process annotation data to remove duplicates
library(dplyr)
dup <- annot[duplicated((annot[ ,c(1,2,4)])), ]
annot <- anti_join(annot,dup)
rm(dup, all.genes)
## GX data rownames -> Gene names
library(mygene)
convertName <- function(df)
{
    rownames(df) <- sub(".*\\|(\\d*)", "\\1", rownames(df))
    temp <- as.data.frame(queryMany(rownames(df),scopes = "entrezgene",fields = "symbol",species = "human"))
    temp <- temp[complete.cases(temp$symbol), ]
    df <- as.data.frame(df)
    df$query <- rownames(df)
    temp <- inner_join(df, temp, by = "query")
    rownames(temp) <- temp$symbol
    df <- temp[ ,1:ncol(df) - 1]
    df <- df[order(rownames(df)), ]
    df <- df[ ,order(colnames(df))]
    return(df)
}

datan <- convertName(datan)
datat <- convertName(datat)

## Remove unpaired cases and unannotated CpG islands
pairNames <- function(methy, expression, annot)
{
    methy <- methy[ ,colnames(methy) %in% colnames(expression)]
    methy <- methy[rownames(methy) %in% annot$IlmnID, ]
    methy <- methy[ ,!duplicated(colnames(methy))]
    methy <- methy[ ,order(colnames(methy))]
    expression <- expression[ ,colnames(expression) %in% colnames(methy)]
    expression <- expression[rownames(expression) %in% annot$Gene_Name, ]
    output <- list(methy, expression)
    return(output)
}
normal <- pairNames(methyn, datan, annot)
methyn <- as.matrix(normal[[1]])
datan <- as.matrix(normal[[2]])
tumor <- pairNames(methyt, datat, annot)
methyt <- as.matrix(tumor[[1]])
datat <- as.matrix(tumor[[2]])
rm(normal, tumor)

## Calculate with glmnet package
library(glmnet)
fitGLM <- function(methy, expression, annot)
{
    GLMdata <- list()
    for(i in 1:nrow(expression))
    {
        # print(i / nrow(expression))
        y <- expression[i, , drop = F]
        CpG <- annot[annot$Gene_Name %in% rownames(y),1]
        x <- t(methy[rownames(methy) %in% CpG, , drop = F])
        if (is.null(dim(x)) | ncol(x) <= 1)
        {
            GLMdata[[i]] <- 0
        }
        else
        {
            GLMdata[[i]] <- cv.glmnet(x,y)
        }
        
    }
    return(GLMdata)
}

normalGLM <- fitGLM(methyn, datan, annot)
save(normalGLM, file = "./normalGLM.RData")
tumorGLM <- fitGLM(methyt, datat, annot)
save(tumorGLM, file = "./tumorGLM.RData")