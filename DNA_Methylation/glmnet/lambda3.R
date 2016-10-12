# Processing of cv.glmnet data.
# Find lambdas that make temp == 3 and see how these CpG islands predict gene expression.

## Load required data
setwd("~/data/DNA_Methylation/Methy_array/RData/Glmnet")
rm(list = ls())
load("normalGLM.RData")
load("tumorGLM.RData")

load("Raw_Data/GLM.RData")

## Define function that greps 3 Xs and store it as a matrix in a list
getNzero <- function(listInList, coefNum)
{
    if (length(listInList) != 10)  # cor.test or just 0
        return ("Remove")
    else
    {
        nzero <- unlist(listInList[6])  # nzero list
        Beta <- as.matrix(listInList$glmnet.fit$beta) # CpG island matrix
        nzero <- names(nzero[nzero == coefNum])  # s(\d*)
        if (length(nzero) >= 1)  # More than 1 set of CpG islands
        {
            nzero <- sub("nzero.", "", nzero)
            nzero <- nzero[length(nzero)]  # Only keep the one with best performance
            CpG <- Beta[ ,colnames(Beta) %in% nzero]
            CpG <- as.matrix(CpG[CpG != 0])  # Same as last if
            colnames(CpG) <- nzero
            return(CpG)
        }
        else  # No lambda to make x == coefNum
            return("Remove")
    }
}

## Get the Matrices
getCpGs <- function(df, coefNum)
{
    CpG_list <- list()
    for (i in 1:length(df))  # List too complicated, use loop rather than apply <- Room for improvement!
    {
        x <- df[[i]]
        CpG_list[[i]] <- getNzero(x, coefNum)
    }
    return(CpG_list)
}

## Linear models of Gene Expression ~ remaining CpG islands
buildLnr <- function(df, methy, expression, annot)
{
    r2df <- data.frame(rownames(expression), rep(0, nrow(expression)))
    colnames(r2df) <- c("Gene_Symbol", "R_Square") # Build data frame to store R^2
    lmodel <- list()  # List to store linear model results
    pb <- txtProgressBar(min = 0, max = length(df), char = "#", style = 3)
    for (i in 1:length(df))
    {
        temp <- df[[i]]
        output <- getR2(i, temp, methy, expression, annot)
        lmodel[[i]] <- output[1]
        r2df[i,2] <- output[2]
        setTxtProgressBar(pb, i)
    }
    r2df <- r2df[complete.cases(r2df$R_Square), ]
    return(list(lmodel, r2df))
}

## Fit linear model and calculate R^2
getR2 <- function(i, temp, methy, expression, annot)
{
    if (class(temp) == "character" | is.na(expression[i,1]))  # # Remove genes that can't be predicted
    {
        return(list(0, NA))
    }
    if (class(temp) == "matrix")  # One or more CpG islands
    {
        x <- methy[rownames(methy) %in% rownames(temp), , drop = F]  # CpG islands
        y <- as.data.frame(expression[i, , drop = F])  # Gene expression data
        fitlm <- as.data.frame(t(rbind(y,x)))
        fitlm <- lm(fitlm)  # Linear model
        rSquare <- summary(fitlm)$r.squared  # R^2
        return(list(fitlm, rSquare))  # Return values
    }
}
## Normal and Tumor lists of Xs
CpG_nlist <- getCpGs(normalGLM, coefNum)
CpG_tlist <- getCpGs(tumorGLM, coefNum)

normal <- buildLnr(CpG_nlist, methyn, datan, annot)
#normalModel <- normal[1]
normalR2 <- normal[[2]]
# time1 <- proc.time()
tumor <- buildLnr(CpG_tlist, methyt, datat, annot)
# time2 <- proc.time() - time1
#tumorModel <- tumor[[1]]
tumorR2 <- tumor[[2]]

save(normalModel, normalR2, tumorModel, tumorR2, file = "./lmModel.RData")
