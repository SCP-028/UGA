# This is for getting to know the methylation level of Tumor Supressor Genes, and see wheather there's strong correlation between methylation and gene expression levels.
setwd("~/data/DNA_Methylation/Methy_array/RData/Glmnet")
gene <- read.table("~/data/DNA_Methylation/Bisulfite_Seq/All_down_exp_TSGs_pan-cancer.txt")
load("../COAD.TN_noNA.RData")
methyn <- datan_final
methyt <- datat_final
rm(datan_final, datat_final)
annot <- read.table("annot_removed.txt", header = T)
library(dplyr)

ConvertToGeneName <- function(df, annotation)
{
    df <- df[rownames(df) %in% annotation$IlmnID, ]
    annotation <- annotation[annotation$IlmnID %in% rownames(df), ]
    annotation <- annotation[order(rownames(annotation)), ]
    df <- df[order(rownames(df)), ]
    df$IlmnID <- rownames(df)
    resultdf <- as.matrix(inner_join(annotation, df, by = "IlmnID"))  # inner_join saves duplicates.
    rownames(resultdf) <- resultdf[ ,1]
    resultdf <- cbind(resultdf[ ,2], resultdf[ ,4], resultdf[ ,grep("TCGA", colnames(resultdf))])
    colnames(resultdf)[1:2] <- c("Gene_Name", "Gene_Group")
    return(resultdf)
}

CpGTTest <- function(x,y)
{
    x <- x[order(rownames(x)), ]
    y <- y[order(rownames(y)), ]
    tResult <- matrix(nrow = nrow(x), ncol = 6, 
                    dimnames = list(rownames(x), c("Gene_Name", "Gene_Group", "Normal","Cancer","P","Methylation")))
    tResult[ ,1:2] <- x[ ,1:2]
    x <- x[ ,grep("TCGA", colnames(x))]
    class(x) <- "numeric"
    y <- y[, grep("TCGA", colnames(y))]
    class(y) <- "numeric"
    tResult[ ,3] <- rowMeans(x)
    tResult[ ,4] <- rowMeans(y)
    for(i in 1:nrow(x))
    {
        VAR <- var.test(x[i, ], y[i, ])$p.value
        if (VAR <= 0.05)
            tResult[i,5] <- t.test(x[i, ], y[i, ], var.equal = F, paired = F)$p.value
        else
            tResult[i,5] <- t.test(x[i, ], y[i, ], var.equal = T, paired = F)$p.value
        if (tResult[i,5] <= 0.05)
            tResult[i,6] <- ifelse (tResult[i,4] >= tResult[i,3], "Hyper", "Hypo")
        else
            tResult[i,6] <- "-"
    }
    return(tResult)
}

annot <- annot[annot$Gene_Name %in% gene$V2, ]
methyn <- ConvertToGeneName(as.data.frame(methyn), annot)
methyt <- ConvertToGeneName(as.data.frame(methyt), annot)

projects <- c("BLCA", "BRCA", "COAD", "HNSC", "KICH", "KIRC", "LUAD", "LUSC", "PRAD", "THCA", "UCEC")
# for (i in projects)
# {
#     temp <- gene[grep(i, gene$V4), ]
#     x <- methyn[methyn[ ,1] %in% gene$V2, ]
#     y <- methyt[methyt[ ,1] %in% gene$V2, ]
#     tResult <- paste("t.result", i, sep = "_")
#     assign(tResult, CpGTTest(methyn, methyt))
# }
for (i in projects)
{
    temp <- gene[grep(i, gene$V4), ]
    TResult <- tresult[tresult$Gene_Name %in% temp$V2, ]
    hyper <- count(TResult[TResult$Methylation == "Hyper", ], "Gene_Group")
    hypo <- count(TResult[TResult$Methylation == "Hypo", ], "Gene_Group")
    countResult <- cbind(hyper, hypo$freq)
    colnames(countResult) <- c("Gene_Group", "Hyper", "Hypo")
    dir.create(file.path(paste("TCGA-", i, sep = "")))
    write.table(TResult, paste("./TCGA-", i, "/t_test.txt", sep = ""))
    write.table(countResult, paste("./TCGA-", i, "/gene_group.txt", sep = ""))
}


157 96