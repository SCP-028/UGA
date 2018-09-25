merge.result <- function(pathways, directory = "./"){
    library(dplyr)
    resultFiles <- list.files(paste(directory, sep = ""),
                            pattern = "\\.txt")
    fileNames <- sub("(.*)\\.txt", "\\1", resultFiles)
    resultFiles <- paste(directory, resultFiles, sep = "")
    df <- data.frame(annotation = names(pathways))
    df$annotation <- as.character(df$annotation)
    pb <- txtProgressBar(min = 0, max = length(resultFiles),
                         style = 3, char = "#")
    for (i in seq_along(resultFiles)){
        temp <- data.table::fread(resultFiles[i], sep = " ",
                                header = F, data.table = F)
        colnames(temp) <- c("annotation", "ES", "R2")
        temp$annotation <- as.character(temp$annotation)
        temp <- temp[order(temp$annotation, decreasing = T), ]
        temp <- temp[!duplicated(temp$annotation), ]
        df <- left_join(df, temp, by = "annotation")
        colnames(df)[(ncol(df) - 1):ncol(df)] <- c(paste(fileNames[i], "ES", sep = "-"), paste(fileNames[i], "R2", sep = "-"))
        setTxtProgressBar(pb, i)
        gc()
    }
    return(df)
}

setwd("~/data/DNA_Methylation/Methy_array/RData")
load("./piano/pathways.list.GOMaco.CORUM.RCT.HGNC.Msigdb.gsc.RData")
pathways <- lapply(gsc[[1]], function(x) sub("(.*)\\|.*", "\\1", x))
pathways <- pathways[!duplicated(names(pathways))]
noR2.df <- merge.result(pathways, directory = "./GSEA/GSEA_no/")
withR2.df <- merge.result(pathways, directory = "./GSEA/GSEA_with/")
save(noR2.df, withR2.df, file = "./GSEAResults.RData")
