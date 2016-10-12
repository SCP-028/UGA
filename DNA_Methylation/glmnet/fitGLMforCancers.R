# Pre-process methy and GX datasets of different cancers,
# and fit GLM model.
library(dplyr)
library(mygene)
library(glmnet)
try(setwd("/home/yi/data/DNA_Methylation/Methy_array/RData/"))
try(setwd("/lustre1/yz73026/array"))

getProject <- function(methy, expression){
#   Return project names that exist in both folders
    x <- sub("(.*)\\.TN.*", "\\1", methy)
    y <- sub("(.*)\\.up.*", "\\1", expression)
    z <- x[x %in% y]
    return(z)
}

convertName <- function(df){
#   Convert expression df names from symbol|entrez to symbol
    library(mygene)
    library(dplyr)
    rownames(df) <- sub(".*\\|(\\d*)", "\\1", rownames(df))
    temp <- as.data.frame(queryMany(rownames(df),
                                    scopes = "entrezgene",
                                    fields = "symbol",
                                    species = "human"))
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

pairNames <- function(methy, expression, annot){
#   Keep cases and genes that appear in both methy and expression,
#   and rows in methy that appear in annot.
    colnames(methy) <- sub("(TCGA)-(..)-(\\w{4}).*",
                            "\\1_\\2_\\3",colnames(methy))
    methy <- methy[complete.cases(methy), ]
    expression <- expression[complete.cases(expression), ]
    methy <- methy[ ,colnames(methy) %in% colnames(expression)]
    methy <- methy[rownames(methy) %in% annot$IlmnID, ]
    methy <- methy[ ,!duplicated(colnames(methy))]
    methy <- methy[ ,order(colnames(methy))]
    expression <- expression[ ,colnames(expression) %in% 
                                colnames(methy)]
    expression <- expression[rownames(expression) %in% 
                                annot$Gene_Name, ]
    output <- list(methy, expression)
    return(output)
}

fitGLM <- function(methy, expression, annot){
#   Use glmnet to fit GLM row by row
    library(glmnet)
    GLMdata <- list()
    pb <- txtProgressBar(min = 0, max = nrow(expression),
                        char = "#", style = 3)
    for(i in 1:nrow(expression)){
        y <- expression[i, , drop = F]
        CpG <- annot[annot$Gene_Name %in% rownames(y),1]
        x <- t(methy[rownames(methy) %in% CpG, , drop = F])
        if (is.null(dim(x)) | ncol(x) <= 1){
            GLMdata[[i]] <- 0
        }
        else {
            GLMdata[[i]] <- cv.glmnet(x,y)
        }
        setTxtProgressBar(pb, i)        
    }
    return(GLMdata)
}

##########
# Find files to load
methyFiles <- list.files("./methylation450", pattern = "\\.RData")
expressionFiles <- list.files("./expression", pattern = "\\.RData")
projects <- getProject(methyFiles, expressionFiles)
methyFiles <- paste("./methylation450/", projects, 
                ".TNMethylation450.RData", sep = "")
expressionFiles <- paste("./expression/", projects, 
                ".update1.genes.normalized.RData", sep = "")
resultFiles <- paste("./glmnet/", projects, sep = "")

# Remove duplicated IlmnID, Gene name, and Gene group
annot <- read.csv("annot_removed", sep = "\t")
annot <- anti_join(annot, annot[duplicated((annot[ ,c(1,2,4)])), ])

for (i in seq_along(projects)){
#   Different tumor types
    load(methyFiles[i])
    if(ncol(datan) != 0 && ncol(datat) != 0){
#       Some data seem to be missing
        methyn <- datan
        methyt <- datat
        load(expressionFiles[i])
        datan <- convertName(datan)
        datat <- convertName(datat)
        normal <- pairNames(methyn, datan, annot)
        methyn <- as.matrix(normal[[1]])
        datan <- as.matrix(normal[[2]])
        tumor <- pairNames(methyt, datat, annot)
        methyt <- as.matrix(tumor[[1]])
        datat <- as.matrix(tumor[[2]])
        rm(normal, tumor)
        normalGLM <- fitGLM(methyn, datan, annot)
        save(normalGLM, file = paste(resultFiles[i],
                                    ".normalGLM.RData", sep = ""))
        tumorGLM <- fitGLM(methyt, datat, annot)
        save(tumorGLM, file = paste(resultFiles[i],
                                    ".tumorGLM.RData", sep = ""))
    }
}
