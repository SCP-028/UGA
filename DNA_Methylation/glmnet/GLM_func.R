rlist <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
   args <- as.list(match.call())
   args <- args[-c(1:2,length(args))]
   length(value) <- length(args)
   for(i in seq(along=args)) {
     a <- args[[i]]
     if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
   }
   x
}

####### BUILD GLM #######
getProject <- function(methy, expression){
#   Return project names that exist in both folders (the 2 args)
    x <- sub("(.*)\\.TN.*", "\\1", methy)
    y <- sub("(.*)\\.up.*", "\\1", expression)
    z <- x[x %in% y]
    return(z)
}

readAnnot <- function(methyPath="./glmnet_TF/annot_removed.csv",
                     tfPath="./glmnet_TF/trrust_rawdata.txt"){
#   Add transcription factor info to annotation file, and
#   clean duplicated records.
    library(dplyr)
    # Read annotation files
    annot <- read.csv(methyPath, header=T, sep="\t")
    tf <- read.csv(tfPath, header=T, sep="\t")
    annot <- anti_join(annot, annot[duplicated((annot[ ,c(1,2,4)])), ])
    # Add level to factor
    x <- annot$Gene_Group
    levels(x) <- c(levels(x), "TFmethy")
    annot$Gene_Group <- x
    tf$PMID <- NULL
    tf <- tf[tf$Gene %in% annot$Gene_Name, ]
    tf <- tf[tf$TF %in% annot$Gene_Name, ]
    # TODO: aggregate TF name and Interaction
    return(list(annot, tf))
}

addTRED <- function(filePath = "./TRED.csv"){
#   Add info from the TRED database, return data frame.
#   filePath: TRED.csv file
    df <- read.csv(filePath, sep = "\t", header = T)
    s <- strsplit(as.character(df$Gene), split = "; ")
    df <- data.frame(TF = rep(df$TF, sapply(s, length)), Gene = unlist(s))
    return(df)
}

addNEPH <- function(filePath = "./Neph2012/"){
#   Add info from the Neph2012 database, return data frame.
#   filePath: folder containing 41 cell types
    pathList <- list.files(filePath, pattern = "\\.genes$", recursive = T)
    pathList <- paste0(filePath, pathList)
    df <- read.csv(pathList[1], header = F, sep = "\t")
    for (i in seq_along(pathList)){
        temp <- read.csv(pathList[i], header = F, sep = "\t")
        df <- merge(df, temp, all = T)
    }
    colnames(df) <- c("TF", "Gene")
    return(df)
}

addENCODE <- function(){
#   Use lists from tftargets to generate data frame of ENCODE results.
    library(tftargets)
    library(mygene)
    df <- data.frame(TF = names(ENCODE),
                     Gene = sapply(ENCODE, function(x) paste(x, collapse = "; ")))
    s <- strsplit(as.character(df$Gene), split = "; ")
    df <- data.frame(TF = rep(df$TF, sapply(s, length)), Gene = unlist(s))
    temp <- as.data.frame(queryMany(df$Gene,
                                    scopes = "entrezgene",
                                    fields = "symbol",
                                    species = "human"))
    df$Gene <- temp$symbol
    df <- df[complete.cases(df$Gene), ]
    return(df)
}

addMARBACH <- function(){
#   Use lists from tftargets to generate data frame of Marbach2016 results.
    library(tftargets)
    df <- data.frame(TF = names(Marbach2016),
                     Gene = sapply(Marbach2016, function(x) paste(x, collapse = "; ")))
    s <- strsplit(as.character(df$Gene), split = "; ")
    df <- data.frame(TF = rep(df$TF, sapply(s, length)), Gene = unlist(s))
    return(df)
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
    methy <- methy[ ,unique(colnames(methy))]
    methy <- methy[ ,order(colnames(methy))]
    expression <- expression[ ,colnames(expression) %in% 
                                colnames(methy)]
    expression <- expression[rownames(expression) %in% 
                                annot$Gene_Name, ]
    return(list(methy, expression))
}

fitGLM <- function(methy, expression, annotation, trans){
#   Use glmnet to fit GLM
#   methy and expression are 2 data.frames, and annotation and trans
#   are the annotation files for matching CpG islands to genes, and 
#   transcription factors to genes, respectively.
    library(glmnet)
    GLMdata <- vector("list", nrow(expression))
    for(i in 1:nrow(expression)){
        # Gene expression
        y <- expression[i, , drop = F]
        # Add transcription factor name(s)
        z <- as.character(trans$TF[trans$Gene == rownames(y)])
        z <- c(rownames(y), z)
        CpG <- as.character(annotation[annotation$Gene_Name %in% z,1])
        x <- t(methy[rownames(methy) %in% CpG, , drop = F])
        if (is.null(dim(x)) | ncol(x) <= 1){
            next
        }
        else {
            df <- list("Gene" = rownames(y),
                       "TF" = z[z != rownames(y)],
                       "CpG site" = CpG,
                       "GLM" = cv.glmnet(x,y, parallel = T))
            GLMdata[[i]] <- df
        }
    }
    names(GLMdata) <- rownames(expression)
    GLMdata <- Filter(Negate(function(x) is.null(unlist(x))),
                      GLMdata)
    return(GLMdata)
}

####### EXTRACT RESULTS #######
buildDf <- function(row.num, col.num = 12){
#   Build data.frame to store all results
    df <- data.frame(matrix(nrow = row.num, ncol = col.num))
    return(df)
}

extractTF <- function(df){
#   Return a list of transcription factors in the same order as the input
#   data frame's gene order.
    df.tf <- lapply(df, function(x) paste(x$TF, collapse = ";"))
    return(df.tf)
}

extractLambda <- function(df, func = "lambda.min"){
#   Return a list of lambda.min or lambda.1st.
    if (func == "lambda.min"){
        df <- lapply(df, function(x) x$GLM$lambda.min)
    }
    if (func == "lambda.1se"){
        df <- lapply(df, function(x) x$GLM$lambda.1se)
    }
    return(df)
}

extractCoef <- function(df, predictor = 2){
#   Return a data frame of %dev and lambda values.
    df <- lapply(df, function(x) x$GLM$glmnet.fit)
    # predictor num, %dev, and lambda
    df <- lapply(df, function(x) cbind(x$df, x$dev.ratio, x$lambda))
    df <- lapply(df, function(x) x[x[ ,1] == predictor, , drop = F])
    df <- lapply(df, function(x) as.data.frame(x))
    resultdf <- data.frame(matrix(nrow = length(df), ncol = 2))
    for (i in seq_along(df)){
        temp <- df[[i]]
        if (nrow(temp) != 0){
            temp <- unlist(temp[nrow(temp), ])
            resultdf[i, ] <- tail(temp, 2)
        }
    }
    colnames(resultdf) <- c(paste0("coef", predictor, "_%dev"),
                            paste0("coef", predictor, "_lambda"))
    rownames(resultdf) <- names(df)
    return(resultdf)
}

getSummary <- function(df){
#   Summarize different model performances with 2->5 predictors.
    result <- buildDf(row.num = length(df), col.num = 4)
    colnames(result) <- c("gene_symbol", "transcription_factor",
                          "lambda_1se", "lambda_min")
    result$gene_symbol <- names(df)
    result$transcription_factor <- extractTF(df)
    result$lambda_1se <- extractLambda(df, func = "lambda.1se")
    result$lambda_min <- extractLambda(df, func = "lambda.min")
    for (i in 2:5){
    #   2 to 5 predictors
        temp <- extractCoef(df, predictor = i)
        result <- cbind.data.frame(result, temp)
    }
    return(result)
}


getCpG <- function(GLMlist, predictor, resultdf){
#   get CpG islands' names (and groups)
    library(glmnet)
    predictor <- paste0("coef", predictor, "_lambda")
    CpG.list <- vector("list", length(GLMlist))
    for (i in seq_along(CpG.list)){
        # Extract model coefficients from GLM objects, s is lambda
        tempLambda <- resultdf[predictor][i,1]
        if (!is.na(tempLambda)) {        
            temp <- as.matrix(coef(GLMlist[[i]]$GLM,
                              s = tempLambda))
            if (ncol(temp) != 0){
                temp <- names(temp[temp[ ,1] != 0, ])
                temp <- temp[temp != "(Intercept)"]
                CpG.list[i] <- paste(temp, collapse = ";")
            }
        }
    }
    return(CpG.list)
}

getGenelist <- function(genename, resultdf){
#   Return a gene's symbol and its (if any) TF.
    x <- resultdf[rownames(resultdf) == genename, ]
    x <- c(x$gene_symbol,
           unlist(strsplit(unlist(x$transcription_factor), ";")))
    return(x)
}

calcGroup <- function(result, resultdf, predictor, annotation){
#   Calculate methylation groups for predictors.
    predictor <- paste0("coef", predictor, "_CpG")
    cpg <- result[ ,2]
    names(cpg) <- result$gene_symbol
    cpg <- unlist(cpg, recursive = F)
    cpg <- lapply(cpg, function(x) unlist(strsplit(x, ";")))
    df <- buildDf(row.num = length(cpg), col.num = 7)
    colnames(df) <- c("1stExon", "3'UTR", "5'UTR", "Body", "TSS1500",
                      "TSS200", "TFmethy")
    rownames(df) <- names(cpg)
    for (i in seq_along(cpg)){
        x <- annot[annot$IlmnID %in% cpg[[i]], ]
        x <- x[x$Gene_Name %in% getGenelist(names(cpg)[i], resultdf), ]
        x$Gene_Group[x$Gene_Name != names(cpg)[i]] <- "TFmethy"
        df[i, ] <- t(as.matrix(summary(x$Gene_Group)))
    }
    df <- merge(result, df, by = "row.names", all.x = T)
    return(df)
}

detailSummary <- function(GLMlist, predictor, resultdf,
                          annotation, transAnnot){
#   GLMlist is normalGLM or tumorGLM, used to extract coef and other info,
#   predictor is a number between 2 and 5,
#   resultdf is for less computation,
#   annot and trans are annotation files.
    result <- buildDf(row.num = length(GLMlist), col.num = 4)
    colnames(result) <- c("gene_symbol", "coef2_CpG", "coef2_%dev",
                          "coef2_lambda")
    rownames(result) <- result$gene_symbol <- names(GLMlist)
    result$"coef2_%dev" <- resultdf$"coef2_%dev"
    result$coef2_lambda <- resultdf$coef2_lambda
    result$coef2_CpG <- getCpG(GLMlist, predictor, resultdf)
    result <- calcGroup(result, resultdf, predictor, annotation)
    result$Row.names <- NULL
    return(result)
}