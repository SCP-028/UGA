rlist <- structure(NA,class="result")
# Usage example: rlist[x, y] <- some.function.that.returns.2.values()
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
#   df: gene expression data.frames
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
    library(glmnet, quietly = T)
    GLMdata <- vector("list", nrow(expression))
    # pb <- txtProgressBar(min = 0, max = nrow(expression),
                         # char = "#", style = 3)
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
        # setTxtProgressBar(pb, i)
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
    library(glmnet, quietly = T)
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
    return(as.character(CpG.list))
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
        x <- x[!duplicated(x[ ,c("IlmnID", "Gene_Name")]), ]
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
    colnames(result) <- c("gene_symbol",
                          paste0("coef", predictor, "_CpG"),
                          paste0("coef", predictor, "_%dev"),
                          paste0("coef", predictor, "_lambda"))
    rownames(result) <- result$gene_symbol <- names(GLMlist)
    result[ ,3] <- resultdf[ ,paste0("coef", predictor, "_%dev")]
    result[ ,4] <- resultdf[ ,paste0("coef", predictor, "_lambda")]
    result[ ,2] <- getCpG(GLMlist, predictor, resultdf)
    result <- calcGroup(result, resultdf, predictor, annotation)
    result$Row.names <- NULL
    return(result)
}

finalGroup <- function(filepath, minDev = 0, existTF = F){
#   Count all results' methylation groups.
#   minDev: minimum value of %dev, remove all values smaller than it.
#   existTF: if TRUE, remove all genes that don't have a regulating TF.
    library(dplyr, quietly = T)
    filepaths <- list.files(filepath, pattern = "GLM\\.RData$",
                            full.names = T)
    filenames <- list.files(filepath, pattern = "GLM\\.RData$")
    cancer_type <- vector("character", 672)  # 12 cancer types
    sample_type <- rep(c(rep("normal", 28), rep("tumor", 28)), 12)
    coef_num <- rep(c(rep(2, 7), rep(3, 7), rep(4, 7), rep(5, 7)), 24)
    methy_group <- rep(c("1stExon", "3'UTR", "5'UTR", "Body", "TSS1500",
                         "TSS200", "TFmethy"), 96)
    group_count <- vector("numeric", 672)
    for (i in seq_along(filenames)){
        load(filepaths[i])
        cancer_type[(28*i-27):(28*i)] <- sub("(^\\w{4}).*", "\\1", filenames[i])
        if (existTF){
            geneName <- results %>%
                            filter(transcription_factor != "") %>%
                            select(gene_symbol)
            geneName <- as.character(geneName$gene_symbol)
            group_count[(28*i-27):(28*i-21)] <- colSums(coef2 %>% 
                                                            filter(gene_symbol %in% geneName) %>%
                                                            filter(`coef2_%dev` >= minDev) %>%
                                                            select(-contains("_")),
                                                        na.rm = T)
            group_count[(28*i-20):(28*i-14)] <- colSums(coef3 %>% 
                                                            filter(gene_symbol %in% geneName) %>% 
                                                            filter(`coef3_%dev` >= minDev) %>%
                                                            select(-contains("_")),
                                                        na.rm = T)
            group_count[(28*i-13):(28*i-7)] <- colSums(coef4 %>% 
                                                            filter(gene_symbol %in% geneName) %>% 
                                                            filter(`coef4_%dev` >= minDev) %>%
                                                            select(-contains("_")),
                                                        na.rm = T)
            group_count[(28*i-6):(28*i)] <- colSums(coef5 %>% 
                                                            filter(gene_symbol %in% geneName) %>% 
                                                            filter(`coef5_%dev` >= minDev) %>%
                                                            select(-contains("_")),
                                                        na.rm = T)
        }
        else{
            group_count[(28*i-27):(28*i-21)] <- colSums(coef2 %>%
                                                            filter(`coef2_%dev` >= minDev) %>%
                                                            select(-contains("_")),
                                                        na.rm = T)
            group_count[(28*i-20):(28*i-14)] <- colSums(coef3 %>%
                                                            filter(`coef3_%dev` >= minDev) %>%
                                                            select(-contains("_")),
                                                        na.rm = T)
            group_count[(28*i-13):(28*i-7)] <- colSums(coef4 %>%
                                                            filter(`coef4_%dev` >= minDev) %>%
                                                            select(-contains("_")),
                                                        na.rm = T)
            group_count[(28*i-6):(28*i)] <- colSums(coef5 %>%
                                                            filter(`coef5_%dev` >= minDev) %>%
                                                            select(-contains("_")),
                                                        na.rm = T)
        }
    }
    result <- data.frame(cancer_type, sample_type, coef_num, methy_group,
                         group_count)
    return(result)
}

plotFig <- function(df, filename){
#   Barplot to get a general picture of the distribution of the data.
#   df: the output of function finalGroup.
    library(ggplot2, quietly = T)
    ggplot(df, aes(x=cancer_type, y=group_count))+
    geom_bar(stat="identity", aes(fill=methy_group), position="dodge")+
    facet_grid(coef_num~sample_type, labbller=label_both)
    ggsave(filename, dpi = 200, width = 20, height = 45, units = "in")
}

retriveCoef <- function(GLM, df, minDev){
#   Extract each gene's correlation coefficient (positive or negative)
#   GLM is normalGLM or tumorGLM list
#   df is coef2 - coef5 data.frame
    library(dplyr, quietly = T)
    library(glmnet, quietly = T)
    df <- df[df[ ,grep("dev", colnames(df))] >= minDev, ]
    x <- df %>%
            filter(!is.na(gene_symbol)) %>%
            select(gene_symbol, contains("lambda"))
    x <- x[rowSums(is.na(x)) == 0, ]
    x <- x[order(x$gene_symbol), ]
    y <- x[ ,2]
    GLM <- lapply(GLM, function(x) x$GLM$glmnet.fit)
    GLM <- GLM[names(GLM) %in% x$gene_symbol]
    GLM <- GLM[order(names(GLM))]
    result <- vector("list", length(GLM))

    # pb <- txtProgressBar(min = 0, max = length(GLM), char = "#", style = 3)
    for (i in seq_along(GLM)){
        temp <- as.matrix(coef(GLM[[i]], s = y[i]))
        result[[i]] <- temp[(temp[ ,1] != 0) &
                            (rownames(temp) != "(Intercept)"), ]
        # setTxtProgressBar(pb,i)
    }
    names(result) <- names(GLM)
    result <- Filter(Negate(function(x) is.null(unlist(x))), result)
    return(result)
}

coefGroup <- function(coef.list, annot, coef.num){
#   Determine which group the CpG island belongs to.
#   coef.list: the output of function retriveCoef.
    library(dplyr, quietly = T)
    gene_symbol <- rep(names(coef.list), each = coef.num)
    ilmnID <- lapply(coef.list, function(x) names(x))
    ilmnID <- unlist(ilmnID)
    names(ilmnID) <- gene_symbol
    subAnnot <- annot %>%
                select(IlmnID, Gene_Name, Gene_Group) %>%
                filter(IlmnID %in% ilmnID)
    coef_value <- unlist(coef.list)
    methy_group <- vector("character", length(ilmnID))
    # pb <- txtProgressBar(min = 0, max = length(ilmnID),
                         # char = "#", style = 3)
    for (i in seq_along(ilmnID)){
        x <- subAnnot %>%
                filter(IlmnID == ilmnID[i]) %>%
                filter(Gene_Name == gene_symbol[i])
        methy_group[i] <- ifelse(nrow(x) == 1,
                                 as.character(x$Gene_Group), "TFmethy")
        # setTxtProgressBar(pb, i)
    }
    result <- data.frame(gene_symbol=gene_symbol,
                         ilmnID=ilmnID,
                         methy_group=methy_group,
                         coef_value=coef_value)
    return(result)
}

coefDist <- function(df){
#   See distribution of function coefGroup results.
    library(ggplot2)
    ggplot(df, aes(x=coef_value))+
    geom_density(aes(group=methy_group, colour=methy_group,
                     fill=methy_group),
                 alpha=0.3)
    # facet_grid(coef_num~sample_type, labbller=label_both)
    # , filename
    # ggsave(filename, dpi = 200, width = 20, height = 45, units = "in")
}