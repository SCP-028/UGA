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

getProject <- function(methy, expression){
#   Return project names that exist in both folders
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
    tf$PMID <- NULL
    tf <- tf[tf$Gene %in% annot$Gene_Name, ]
    tf <- tf[tf$TF %in% annot$Gene_Name, ]
    # TODO: aggregate TF name and Interaction
    return(list(annot, tf))
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
                       "GLM" = cv.glmnet(x,y))
            GLMdata[[i]] <- df
        }
    }
    names(GLMdata) <- rownames(expression)
    GLMdata <- Filter(Negate(function(x) is.null(unlist(x))),
                      GLMdata)
    return(GLMdata)
}
