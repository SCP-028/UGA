subsetBy <- function(df, project = NULL, type = c("normal", "tumor"),
                      coef = NULL, annotation = NULL, cutoff = NULL,
                      score = c("ES", "R2")){
#   Return a subset of the original data frame
    if (!is.null(project)){
        project <- paste0("(", project, ")", collapse = "|")
        tag.project <- grep(ignore.case = T, project, colnames(df))
        df <- df[ ,c(1, tag.project)]
    }
    type <- paste0("(", type, ")", collapse = "|")
    score <- paste0("(", score, ")", collapse = "|")
    df <- df[ ,c(1, grep(ignore.case = T, type, colnames(df)))]
    df <- df[ ,c(1, grep(ignore.case = F, score, colnames(df)))]
    if (!is.null(coef)){
        coef <- paste0(coef, collapse = "")
        coef <- paste0("[", coef, "]", "-")
        tag.coef <- grep(ignore.case = T, coef, colnames(df))
        df <- df[ ,c(1, tag.coef)]
    }
    if (!is.null(annotation)){
        annotation <- paste0("(", annotation, ")", collapse = "|")
        tag.annotation <- grep(ignore.case = T, annotation, df$annotation)
        df <- df[tag.annotation, ]
    }
    if (!is.null(cutoff)){
#   cutoff = abs(tumor - normal)
        
    }
    return(df)
}

getGenes <- function(keyword, gsc = NULL){
#   Get gene symbols from specific pathways
    if (missing(keyword)){
        stop("Please input the keyword(s) you want to extract.")
    }
    if (missing(gsc) | length(gsc) != 2){
        load("~/data/DNA_Methylation/Methy_array/RData/piano/pathways.list.GOMaco.CORUM.RCT.HGNC.Msigdb.gsc.RData")
    }
    pathways <- gsc[[1]]
    keyword <- paste0("(", keyword, ")", collapse = "|")
    pathways <- pathways[grep(keyword, names(pathways))]
    pathways <- unlist(pathways)
    pathways <- data.frame( pathway = sub("(.*)\\..*", "\\1", names(pathways)),
                            uniprotKB = sub(".*\\.(.*)", "\\1", names(pathways)),
                            entrez = sub(".*\\|(.*)", "\\1", pathways),
                            symbol = sub("(.*)\\|.*", "\\1", pathways))
    rownames(pathways) <- 1:nrow(pathways)
    try(return(pathways))
}
