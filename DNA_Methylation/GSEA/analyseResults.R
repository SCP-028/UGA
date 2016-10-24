subsetBy <- function(df, project = NULL, type = c("normal", "tumor"),
                      annotation = NULL, score = c("ES", "R2"),
                      coef = NULL){
#   Return a subset of the original data frame
    if (!is.null(project)){
        project <- paste0("(", project, ")", collapse = "|")
        tag.project <- grep(ignore.case = T, project, colnames(df))
        df <- df[ ,c(1, tag.project)]
    }
    type <- paste0("(", type, ")", collapse = "|")  # Choose type
    df <- df[ ,c(1, grep(ignore.case = T, type, colnames(df)))]
    score <- paste0("(", score, ")", collapse = "|")  # Choose score type
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
                            symbol = sub("(.*)\\|.*", "\\1", pathways) )
    rownames(pathways) <- 1:nrow(pathways)
    try(return(pathways))
}

compareSameProj <- function(df, project, cutoff = 0.5, type){
#   Compare differences between normal and tumor samples
#   cutoff = abs(tumor - normal)
#   type: tumor - Tumor better than normal
#         normal - Normal better than tumor
    df <- subsetBy(df, project = project)
    tag.ES <- grep("ES", colnames(df))
    x.ES <- grep( paste(type, "ES", sep = ".*"), colnames(df) )
    y.ES <- tag.ES[!(tag.ES %in% x.ES)]  # Always x > y
    tag.R2 <- grep("R2", colnames(df))
    x.R2 <- grep( paste(type, "R2", sep = ".*"), colnames(df) )
    y.R2 <- tag.R2[!(tag.R2 %in% x.R2)]
    message( "Level 1: don't have Enrichment score\n",
             "Level 2: have a ES <= 0\n",
             "Level 3: have a positive ES, but with a corresponding P value > 0.05\n",
             "Level 4: Well predicted (somewhat, depends on the ES" )
    tag.level.1 <- 
    df[ ,x.R2] <- apply(df[ ,x.R2], 1:2, function(x) 
                        ifelse(is.na(x), NA, ifelse(
                        (x > 0) && (x <= 0.05), x, NA)) )
    df[ ,x.ES] <- apply(df[ ,x.ES], 1:2, function(x)
                        ifelse(is.na(x), NA, ifelse(
                        x > 0, x, NA)) )
    
}