name.convert <- function(df, regex='', genelist='',
                         description=T) {
    #' Convert ensembl to gene symbol.
    #'
    #' Require dplyr to work (uses inner_join function).
    #' Require biomaRt to work (retrieve annotation).
    #'  
    #' @param df The data frame whose rownames are to be converted.
    #' @param regex The regular expression for extracting genes of interest.
    #' @param genelist The gene symbols of genes of interest.
    library(dplyr)
    library(biomaRt)
    # Check rownames of df
    if ("ensembl" %in% colnames(df)) {
        rownames(df) <- as.character(df$ensembl)
        df$ensembl <- NULL
    }
    # Retrieve biomaRt annotation data frame
    message("Downloading biomaRt manifest...")
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    ensembl <- sub("(.*)\\..*$", "\\1", rownames(df))
    message("Retrieving data...")
    nameAnnot <- getBM(filters= "ensembl_gene_id",
                       attributes= c("ensembl_gene_id","hgnc_symbol", "description"),
                       values=ensembl, mart= mart)
    # get collagen gene names
    if (regex != '') {
        nameAnnot <- nameAnnot[grep(regex, nameAnnot$hgnc_symbol), ]
    }
    else if (all(genelist != '')) {
        genelist <- paste0(sub("(.*)", "(^\\1)", genelist), collapse="|")
        nameAnnot <- nameAnnot[grep(genelist, nameAnnot$hgnc_symbol), ]
        
    }
    else {
        message("No regex or gene list provided, returning all genes...")
    }
    # nameAnnot <- nameAnnot[nameAnnot$ensembl_gene_id %in% ensembl, ]
    df$ensembl_gene_id <- ensembl
    df <- inner_join(nameAnnot, df, by = "ensembl_gene_id")
    df <- df[order(rowSums(df[ ,4:ncol(df)]), decreasing=T), ]
    df <- df[!duplicated(df$hgnc_symbol), ]
    rownames(df) <- df$hgnc_symbol
    df$hgnc_symbol <- NULL
    if (description) {
        return(df)
    }
    else {
        return(df[ ,3:ncol(df)])
    }
}


edgeR.de.test <- function(df1, df2) {
    #' Use edgeR to perform gene differential expression analysis.
    #'
    #' Require edgeR and dplyr to work.
    #' 
    #' @param df1 First data frame / matrix. Must be counts value!!
    #' @param df2 Second. The result comes as df2:df1.
    library(edgeR)
    df1 <- df1[order(rownames(df1)), ]
    df2 <- df2[order(rownames(df2)), ]
    group <- c(rep("i_iii", ncol(df1)), rep("iv", ncol(df2)))
    x <- cbind(df1, df2)
    y <- DGEList(counts=x, group=group)
    y <- estimateDisp(y)
    ## logFC logCPM PValue
    et <- exactTest(y)$table
    colnames(et) <- c("log2_fold_change", "log2_CPM", "p_value")
    et$symbol <- rownames(et)
    rownames(et) <- NULL
    return(et)
}


transform.for.plot <- function(df, stage='stage_i_iii') {
    #' Melt data frame for ggplot.
    #'
    #' @param df The data frame prepared to melt.
    #' @param stage The tumor stage of `df`.
    library(reshape2)
    df$symbol <- rownames(df)
    df <- melt(df, id.vars = 'symbol', variable.name = 'sample', value.name = 'log2_FPKM')
    df$log2_FPKM <- log2(as.numeric(as.character(df$log2_FPKM)))
    df$stage <- stage
    return(df)
}


get.DEG <- function(df, result, p_value=0.01, logFC=1, logCPM=1,
                    FPKM=0, up=F, strict=F) {
    #' Filter the edgeR results to get only DEGs.
    #'
    #' @param df [data.frame] The FPKM values .
    #' @param result [data.frame] The result given by edgeR.
    #' @param p_value [float] The p-value cutoff.
    #' @param logFC [int] the log2(Fold Change) cutoff.
    #' @param logCPM [int] the log2(Count Per Million) cutoff.
    #' @param FPKM [int] the filtering value for df, useless
    #' if strict=F.
    #' @param up [bool] Whether to keep the down-regulated
    #' genes.
    #' @param strict [bool] Filtering rows in df with too many
    #' zeros.
    #' @return [data.frame] a filtered df.
    if (strict) {
        df <- df[(rowSums(df <= FPKM) <= (ncol(df) / 2)), ]
    }
    result <- result[(abs(result$log2_fold_change) >= logFC) &
                     (result$p_value <= p_value) &
                     (result$log2_CPM >= logCPM), ]
    if (up) {
        result <- result[result$log2_fold_change >= logFC, ]
    }
    df <- df[rownames(df) %in% result$symbol, ]
    df$hgnc_symbol <- rownames(df)
    return(df)
}