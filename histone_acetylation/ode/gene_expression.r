setwd("C:/Users/jzhou/Desktop/expression_FPKM")
filename <- list.files("./", pattern="RData$")
projects <- sub("\\.RData", "", filename)
mitochondria <- c("PDHA1", "CS", "ACOT13", "ACO2",
                  "ACSS1", "PC", "OGDH", "IDH2")
cytosol <- c("ACOT12", "ACLY", "HMGCS1", "HDAC1",
             "HDAC2", "HDAC3", "ACSS2", "SLC16A3",
             "KAT2A", "ACACA")
genelist <- c(mitochondria, cytosol)
genelist <- genelist[order(genelist)]

retrieve.ensembAnnot <- function(df, codingOnly=F) {
    #' Prepare for converting ensembl to gene symbol.
    #'
    #' Require biomaRt to work (retrieve annotation).
    #'  
    #' @param df The data frame whose rownames are to be converted.
    #' @param codingOnly Keep only the protein-coding genes and ditch the rest.
    #'
    #' @return A data.frame of matching ensembl IDs and symbols.
    library(biomaRt)
    
    if ("ensembl" %in% colnames(df)) {
        rownames(df) <- as.character(df$ensembl)
    }
    # Retrieve biomaRt annotation data frame
    ensembl <- sub("(.*)\\..*$", "\\1", rownames(df))
    message("Downloading biomaRt manifest...")
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    message("Retrieving annotation data...")
    ensemblAnnot <- getBM(filters= "ensembl_gene_id",
                          attributes= c("ensembl_gene_id", "hgnc_symbol",
                                        "gene_biotype", "description"),
                          values=ensembl, mart= mart)
    ensemblAnnot <- ensemblAnnot[ensemblAnnot$hgnc_symbol != '', ]
    if (codingOnly) {
        ensemblAnnot <- ensemblAnnot[ensemblAnnot$gene_biotype == 'protein_coding', ]
    }
    ensemblAnnot$gene_biotype <- NULL
    return(ensemblAnnot)
}


name.convert <- function(df, ensemblAnnot, codingOnly=F,
                         regex='', genelist='', description=F,
                         libSize=F) {
    #' Convert ensembl to gene symbol.
    #'
    #' Require dplyr to work (uses inner_join function).
    #'  
    #' @param df The data frame whose rownames are to be converted.
    #' @param ensemblAnnot Annotation data.frame from `retrieve.ensembAnnot`.
    #' @param codingOnly Keep only the protein-coding genes and ditch the rest.
    #' @param regex The regular expression for extracting genes of interest.
    #' @param genelist The gene symbols of genes of interest.
    #' @param description Whether to keep description & ensembl of genes.
    #' @param libSize Keep library size for DEG test or not.
    #'
    #' @return A data.frame df with its rownames converted.
    library(dplyr)
    # Check rownames of df
    if ("ensembl" %in% colnames(df)) {
        rownames(df) <- as.character(df$ensembl)
        df$ensembl <- NULL
    }
    ensembl <- sub("(.*)\\..*$", "\\1", rownames(df))
    library_size <- colSums(df)
    # get collagen gene names
    if (regex != '') {
        message("Regex found, filtering...")
        ensemblAnnot <- ensemblAnnot[grep(regex, ensemblAnnot$hgnc_symbol), ]
    }
    else if (all(genelist != '')) {
        message("Gene list found, filtering...")
        genelist <- paste0(sub("(.*)", "(^\\1$)", genelist), collapse="|")
        ensemblAnnot <- ensemblAnnot[grep(genelist, ensemblAnnot$hgnc_symbol), ]
    }
    else {
        message("No regex or gene list provided, returning all genes...")
    }
    df$ensembl_gene_id <- ensembl
    df <- inner_join(ensemblAnnot, df, by = "ensembl_gene_id")
    df <- df[order(rowSums(df[ ,grep("TCGA", colnames(df))]), decreasing=T), ]
    df <- df[!duplicated(df$hgnc_symbol), ]
    rownames(df) <- df$hgnc_symbol
    df$hgnc_symbol <- NULL
    if (libSize) {
        if (description) {
            return(list(exp=df, libSize=library_size))
        }
        else {
            return(list(exp=df[ ,grep("(TCGA)|(TARGET)", colnames(df))],
                        libSize=library_size))
        }
    }
    else {
        if (description) {
            return(df)
        }
        else {
            return(df[ ,grep("TCGA", colnames(df))])
        }
    }
}

load(filename[1])
ensembl <- retrieve.ensembAnnot(datan)
result <- data.frame(row.names = genelist)
for (i in seq_along(projects)) {
    load(filename[i])
    if (ncol(datan) >= 5 & ncol(datat) >= 5) {
        datan <- name.convert(datan, ensembl, genelist=genelist)
        datat <- name.convert(datat, ensembl, genelist=genelist)
        if (nrow(datan) == length(genelist) & nrow(datat) == length(genelist)) {
            dfn <- rowMeans(datan)
            dfn <- dfn[order(names(dfn))]
            dft <- rowMeans(datat)
            dft <- dft[order(names(dft))]
            result[paste(projects[i], "normal", sep="_")] <- dfn
            result[paste(projects[i], "tumor", sep="_")] <- dft
        }
        else {
            message(
                paste(projects[i], "Missing",
                      length(genelist) - min(nrow(datan), nrow(datat)),
                      "genes!")
            )
        }
    }
}
