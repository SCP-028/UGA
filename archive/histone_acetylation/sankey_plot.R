library(dplyr)
library(networkD3)
setwd("~/data/expression_FPKM")
projects <- list.files("./", pattern = "RData$")
coefs <- read.csv("../genes related to acetate.csv")
colnames(coefs) <- c("symbol", "vmax")
stdf <- read.csv("../mitochondrial.csv")

name.convert <- function(df, ensemblAnnot, codingOnly=F,
                         regex='', genelist='', description=F) {
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
  #'
  #' @return A data.frame df with its rownames converted.
  library(dplyr)
  # Check rownames of df
  if ("ensembl" %in% colnames(df)) {
    rownames(df) <- as.character(df$ensembl)
    df$ensembl <- NULL
  }
  ensembl <- sub("(.*)\\..*$", "\\1", rownames(df))
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
  if (description) {
    return(df)
  }
  else {
    return(df[ ,grep("TCGA", colnames(df))])
  }
}


sankey.plot <- function(df, vmax, stdf) {
  #' Plot a sankey network.
  #' 
  #' @param df The data frame of FPKM values.
  #' @param vmax A data frame containing all Vmax values.
  #' @param stdf Source -> Target.
  library(dplyr)
  df <- rowMeans(df)
  vmax <- vmax[vmax$symbol %in% names(df), ]
  vmax <- vmax[match(names(df), vmax$symbol), ]
  df <- df * vmax$vmax
  df <- as.data.frame(df)
  colnames(df) <- "value"
  df$symbol <- rownames(df)
  
  df <- inner_join(stdf, df, by = "symbol")
  df <- rbind.data.frame(data.frame(source=df$source, target=df$symbol, value=df$value),
                           data.frame(source=df$symbol, target=df$target, value=df$value))
  nodes <- as.data.frame(unique(c(as.character(df$source), as.character(df$target))))
  colnames(nodes) <- "name"
  df$source <- match(df$source, nodes$name) - 1
  df$target <- match(df$target, nodes$name) - 1
  sankeyNetwork(Links = df, Nodes = nodes, Source = "source", Target = "target", Value = "value",
                NodeID = "name", fontSize = 20, nodeWidth = 30)
}


for(i in seq_along(projects)) {
  load(projects[i])
  datan <- name.convert(datan, annot, genelist = coefs$genes.related.to.acetate)
  datat <- name.convert(datat, annot, genelist = coefs$genes.related.to.acetate)
  sankey.plot(datan, coefs, stdf)
  sankey.plot(datat, coefs, stdf)
}