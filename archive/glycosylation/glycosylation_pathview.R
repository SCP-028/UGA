library(edgeR)
library(dplyr)
library(ggplot2)

setwd("C:/Users/yz73026/Desktop/")

network <- read.csv("./glycosylation/glycosylation_pathway.csv")
colnames(network) <- c("source", "interaction", "target", "description")
genelist <- unique(c(as.character(network$source[network$interaction == "rc"]),
                     as.character(network$target[network$interaction == "cr"]),  # enzymes in reactions
                     as.character(network$target[network$interaction == "pp"]),  # protein-protein interaction
                     as.character(network$source[network$interaction == "rp"])))  # glycosylation pool enzymes
compound <- unique(c(as.character(network$source[!network$source %in% genelist]),
                     as.character(network$target[!network$target %in% genelist])))

projects <- c("BRCA", "COAD", "LUAD", "LUSC", "BLCA", "LIHC", "HNSC", "KICH", "THCA", "PRAD")
filenames <- paste0("./expression_count/", projects, ".RData")


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
      return(list(exp=df[ ,grep("TCGA", colnames(df))],
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


edgeR.de.test <- function(df1, df2, group1, group2, sepLibSize=F) {
  #' Use edgeR to perform gene differential expression analysis.
  #'
  #' Require edgeR and dplyr to work.
  #'
  #' @param df1 First data frame / matrix. Must be counts value!!
  #' @param df2 Second. The result comes as df2:df1.
  #' @param group1 Name of first group.
  #' @param group2 Name of second group.
  #' @param sepLibSize If `name.convert` has a separate library size,
  #'                   set this to TRUE.
  library(edgeR)
  if(sepLibSize) {
    libSize <- c(df1$libSize, df2$libSize)
    df1 <- df1$exp[order(rownames(df1$exp)), ]
    df2 <- df2$exp[order(rownames(df2$exp)), ]
    group <- c(rep(group1, ncol(df1)), rep(group2, ncol(df2)))
    design <- model.matrix(~0+group)
    x <- cbind(df1, df2)
    y <- DGEList(counts=x, group=group, lib.size=libSize)
  }
  else {
    df1 <- df1[order(rownames(df1)), ]
    df2 <- df2[order(rownames(df2)), ]
    group <- c(rep(group1, ncol(df1)), rep(group2, ncol(df2)))
    design <- model.matrix(~0+group)
    x <- cbind(df1, df2)
    y <- DGEList(counts=x, group=group)
  }
  y <- y[rowSums(cpm(y) > 1) >= 2, , keep.lib.sizes=F]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design=design)
  ## logFC logCPM PValue
  et <- exactTest(y)$table
  colnames(et) <- c("log2_fold_change", "log2_CPM", "p_value")
  et$symbol <- rownames(et)
  rownames(et) <- NULL
  return(et)
}


result <- data.frame()
unlisted <- data.frame(symbol=genelist)

for (i in seq_along(projects)) {
  load(filenames[i])
  datan <- name.convert(datan, annot, genelist=genelist, libSize = T)
  datat <- name.convert(datat, annot, genelist=genelist, libSize = T)
  df <- edgeR.de.test(datan, datat, "control", "tumor", sepLibSize = T)
  # Some unfound genes
  df <- left_join(unlisted, df, by="symbol")
  df[is.na(df)] <- 0
  df$project <- projects[i]
  result <- rbind.data.frame(result, df)
}

nodes <- data.frame(symbol=c(genelist, compound), type=c(rep("enzyme", length(genelist)), rep("compound", length(compound))))
nodes <- left_join(nodes, result, by="symbol")
nodes[is.na(nodes)] <- 0

for (i in seq_along(projects)) {
  df <- nodes[nodes$type == "compound" | nodes$project == projects[i], ]
  write.csv(df, file=paste0("./glycosylation/nodes_", projects[i], ".sif"), quote=F, row.names = F)
}

write.csv(nodes, file="./glycosylation/cytoscape_nodes.sif", quote = F, row.names = F)
