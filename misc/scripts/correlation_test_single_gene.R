library(ggplot2)

the_gene <- c("BCAT1", "BCAT2")
setwd("C:/Users/yz73026/Desktop/")
load("./differentially_expressed_genes.RData")  # result
deg <- result[(abs(result$log2_fold_change) >= 1) & (result$p_value <= 0.001) & (abs(result$log2_CPM) >= 1), ]

filenames <- list.files("./expression_FPKM", pattern="RData$", full.names=T)
projects <- sub("^.*/([^/]*).RData$", "\\1", filenames)

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


nice.heatmap <- function(df, title) {
    #' Plot a nice-looking heatmap.
    #'
    #' @param df The melted data.frame used by ggplot.
    #' @param title The title of the plot.
    #'
    #' @return A ggplot list.
    library(ggplot2)
    ggplt <- ggplot(df, aes(Var1, Var2, fill=value))+
        geom_tile(color="white")+
        ggtitle(title)+
        scale_fill_gradient2(low="#6D9EC1", high="#E46726", mid="white",
                             midpoint=0, limit=c(-1, 1), space="Lab",
                             name="Spearman\nCorrelation")+
        coord_fixed()+
        theme_minimal()+
        theme(
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            panel.border=element_blank(),
            text=element_text(family="HelveticaNeue"),
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
            legend.justification=c(0, 1),
            legend.position=c(0, 1),
            legend.direction="horizontal")+
        guides(fill=guide_colorbar(barwidth=7, barheight=1,
                                   title.position="top",
                                   title.hjust=0.5))
    return(ggplt)
}


prep.cor.heatmap <- function(gene, genelist, df, cutoff=0.7, sample="") {
  #' Calculate spearman correlation coefficients and cluster.
  #'
  #' @param gene The gene you want to test.
  #' @param genelist The list of genes you want to test against.
  #' @param df The whole gene expression dataset.
  #' @param cutoff Calculated spearman correlation cutoff.
  #' @param sample Normal or tumor samples, used for later plotting.
  #'
  #' @return A melted data.frame of the correlation coefficients.
  library(reshape2)
  gene <- t(df[rownames(df) == gene, ])
  genelist <- t(df[rownames(df) %in% genelist, ])
  df <- cor(gene, genelist, method="spearman")
  df <- melt(df, na.rm=T)
  df$sample <- sample
  colnames(df) <- c("Var1", "Var2", "value", "sample")
  df <- df[df$value >= cutoff, ]
  return(df)
}


result <- data.frame()
for (i in seq_along(projects)) {
    load(filenames[i])
    if (ncol(datan) >= 10 & ncol(datat) >= 10) {
        genelist <- deg$symbol[deg$project == projects[i]]
        datan <- name.convert(datan, annot)
        datat <- name.convert(datat, annot)
        for (j in seq_along(the_gene)) {
            dfn <- prep.cor.heatmap(gene=the_gene[j], genelist=genelist, df=datan,
                                    sample=paste0(projects[i], "_control_", the_gene[j]))
            dft <- prep.cor.heatmap(gene=the_gene[j], genelist=genelist, df=datat,
                                    sample=paste0(projects[i], "_tumor_", the_gene[j]))
            result <- rbind.data.frame(result, dfn)
            result <- rbind.data.frame(result, dft)
        }
    }
}

if (cluster) {
dd <- as.dist((1 - df) / 2)
hc <- hclust(dd)
df <- df[hc$order, hc$order]
}