library(ggpubr)

setwd("C:/Users/yz73026/Desktop/")

filenames <- list.files("./expression_FPKM", pattern="RData$", full.names=T)
projects <- sub("^.*/([^/]*).RData$", "\\1", filenames)
genelist <- c("SLC7A1", "ODC1", "SLC3A2", "GATM", "PYCR1", "PYCR2", "PYCR3")

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


prep.cor.heatmap <- function(df, genelist, cluster=F,
                             sample=c("normal", "tumor")) {
  #' Calculate spearman correlation coefficients and cluster. 
  #'
  #' @param df The input values.
  #' @param genelist Used for arranging the order of the genes.
  #' @param cluster Whether to cluster the genes according to coefficients.
  #' @param sample Normal or tumor samples, used for later plotting.
  #'
  #' @return A melted data.frame of the correlation coefficients.
  library(reshape2)
  df <- df[match(genelist, rownames(df)), ]
  df <- df[apply(df, 1, var) != 0, ]
  df <- cor(t(df), method="spearman")
  if (cluster) {
    dd <- as.dist((1 - df) / 2)
    hc <- hclust(dd)
    df <- df[hc$order, hc$order]    
  }
  df[upper.tri(df)] <- NA
  df <- melt(df, na.rm=T)
  df$sample <- sample
  colnames(df) <- c("Var1", "Var2", "value", "sample")
  return(df)
}


for (i in seq_along(projects)) {
    load(filenames[i])
    if (ncol(datan) >= 10 & ncol(datat) >= 10) {
      datan <- name.convert(datan, annot)
      datat <- name.convert(datat, annot)
      dfn <- prep.cor.heatmap(datan, genelist, sample = "control")
      dft <- prep.cor.heatmap(datat, genelist, sample = "tumor")
      pn <- nice.heatmap(dfn, paste0(projects[i], "_control"))
      pt <- nice.heatmap(dft, paste0(projects[i], "_tumor"))
      p <- ggarrange(pn, pt, ncol=2)
      ggsave(filename=paste0(projects[i], ".tiff"), plot=p, device="tiff", width=16, height=9, units="in")
    }
}
