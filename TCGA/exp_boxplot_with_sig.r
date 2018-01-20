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


transform.data <- function(df, project) {
    library(reshape2)
    # df <- df[rowMeans(df) >= 1, ]
    df$symbol <- rownames(df)
    df <- melt(df, id.vars = "symbol")
    df$project <- project
    return(df)
}


exp.boxplot <- function(datan, datat, symbol) {
    library(dplyr)
    library(ggplot2)
    library(ggpubr)

    if (nrow(datan) != 0 & nrow(datat) != 0){
        datan$sample <- "control"
        datat$sample <- "tumor"
        df <- rbind.data.frame(datan, datat)
        colnames(df) <- c("symbol", "barcode", "log2FPKM", "project", "sample")
        df <- df[df$symbol == symbol, ]
        df$log2FPKM <- log2(df$log2FPKM + 0.1)
        # Calculate mean of each group
        df <- df %>% group_by(symbol) %>% mutate(meanFPKM=mean(log2FPKM))
        # Each group is compared to the mean value
        p <- ggplot(df, aes(sample, log2FPKM))+
                geom_boxplot()+
                facet_wrap(facets = "project")+
                stat_compare_means(ref.group="control", label="p.signif")+
                geom_hline(aes(yintercept = meanFPKM), linetype = 2)+
                ggtitle(symbol)
        message(paste0(symbol, " boxplot finished."))
        return(p)
    }
    else {
        message(paste0(symbol, " not found in table."))
    }
}

########## Fill in this part ##########
FILENAMES <- list.files("./expression_FPKM", pattern="RData$", full.names=T)
PROJECTS <- sub("^.*/([^/]*).RData$", "\\1", FILENAMES)
GENELIST <- c("NR1H4", "SOCS1", "SOCS2")
#######################################

dfn <- data.frame()
dft <- data.frame()

for (i in seq_along(PROJECTS)) {
    load(FILENAMES[i])
    if (ncol(datan) >= 10 & ncol(datat) >= 10) {
        datan <- name.convert(datan, annot, genelist=GENELIST)
        datat <- name.convert(datat, annot, genelist=GENELIST)
        datan <- transform.data(datan, PROJECTS[i])
        datat <- transform.data(datat, PROJECTS[i])
        datan <- datan[datan$symbol %in% datat$symbol, ]
        datat <- datat[datat$symbol %in% datan$symbol, ]
        if (nrow(datan) != 0 & nrow(datat) != 0) {
            dfn <- rbind.data.frame(dfn, datan)
            dft <- rbind.data.frame(dft, datat)
        }
    }
}

for (i in seq_along(GENELIST)) {
   p <- exp.boxplot(dfn, dft, GENELIST[i])
   ggsave(filename=paste0(GENELIST[i], ".tiff"),
               plot=p, device="tiff", width=16, height=9, units="in", dpi=200)
}
