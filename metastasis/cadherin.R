try(setwd("C:/Users/jzhou/Desktop/expression_FPKM"), silent=T)
try(setwd("/home/yizhou/data/expression_FPKM"), silent=T)

annot <- read.csv("./annotation/annot.tsv")
filen <- list.files("./", pattern = "normal.csv$")
filet <- list.files("./", pattern = "tumor.csv$")
project <- intersect(sub("(.*)_.*$", "\\1", filen),
                     sub("(.*)_.*$", "\\1", filet))
rm(filen, filet)


name.convert <- function(df, ensemblAnnot='', codingOnly=F,
                         regex='', genelist='', description=F) {
    #' Convert ensembl to gene symbol.
    #'
    #' Require dplyr to work (uses inner_join function).
    #' Require biomaRt to work (retrieve annotation).
    #'  
    #' @param df The data frame whose rownames are to be converted.
    #' @param ensemblAnnot Saved annotation data.frame, can be used to save time.
    #' @param codingOnly Keep only the protein-coding genes and ditch the rest.
    #' @param regex The regular expression for extracting genes of interest.
    #' @param genelist The gene symbols of genes of interest.
    #'
    #' @return A data.frame df with its rownames converted.
    library(dplyr)
    library(biomaRt)
    # Check rownames of df
    if ("ensembl" %in% colnames(df)) {
        rownames(df) <- as.character(df$ensembl)
        df$ensembl <- NULL
    }
    # Retrieve biomaRt annotation data frame
    ensembl <- sub("(.*)\\..*$", "\\1", rownames(df))
    if (all(ensemblAnnot == '')) {
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
    }
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
    # ensemblAnnot <- ensemblAnnot[ensemblAnnot$ensembl_gene_id %in% ensembl, ]
    df$ensembl_gene_id <- ensembl
    df <- inner_join(ensemblAnnot, df, by = "ensembl_gene_id")
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


transform_data <- function(df, annot, ensAnnot) {
    library(reshape2)
    df <- name.convert(df, ensemblAnnot = ensAnnot, regex = "(^P?CDH[AB]?\\d{1,2}$)|(^DSG[1234]$)")
    df <- df[rowMeans(df) >= 1, ]
    df$symbol <- rownames(df)
    df <- melt(df, id.vars = "symbol")
    df$stage <- 0
    stage1 <- annot$barcode[grepl("(^i$)|(\\si[abc]?$)|(1)", annot$tumor_stage) &
                                annot$project == project[i]]
    stage2 <- annot$barcode[grepl("(^ii$)|(\\si{2}[abc]?$)|(2)", annot$tumor_stage) &
                                annot$project == project[i]]
    stage3 <- annot$barcode[grepl("(^iii$)|(\\si{3}[abc]?$)|(3)", annot$tumor_stage) &
                                annot$project == project[i]]
    stage4 <- annot$barcode[grepl("(^iv$)|(\\siv[abc]?$)|(4)", annot$tumor_stage) &
                                annot$project == project[i]]
    df$stage[df$variable %in% stage1] <- "i"
    df$stage[df$variable %in% stage2] <- "ii"
    df$stage[df$variable %in% stage3] <- "iii"
    df$stage[df$variable %in% stage4] <- "iv"
    df <- df[df$stage != 0, ]
    return(df)
}

library(ggplot2)
library(ggpubr)

for (i in seq_along(project)) {
    datan <- data.table::fread(paste0(project[i], "_normal.csv"), header = T,
                               stringsAsFactors = F, sep = "\t", data.table = F)
    datat <- data.table::fread(paste0(project[i], "_tumor.csv"), header = T,
                               stringsAsFactors = F, sep = "\t", data.table = F)
    if (ncol(datan) <= 5 | ncol(datat) <= 5) {
        message(paste0(project[i], " ignored because it has fewer than 5 cases."))
    }
    else {
        datan <- transform_data(datan, annot, ensAnnot)
        datat <- transform_data(datat, annot, ensAnnot)
        datan <- datan[datan$symbol %in% datat$symbol, ]

        if (nrow(datan) != 0 & nrow(datat) != 0){
            datan$sample <- "control"
            datat$sample <- "tumor"
            df <- rbind.data.frame(datan, datat)
            colnames(df) <- c("symbol", "barcode", "log2FPKM", "stage", "sample")
            df$stage[df$sample == "control"] <- "control"
            df$log2FPKM <- log2(df$log2FPKM + 0.1)
            df$symbol <- factor(df$symbol, levels=(unique(df$symbol[order(df$symbol,
                                                                          nchar(df$symbol))])))
            genelist <- levels(df$symbol)
            # Remove genes that ANOVA says there's no difference between groups
            for (j in seq_along(genelist)) {
                temp <- df[df$symbol == genelist[j], ]
                temp <- summary(aov(log2FPKM ~ stage, data=temp))
                if (temp[[1]]$`Pr(>F)`[1] > 0.01) {
                    df <- df[df$symbol != genelist[j], ]
                }
            }
            # symbol variable value stage sample
            # Calculate mean of each group
            df <- df %>% group_by(symbol) %>% mutate(meanFPKM=mean(log2FPKM))
            # Each group is compared to the mean value
            ggplot(df, aes(stage, log2FPKM))+
                geom_boxplot()+
                facet_wrap(facets = "symbol", switch="both")+
                stat_compare_means(ref.group=".all.", label="p.signif", method="t.test")+
                geom_hline(aes(yintercept = meanFPKM), linetype = 2)+
                ggtitle(project[i])
            ggsave(paste0("../CDH/", project[i], ".tiff"), device = "tiff", width = 16,
                   height = 9, units = "in")
            message(paste0(project[i], " finished."))
        }
        else {
            message(paste0(project[i], " ignored because no CDH gene was found."))
        }
        
    }
}
