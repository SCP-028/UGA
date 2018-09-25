# Use results from DOI:10.1038/nature11412 to get triple-negative samples
# Find differentially expressed genes with edgeR
# Plot FPKM values into a heatmap
try(setwd("C:/Users/jzhou/Desktop/expression_count"), silent=T)
try(setwd("/home/yizhou/data/expression_count"), silent=T)
source("../DEG_func.R")  # edgeR.de.test, retrieve.ensembAnnot
## Get DEGs of the sialic acid & cadherin genes in triple negative BRCA
negAnnot <- data.table::fread("../annotation/triple_negative.csv", data.table = F,
                           stringsAsFactors = F, skip = 1, header = T)
negAnnot <- negAnnot$`Complete TCGA ID`[(negAnnot$`ER Status` == "Negative") &
                                        (negAnnot$`PR Status` == "Negative") &
                                        (negAnnot$`HER2 Final Status` == "Negative")]
load("BRCA.RData")
datan <- datan[ ,sub("(.*)-.*$", "\\1", colnames(datan)) %in% negAnnot]
datat <- datat[ ,sub("(.*)-.*$", "\\1", colnames(datat)) %in% negAnnot]
load("../sialic_acid/Sialic_acids_marker.RData")  # sialic acid marker genes
cadGenes <- c("P?CDH[AB]?\\d{1,2}", "DSG[1234]")  # cadherin marker genes
genelist <- data.frame(symbol=c(sia.sia, sia.gal, sia.galnac,  # sialiac transferases
                                siaP,  # sialic acid related proteins
                                siaEnzy, sialidases,  # sialic acids synthesis & sialidases
                                cadGenes),
                       category=c(rep("transferase_sialic", length(sia.sia)),
                                  rep("transferase_gal", length(sia.gal)),
                                  rep("transferase_galnac", length(sia.galnac)),
                                  rep("sialic_acid_related", length(siaP)),
                                  rep("sialic_acid_synthesis", length(siaEnzy)),
                                  rep("sialidase", length(sialidases)),
                                  rep("cadherin", length(cadGenes))))
genelist <- genelist[!duplicated(genelist$symbol), ]


name.convert <- function(df, ensemblAnnot, regex='', genelist='', description=F) {
    #' Convert ensembl to gene symbol.
    #'
    #' Require dplyr to work (uses inner_join function).
    #'  
    #' @param df The data frame whose rownames are to be converted.
    #' @param ensemblAnnot Annotation data.frame from `retrieve.ensembAnnot`.
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
    return(df[ ,3:5])
}


deg <- edgeR.de.test(datan, datat, "normal", "tumor")
deg <- deg[(deg$p_value <= 0.05) & (abs(deg$log2_fold_change) >= 1), ]
ensemblAnnot <- retrieve.ensembAnnot(deg)
deg <- name.convert(deg, ensemblAnnot, genelist=genelist$symbol)
deg$symbol <- rownames(deg)
deg <- left_join(deg, genelist, by="symbol")
deg$category[is.na(deg$category)] <- "cadherin"

## Boxplot using FPKM values
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsci)
try(setwd("C:/Users/jzhou/Desktop/expression_FPKM"), silent=T)
try(setwd("/home/yizhou/data/expression_FPKM"), silent=T)
load("BRCA.RData")
datan <- datan[ ,sub("(.*)-.*$", "\\1", colnames(datan)) %in% negAnnot]
datat <- datat[ ,sub("(.*)-.*$", "\\1", colnames(datat)) %in% negAnnot]
stageAnnot <- data.table::fread("./annotation/annot.tsv", data.table=F,
                                stringsAsFactors=F, header=T)

name.convert <- function(df, ensemblAnnot, codingOnly=F, regex='', genelist='', description=F, libSize=F) {
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
            return(list(df, library_size))
        }
        else {
            return(list(df[ ,grep("TCGA", colnames(df))], library_size))
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


transform.data <- function(df, annot) {
    #' Separate different stages, and melt for ggplot.
    #'
    #' Require reshape2 to work (uses melt function).
    #'  
    #' @param df The data frame whose rownames are to be converted.
    #' @param annot Stage information.
    #'
    #' @return A melted data.frame with different stages.
    library(reshape2)
    df$stage <- 0
    stage1 <- annot$barcode[grepl("(^i$)|(\\si[abc]?$)|(1)", annot$tumor_stage) &
                                annot$project == "BRCA"]
    stage2 <- annot$barcode[grepl("(^ii$)|(\\si{2}[abc]?$)|(2)", annot$tumor_stage) &
                                annot$project == "BRCA"]
    stage3 <- annot$barcode[grepl("(^iii$)|(\\si{3}[abc]?$)|(3)", annot$tumor_stage) &
                                annot$project == "BRCA"]
    stage4 <- annot$barcode[grepl("(^iv$)|(\\siv[abc]?$)|(4)", annot$tumor_stage) &
                                annot$project == "BRCA"]
    df$stage[df$variable %in% stage1] <- "i"
    df$stage[df$variable %in% stage2] <- "ii"
    df$stage[df$variable %in% stage3] <- "iii"
    df$stage[df$variable %in% stage4] <- "iv"
    df$stage[grepl(".*-1\\d[ABC]$", df$variable)] <- "normal"
    df <- df[df$stage != 0, ]
    return(df)
}


datan <- name.convert(datan, ensemblAnnot)
datat <- name.convert(datat, ensemblAnnot)
datan <- datan[rownames(datan) %in% deg$symbol, ]
datat <- datat[rownames(datat) %in% deg$symbol, ]
datan$symbol <- rownames(datan)
datat$symbol <- rownames(datat)
datan <- melt(datan, id.vars="symbol")
datat <- melt(datat, id.vars="symbol")
datan$stage <- "normal"
datat$stage <- "tumor"
df <- rbind.data.frame(datan, datat)
df <- left_join(df, deg, by="symbol")
df$value <- log2(df$value + 0.1)
df <- transform.data(df, stageAnnot)  # stages
df$stage <- factor(df$stage, levels=c("normal", "i", "ii", "iii", "iv"))
df <- df %>% group_by(symbol) %>% mutate(meanFPKM=mean(value))
df <- df[order(df$p_value), ]
df$symbol <- factor(df$symbol, levels=unique(df$symbol))
plt <- ggboxplot(df, x="stage", y="value", facet.by = "symbol",
                 color="category", notch=F,
                 title="Triple Negative Breast Cancer")
p <- ggpar(plt, legend="right", palette="npg")+
     stat_compare_means(ref.group="normal", label="p.signif", label.y=7)+
     geom_hline(aes(yintercept = meanFPKM), linetype = 2)+
     theme(
         axis.title=element_blank(),
         text=element_text(family="HelveticaNeue"),
         )
ggsave(p, filename="../triple_negative.tiff", device="tiff", units="in",
       width=16, height=9)
