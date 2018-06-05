rm(list=ls())
library(glue)
library(tidyverse)
library(foreach)
library(doSNOW)
# source("https://bioconductor.org/biocLite.R")
# biocLite("reactome.db")
library(fgsea)
source("./genelist.R")

setwd("~/data")
projects <- c(
    "BLCA", "BRCA", "COAD", "HNSC", "KICH", "KIRC", "LUAD",
    "LIHC","LUSC", "PRAD", "READ", "SKCM", "STAD", "THCA"
)
filenames <- glue("./expression_FPKM/{projects}.RData")
genelists <- load_aa_data()
pathway_names <- names(genelists)

cl <- makeCluster(length(projects))
registerDoSNOW(cl)


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
    require(dplyr)
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
    df <- df[order(rowSums(df[ ,grep("(TCGA)|(TARGET)", colnames(df))]), decreasing=T), ]
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
            return(df[ ,grep("(TCGA)|(TARGET)", colnames(df))])
        }
    }
}


result <- foreach (i = seq_along(projects)) %dopar% {
    load(filenames[i])
    datat <- as.data.frame(t(name.convert(datat, annot)))
    df <- data.frame()
    for (j in seq_along(genelists)) {
        genelist <- genelists[[j]]
        for (k in seq_along(genelist)) {
            tmp <- t(cor(
                as.matrix(select(datat, genelist[k])),
                as.matrix(select(datat, -one_of(genelist[k]))),
                method="pearson"))
            tmp <- as.data.frame(tmp)
            colnames(tmp) <- "value"
            tmp$gene <- rownames(tmp)
            tmp$enzyme <- genelist[k]
            tmp$pathway <- names(genelists)[j]
            df <- rbind.data.frame(df, tmp)
        }
    }
    df %>%
      filter(!is.na(value)) %>%
      arrange(pathway, enzyme, desc(abs(value))) %>%
      select(pathway, enzyme, gene, value)
}
names(result) <- projects
save(result, file = "aa_paper_correlation.RData")

df <- data.frame()
for (i in seq_along(result)) {
  tmp <- result[[i]] %>%
    group_by(pathway, enzyme) %>%
    top_n(n=20, wt=abs(value)) %>%
    mutate(cancer=names(result)[i]) %>%
    select(cancer, pathway, enzyme, gene, value)
  df <- rbind.data.frame(df, tmp)
}
write_excel_csv(df, path="./aa_paper_top_correlation.csv")

all_genes <- as.character(annot$entrezgene[!is.na(annot$entrezgene)])
pathways <- reactomePathways(unique(all_genes))
for (i in seq_along(projects)) {
    for (j in seq_along(pathway_names)) {
        tmp <- df[df$cancer == projects[i] & df$pathway == pathway_names[j], ]
        cor_genes <- unique(tmp$gene)
        cor_annot <- annot[annot$hgnc_symbol %in% cor_genes, ] %>%
        select(hgnc_symbol, entrezgene) %>%
        rename(gene = hgnc_symbol) %>%
        left_join(tmp, by="gene") %>%
        select(entrezgene, value) %>%
        arrange(abs(value)) %>%
        distinct(entrezgene, .keep_all = T)
        ranks <- cor_annot$value
        names(ranks) <- as.character(cor_annot$entrezgene)
        gseaRes <- fgsea(pathways, ranks, minSize=10, maxSize=500,nperm=10000)
        gseaRes <- gseaRes %>%
        filter(pval <= 0.05) %>%
        arrange(desc(ES))
        data.table::fwrite(gseaRes, file=glue("{projects[i]}_{pathway_names[j]}_gsea.tsv"), sep="\t", sep2=c("", " ", ""))
    }
}

# sessionInfo()
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux 9 (stretch)

# Matrix products: default
# BLAS/LAPACK: /usr/lib/libopenblasp-r0.2.19.so

# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C              LC_PAPER=en_US.UTF-8       LC_NAME=C
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base

# other attached packages:
#  [1] fgsea_1.4.1     Rcpp_0.12.16    doSNOW_1.0.16   snow_0.4-2      iterators_1.0.9 foreach_1.4.4   forcats_0.3.0
#  [8] stringr_1.3.0   purrr_0.2.4     readr_1.1.1     tidyr_0.8.0     tibble_1.4.2    ggplot2_2.2.1   tidyverse_1.2.1
# [15] dplyr_0.7.4     glue_1.2.0

# loaded via a namespace (and not attached):
#  [1] reshape2_1.4.3      haven_1.1.1         lattice_0.20-35     colorspace_1.3-2    rlang_0.2.0         pillar_1.2.2
#  [7] foreign_0.8-70      BiocParallel_1.12.0 modelr_0.1.1        readxl_1.1.0        bindrcpp_0.2.2      bindr_0.1.1
# [13] plyr_1.8.4          munsell_0.4.3       gtable_0.2.0        cellranger_1.1.0    rvest_0.3.2         codetools_0.2-15
# [19] psych_1.8.3.3       parallel_3.4.1      broom_0.4.4         scales_0.5.0        jsonlite_1.5        gridExtra_2.3
# [25] fastmatch_1.1-0     mnormt_1.5-5        hms_0.4.2           stringi_1.2.2       grid_3.4.1          cli_1.0.0
# [31] tools_3.4.1         magrittr_1.5        lazyeval_0.2.1      crayon_1.3.4        pkgconfig_2.0.1     data.table_1.11.0
# [37] xml2_1.2.0          lubridate_1.7.4     assertthat_0.2.0    httr_1.3.1          rstudioapi_0.7      R6_2.2.2
# [43] nlme_3.1-137        compiler_3.4.1
