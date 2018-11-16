# Lysine degradation pathway.
library(glue)
setwd("C:/Users/yz73026/Desktop")
GENELIST <- c(
    "AASS",  # Aminoadipate-Semialdehyde Synthase
    "ALDH7A1",  # Alpha-AASA Dehydrogenase
    "AADAT",  # Aminoadipate Aminotransferase
    "OGDH",  # Oxoglutarate Dehydrogenase
    "OGDHL",  # Oxoglutarate Dehydrogenase Like
    "DLST",  # Dihydrolipoamide S-Succinyltransferase
    "GCDH",  # Glutaryl-CoA Dehydrogenase
    "ECHS1",  # Enoyl-CoA Hydratase, Short Chain 1
    "EHHADH",  # Enoyl-CoA Hydratase And 3-Hydroxyacyl CoA Dehydrogenase
    "HADHA",  # Hydroxyacyl-CoA Dehydrogenase (Trifunctional Protein), Alpha Subunit
    "HADH",  # Hydroxyacyl-CoA Dehydrogenase
    "ACAT1",  # Acetyl-CoA Acetyltransferase 1
    "ACAT2"  # Acetyl-CoA Acetyltransferase 2
)

stageAnnot <- read.csv("./expression_FPKM/annotation/annot.csv")


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
    df <- df[order(df$hgnc_symbol), ]
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


transform.data <- function(df, stageAnnot, project) {
    #' Separate different stages, and melt for ggplot.
    #'
    #' Require reshape2 to work (uses melt function).
    #'
    #' @param df The data frame whose rownames are to be converted.
    #' @param stageAnnot Stage information.
    #' @param project The TCGA project name.
    #'
    #' @return A melted data.frame with different stages.
    library(reshape2)
    # df <- df[rowMeans(df) >= 1, ]
    df$symbol <- rownames(df)
    df <- melt(df, id.vars = "symbol")
    df$stage <- 0
    stage1 <- stageAnnot$barcode[grepl("(^i$)|(\\si[abc]?$)|(1)", stageAnnot$tumor_stage) &
                                stageAnnot$project == project]
    stage2 <- stageAnnot$barcode[grepl("(^ii$)|(\\si{2}[abc]?$)|(2)", stageAnnot$tumor_stage) &
                                stageAnnot$project == project]
    stage3 <- stageAnnot$barcode[grepl("(^iii$)|(\\si{3}[abc]?$)|(3)", stageAnnot$tumor_stage) &
                                stageAnnot$project == project]
    stage4 <- stageAnnot$barcode[grepl("(^iv$)|(\\siv[abc]?$)|(4)", stageAnnot$tumor_stage) &
                                stageAnnot$project == project]
    df$stage[df$variable %in% stage1] <- "i"
    df$stage[df$variable %in% stage2] <- "ii"
    df$stage[df$variable %in% stage3] <- "iii"
    df$stage[df$variable %in% stage4] <- "iv"
    df <- df[df$stage != 0, ]
    return(df)
}


merge.data <- function(datan, datat, proj, stageAnnot) {
    #' Convert ensembl to gene symbol.
    #'
    #' Require dplyr to work (uses inner_join function).
    #' Require ggplot2 & ggpubr to finish the boxplots.
    #'
    #' @param datan The data frame for normal samples.
    #' @param datat The data frame for tumor samples.
    #' @param proj The current cancer project.
    #' @param stageAnnot Stage information.
    #'
    #' @return Genes used in the boxplot.

    library(dplyr)
    library(ggplot2)
    library(ggpubr)

    datan <- transform.data(datan, stageAnnot, proj)
    datat <- transform.data(datat, stageAnnot, proj)
    datan <- datan[datan$symbol %in% datat$symbol, ]
    datat <- datat[datat$symbol %in% datan$symbol, ]
    if (nrow(datan) != 0 & nrow(datat) != 0) {
        datan$sample <- "control"
        datat$sample <- "tumor"
        df <- rbind.data.frame(datan, datat)
        colnames(df) <- c("symbol", "barcode", "FPKM", "stage", "sample")
        df$stage[df$sample == "control"] <- "control"
        # df$log2FPKM <- log2(df$log2FPKM + 0.1)
        df$symbol <- factor(df$symbol, levels=(unique(df$symbol[order(nchar(df$symbol),
                                                                        df$symbol)])))
        df$project <- proj
        return (df)
    }
    else {
        return(NULL)
    }

}


filenames <- list.files("./expression_FPKM/", pattern=".RData$", full.names = T)
projects <- sub(".*/([A-Z-]+)\\.RData", "\\1", filenames)
result <- data.frame()

for (i in seq_along(projects)) {
    load(filenames[i])
    if (ncol(datan) >= 10 & ncol(datat) >= 10) {
        datan <- name.convert(datan, annot, genelist=GENELIST)
        datat <- name.convert(datat, annot, genelist=GENELIST)
        df <- merge.data(datan, datat, projects[i], stageAnnot)
        result <- rbind.data.frame(result, df)
    }
}

# Modify range of y axis
result$FPKM <- log2(result$FPKM + 0.1)
# Calculate mean of each group
df <- result %>%
    group_by(project, stage, symbol) %>%
    mutate(meanFPKM=mean(FPKM))
projects <- projects[projects %in% df$project]

setwd("./lysine")
for (i in seq_along(projects)) {
    temp <- df[df$project == projects[i], ]
    temp$stage <- factor(temp$stage, levels=c("control", "i", "ii", "iii", "iv"))
    temp$symbol <- factor(temp$symbol, levels=GENELIST)
    temp$meanFPKM[temp$stage != "control"] <- NA
    # Violin plot
    p <- facet(ggviolin(
        temp, x="stage", y="FPKM", fill="stage",
        add="boxplot", add.params = list(color="black", fill = "white")
         ), facet.by="symbol", nrow=2, ncol=7)
    # Add t-test significance
    p <- p+
         stat_compare_means(
             comparisons=list(c("i", "control"), c("ii", "control"),
                              c("iii", "control"), c("iv", "control")),
             label="p.signif"
         )
    # Add horizontal line of the mean of control group
    p <- p+geom_hline(aes(yintercept = meanFPKM), linetype = 2)
    # Aesthetic settings
    p <- ggpar(
        p, xlab="Stages", ylab="log2(FPKM + 0.1)", title=projects[i],
        legend="right", legend.title="", palette="npg", ggtheme=theme_light())
    ggsave(filename=glue("{projects[i]}.tiff"),
           plot=p, device="tiff", width=28, height=8, units="in", dpi=300)
}

# sessionInfo()
# R version 3.4.3 (2017-11-30)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows >= 8 x64 (build 9200)

# Matrix products: default

# locale:
# [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252
# [3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C
# [5] LC_TIME=English_United States.1252

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base

# other attached packages:
# [1] reshape2_1.4.3       ggpubr_0.1.6         magrittr_1.5         ggplot2_2.2.1
# [5] dplyr_0.7.4          glue_1.2.0           RevoUtils_10.0.7     RevoUtilsMath_10.0.1

# loaded via a namespace (and not attached):
#  [1] Rcpp_0.12.15     bindr_0.1        munsell_0.4.3    colorspace_1.3-2 R6_2.2.2
#  [6] rlang_0.2.0      stringr_1.3.0    plyr_1.8.4       tools_3.4.3      grid_3.4.3
# [11] gtable_0.2.0     digest_0.6.15    lazyeval_0.2.1   assertthat_0.2.0 tibble_1.4.2
# [16] ggsignif_0.4.0   bindrcpp_0.2     ggsci_2.8        purrr_0.2.4      labeling_0.3
# [21] stringi_1.1.6    compiler_3.4.3   pillar_1.1.0     scales_0.5.0     pkgconfig_2.0.1
