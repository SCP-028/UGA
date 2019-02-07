library(DESeq2)

library(ggplot2)
library(EnhancedVolcano)
library(BiocParallel)
register(SnowParam((10)))

TOP_N_IN_VOLCANO_ADJP <- 30
VOLCANO_PADJ_CUTOFF <- 1e-15
TOP_N_IN_VOLCANO_FC <- 15

setwd("~/data/TCGA")
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
finished_proj <- unlist(lapply(list.files("./DEA/dds"),
                               function(x) strsplit(x, "\\.")[[1]][[1]]))
DEA_projects <- projects[!projects %in% finished_proj]

# Perform DESeq2 on input projects with missing RData files
if (length(DEA_projects) > 0)
{
  annot <- data.table::fread("./annotation/counts_annotation.csv")
  annot_normal <- annot[annot$sample_type == "Solid Tissue Normal", ]
  annot_tumor <- annot[!annot$fileID %in% annot_normal$fileID, ]

  transform.data <- function(df, normalAnnot, stageAnnot, project) {
    #' Separate different stages of a data.frame based on the column names.
    #'
    #' @param df The data frame with symbols as rownames, and TCGA barcode as column names.
    #' @param normalAnnot Annotation of normal samples.
    #' @param stageAnnot Stage information of tumor samples.
    #' @param project The TCGA project name.
    #'
    #' @return A list of 5 character lists: ("N", "I", "II", "III", "IV").
    barcode_N <- normalAnnot$barcode[normalAnnot$project == project]
    barcode_I <- stageAnnot$barcode[grepl("(^i$)|(\\si[abc]?$)|(1)", stageAnnot$tumor_stage, ignore.case = T) &
                                      stageAnnot$project == project]
    barcode_II <- stageAnnot$barcode[grepl("(^ii$)|(\\si{2}[abc]?$)|(2)", stageAnnot$tumor_stage, ignore.case = T) &
                                       stageAnnot$project == project]
    barcode_III <- stageAnnot$barcode[grepl("(^iii$)|(\\si{3}[abc]?$)|(3)", stageAnnot$tumor_stage, ignore.case = T) &
                                        stageAnnot$project == project]
    barcode_IV <- stageAnnot$barcode[grepl("(^iv$)|(\\siv[abc]?$)|(4)", stageAnnot$tumor_stage, ignore.case = T) &
                                       stageAnnot$project == project]
    return(list(N=barcode_N,
                I=barcode_I,
                II=barcode_II,
                III=barcode_III,
                IV=barcode_IV
    ))
  }

  for (i in seq_along(DEA_projects))
  {
    print(paste0("Working on project ", DEA_projects[i]))
    df <- data.table::fread(paste0("./counts/", DEA_projects[i], ".counts.csv"), data.table = F, stringsAsFactors = F)
    if ("Ensembl" %in% colnames(df)) {
      rownames(df) <- df$Ensembl
      df$Ensembl <- NULL
    }
    df <- head(df, -5)  # remove  __no_feature, __ambiguous,  __too_low_aQual, __not_aligned, and __alignment_not_unique
    df_stages <- transform.data(df, annot_normal, annot_tumor, DEA_projects[i])  # list(N, I, II, III, IV), each item is a list of barcodes
    sample_num <- lapply(df_stages, length)  # number of samples in each group
    samples <- unlist(df_stages)
    df <- df[ ,samples]
    samp_cond <- vector("character")
    for (j in seq_along(sample_num)) {
      samp_cond <- c(samp_cond, rep(names(sample_num)[j], sample_num[j]))
    }
    coldata <- matrix(samp_cond, dimnames = list(samples, "condition"))

    dds <- DESeqDataSetFromMatrix(countData = df,
                                  colData = coldata,
                                  design = ~ condition)
    dds$condition <- relevel(dds$condition, ref = samp_cond[1])
    dds <- dds[rowSums(counts(dds)) >= 10,]
    print(paste0("Starting DESeq2 for ", DEA_projects[i],
                 " with ", nrow(counts(dds)), " genes."))
    dds <- DESeq(dds, parallel = T)
    save(dds, file=paste0("./DEA/dds/", DEA_projects[i], ".RData"))
    gc()
  }
}

# rm(list=grep("^projects$", ls(), invert = T, value = T))

for (i in seq_along(projects))
{
  load(paste0("./DEA/dds/", projects[i], ".RData"))
  # PCA plot
  print(paste0("PCA plot for project ", projects[i]))
  vsd <- vst(dds, blind=FALSE)
  p <- plotPCA(vsd, intgroup=c("condition"))+
    ggtitle(projects[i])
  p$layers[[1]]$aes_params$alpha <- 0.5  # hack for adding transparency
  ggsave(filename=paste0("./DEA/PCA/", projects[i], ".pdf"),
         plot = p, device = "pdf", height = 6, width = 6, units = "in", dpi = "retina")
  # Comparison between each stage
  ## matrix with each column as one experiment, e.g. N vs. I
  comparisons <- utils::combn(c("N", "I", "II", "III", "IV"), 2)
  for (j in seq(ncol(comparisons)))
  {
    print(paste0("Getting results for ", comparisons[2,j], " vs. ", comparisons[1,j]))
    res <- results(dds, contrast = c("condition",comparisons[2,j], comparisons[1,j]))
    res <- res[order(res$padj), ]
    save(res, file = paste0("./DEA/res/", projects[i], "_",
                            comparisons[2,j], "_vs_",
                            comparisons[1,j], ".RData"))

    # Volcano plot of protein coding genes
    print(paste0("Volcano plots for ", comparisons[2,j], " vs. ", comparisons[1,j]))
    rownames(res) <- sapply(strsplit(rownames(res), "\\."),
                            function(x) x[[1]][1])
    res <- data.frame(res@listData, row.names = rownames(res),
                      stringsAsFactors = F)
    res$ENSEMBL <- rownames(res)
    eg <- clusterProfiler::bitr(rownames(res), fromType="ENSEMBL", toType="SYMBOL",
                                OrgDb="org.Hs.eg.db")
    res <- dplyr::inner_join(res, eg, by="ENSEMBL")
    # Mark TOP_N_IN_VOLCANO genes
    selected_symbols <- dplyr::top_n(res, -TOP_N_IN_VOLCANO_ADJP, padj)$SYMBOL
    res_fc <- res[!is.na(res$padj) & res$padj <= VOLCANO_PADJ_CUTOFF, ]
    res_fc <- res_fc[order(res_fc$log2FoldChange), ]
    selected_symbols <- unique(c(
      selected_symbols,
      head(res_fc$SYMBOL, TOP_N_IN_VOLCANO_FC),
      tail(res_fc$SYMBOL, TOP_N_IN_VOLCANO_FC)
      ))
    p <- EnhancedVolcano(res, lab = res$SYMBOL,
                         x = "log2FoldChange",  y = "padj",
                         selectLab = selected_symbols,
                         pCutoff = VOLCANO_PADJ_CUTOFF, FCcutoff = 2.0,
                         colAlpha = 0.5,
                         transcriptLabSize = 3.0,
                         legendPosition = "bottom",
                         legendLabSize = 10, legendIconSize = 3.0,
                         DrawConnectors = T, widthConnectors = 0.2,
                         border = "full",
                         borderWidth = 1.5, borderColour = "black",
                         gridlines.major = FALSE,
                         gridlines.minor = FALSE)
    save_width <- ceiling(max(abs(res$log2FoldChange)))
    ggsave(filename=paste0("./DEA/volcano/", projects[i],
                           "_",comparisons[2,j], "_vs_",
                           comparisons[1,j], ".pdf"),
           plot = p, device = "pdf", width = max(10, save_width), height = 10, units = "in", dpi = "retina")
  }
}

# sessionInfo()
# R version 3.5.1 (2018-07-02)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux 9 (stretch)
#
# Matrix products: default
# BLAS: /usr/lib/openblas-base/libblas.so.3
# LAPACK: /usr/lib/libopenblasp-r0.2.19.so
#
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
#
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] dplyr_0.7.8                 clusterProfiler_3.10.0      EnhancedVolcano_1.0.0       ggrepel_0.8.0               ggplot2_3.1.0               DESeq2_1.22.1
# [7] SummarizedExperiment_1.12.0 DelayedArray_0.8.0          BiocParallel_1.16.0         matrixStats_0.54.0          Biobase_2.42.0              GenomicRanges_1.34.0
# [13] GenomeInfoDb_1.18.1         IRanges_2.16.0              S4Vectors_0.20.1            BiocGenerics_0.28.0
#
# loaded via a namespace (and not attached):
#   [1] fgsea_1.8.0            colorspace_1.3-2       ggridges_0.5.1         qvalue_2.14.0          htmlTable_1.12         XVector_0.22.0         base64enc_0.1-3
# [8] rstudioapi_0.8         farver_1.0             urltools_1.7.1         bit64_0.9-7            AnnotationDbi_1.44.0   xml2_1.2.0             splines_3.5.1
# [15] GOSemSim_2.8.0         geneplotter_1.60.0     knitr_1.20             Formula_1.2-3          jsonlite_1.5           annotate_1.60.0        cluster_2.0.7-1
# [22] GO.db_3.7.0            ggforce_0.1.3          compiler_3.5.1         httr_1.3.1             rvcheck_0.1.1          backports_1.1.2        assertthat_0.2.0
# [29] Matrix_1.2-15          lazyeval_0.2.1         tweenr_1.0.0           acepack_1.4.1          htmltools_0.3.6        prettyunits_1.0.2      tools_3.5.1
# [36] bindrcpp_0.2.2         igraph_1.2.2           gtable_0.2.0           glue_1.3.0             GenomeInfoDbData_1.2.0 reshape2_1.4.3         DO.db_2.9
# [43] fastmatch_1.1-0        Rcpp_1.0.0             enrichplot_1.2.0       ggraph_1.0.2           stringr_1.3.1          XML_3.98-1.16          DOSE_3.8.0
# [50] europepmc_0.3          zlibbioc_1.28.0        MASS_7.3-51.1          scales_1.0.0           hms_0.4.2              RColorBrewer_1.1-2     yaml_2.2.0
# [57] memoise_1.1.0          gridExtra_2.3          UpSetR_1.3.3           triebeard_0.3.0        rpart_4.1-13           latticeExtra_0.6-28    stringi_1.2.4
# [64] RSQLite_2.1.1          genefilter_1.64.0      checkmate_1.8.5        rlang_0.3.0.1          pkgconfig_2.0.2        bitops_1.0-6           lattice_0.20-38
# [71] purrr_0.2.5            bindr_0.1.1            htmlwidgets_1.3        cowplot_0.9.3          bit_1.1-14             tidyselect_0.2.5       plyr_1.8.4
# [78] magrittr_1.5           R6_2.3.0               Hmisc_4.1-1            DBI_1.0.0              pillar_1.3.0           foreign_0.8-71         withr_2.1.2
# [85] units_0.6-1            survival_2.43-1        RCurl_1.95-4.11        nnet_7.3-12            tibble_1.4.2           crayon_1.3.4           viridis_0.5.1
# [92] progress_1.2.0         locfit_1.5-9.1         grid_3.5.1             data.table_1.11.8      blob_1.1.1             digest_0.6.18          xtable_1.8-3
# [99] tidyr_0.8.2            gridGraphics_0.3-0     munsell_0.5.0          ggplotify_0.0.3        viridisLite_0.3.0

