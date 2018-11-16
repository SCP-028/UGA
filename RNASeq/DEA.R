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
         plot = p, device = "pdf",
         height = 6, width = 6, units = "in", dpi = "retina")
  # Comparison between each stage
  ## matrix with each column as one experiment, e.g. N vs. I
  comparisons <- utils::combn(c("N", "I", "II", "III", "IV"), 2)
  for (j in seq(ncol(comparisons)))
  {
    print(paste0("Getting results for ", comparisons[2][j], " vs. ", comparisons[1][j]))
    res <- results(dds, contrast = c("condition",
                                     comparisons[2][j], comparisons[1][j])
                   )
    res <- res[order(res$padj), ]
    save(res, file = paste0("./DEA/res/", projects[i], "_",
                            comparisons[2][j], "_vs_",
                            comparisons[1][j], ".RData"))

    # Volcano plot of protein coding genes
    print(paste0("Volcano plots for ", comparisons[2][j], " vs. ", comparisons[1][j]))
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
    ggplot2::ggsave(filename=paste0("./DEA/volcano/", projects[i],
                                    "_",comparisons[2][j], "_vs_",
                                    comparisons[1][j], ".pdf"),
                    plot = p, device = "pdf",
                    width = save_width, height = 10, units = "in", dpi = "retina")
  }
}



