library(fgsea)

setwd("~/data/mice_liver")
set.seed(0)

# Gene sets are downloaded from http://software.broadinstitute.org/gsea/msigdb/collections.jsp#H
# C2, C5 and C6 gene sets are used
pathway_files <- list.files("~/data/annotation/MSigDB", full.names = T)
if ()
pathways <- vector("list")
for (f in pathway_files)
{
  pathways <- c(pathways, gmtPathways(f))
}

up_files <- list.files("./DEA_result", pattern=".*up.*")
down_files <- sub("up", "down", up_files)

for (i in seq_along(up_files))
{
  save_name <- sub("_up", "", strsplit(up_files[i], "\\.")[[1]][1])
  df <- data.table::fread(paste0("./DEA_result/", up_files[i]))
  df_ <- data.table::fread(paste0("./DEA_result/", down_files[i]))
  df <- rbind.data.frame(df, df_)
  # Filter for genes that are significantly differentially expressed
  df <- df[!is.na(df$padj) & df$padj <= 0.01, ]

  # Use log2(fold change) as gene-level stat
  geneList <- df$log2FoldChange
  names(geneList) <- toupper(df$`gene name`)

  ans <- fgsea(pathways=pathways, stats=geneList,
               nperm = 10e5, minSize = 10, maxSize = 500,
               nproc = 20)
  ans <- ans[order(ans$padj), ]

  # Plot top pathways
  topPathwaysUp <- ans[ES > 0][head(order(padj), n=10), pathway]
  topPathwaysDown <- ans[ES < 0][head(order(padj), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  p <- plotGseaTable(pathways[topPathways], geneList, ans,
                gseaParam = 0.5)
  ggplot2::ggsave(plot=p, filename=paste0("plot/", save_name, ".top.png"),
                  device = "png", width = 6, height=12, units = "in", dpi = "retina")

  # Plot collapsed pathways
  collapsedPathways <- collapsePathways(ans[order(pval)][padj < 0.01],
                                        pathways, geneList)
  mainPathways <- ans[pathway %in% collapsedPathways$mainPathways][
    order(-NES), pathway]
  p <- plotGseaTable(pathways[mainPathways], geneList, ans,
                gseaParam = 0.5)
  ggplot2::ggsave(plot=p, filename=paste0("plot/", save_name, ".collapsed.png"),
                  device = "png", width = 6, height=20, units = "in", dpi = "retina")

  # Save results
  data.table::fwrite(ans, file=paste0("GSEA/", save_name, ".tsv"),
                     sep="\t", sep2=c("", ",", ""))
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
#   [1] stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] fgsea_1.8.0 Rcpp_1.0.0
#
# loaded via a namespace (and not attached):
#   [1] rstudioapi_0.8      bindr_0.1.1         magrittr_1.5        tidyselect_0.2.5    munsell_0.5.0       BiocParallel_1.16.0 lattice_0.20-38     colorspace_1.3-2
# [9] R6_2.3.0            rlang_0.3.0.1       fastmatch_1.1-0     plyr_1.8.4          dplyr_0.7.8         tools_3.5.1         parallel_3.5.1      grid_3.5.1
# [17] data.table_1.11.8   gtable_0.2.0        yaml_2.2.0          lazyeval_0.2.1      assertthat_0.2.0    tibble_1.4.2        crayon_1.3.4        Matrix_1.2-15
# [25] bindrcpp_0.2.2      gridExtra_2.3       purrr_0.2.5         ggplot2_3.1.0       glue_1.3.0          labeling_0.3        compiler_3.5.1      pillar_1.3.0
# [33] scales_1.0.0        pkgconfig_2.0.2
