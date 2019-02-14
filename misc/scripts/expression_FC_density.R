library(tidyverse)
library(ggpubr)

# Load datasets --------------------------------------------------------------
rm(list=ls())
projects <- c("LIHC", "LUAD", "ESCA", "STAD", "LUSC", "COAD", "HNSC",
              "KIRP","BLCA", "KIRC", "BRCA", "KICH", "THCA", "PRAD")
projects <- projects[order(projects)]
projects <- paste0("TCGA-", projects)

annot <- data.table::fread(
  "/home/yi/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  select(ensembl = ensembl_gene_id,
         symbol = external_gene_name
  )
# Change these for other gene lists ==========================================
kinases <- data.table::fread(
  "/home/yi/storage/data/metabolic_reprogramming/human_kinases.txt",
  header = F)
kinases <- as.character(unique(kinases$V1))

phosphatases <- data.table::fread(
  "/home/yi/storage/data/metabolic_reprogramming/DEPOD_201410_human_phosphatases.csv"
)
phosphatases <- as.character(unique(phosphatases$`Gene symbol`))

# The gene.lists variable will be used later for plotting
gene.lists <- list(kinases = kinases, phosphatases = phosphatases)



# Don't have to modify anything below this line ==============================

# Density Plots --------------------------------------------------------------
for (i in seq_along(projects)) {
  # Load DEA result ----------------------------------------------------------
  deg <- data.table::fread(
    str_glue("/home/yi/storage/data/TCGA/TumorVsNormal/csv/{projects[i]}.csv")
  ) %>%
    # Ensembl ID to gene symbol
    rename(ensembl = V1) %>%
    mutate(ensembl = gsub("\\.\\d+$", "", ensembl)) %>%
    inner_join(annot, by = "ensembl") %>%
    # For genes with the same symbol, take the one with a lower padj
    filter(padj <= 0.01) %>%
    arrange(padj) %>%
    distinct(symbol, .keep_all = T) %>%
    mutate(
      Normal = baseMean,
      Tumor = baseMean * (2 ^ log2FoldChange)) %>%
    select(symbol, Normal, Tumor, log2FoldChange) %>%
    # Long to wide for plotting
    gather(key = sample, value = value, -symbol)

  # Density plots with expression on the x-axis ------------------------------
  for (j in seq_along(gene.lists)) {
    df <- deg %>% filter(symbol %in% gene.lists[[j]])

    # Plot Tumor vs. Normal
    p1 <- ggdensity(df[df$sample != "log2FoldChange", ],
                    x = "value", add = "mean", rug = TRUE,
                    color = "sample")
    p1 <- ggpar(p1, palette = "jco", legend = "right",
                title = str_glue(
                  "{projects[i]} Tumor vs. Normal for {names(gene.lists)[[j]]}")
    )
    # Plot log2FoldChange
    p2 <- ggdensity(df[df$sample == "log2FoldChange", ],
                    x = "value", add = "mean", rug = TRUE)
    p2 <- ggpar(p2, legend = "right",
                title = str_glue(
                  "{projects[i]} log2(FC) for {names(gene.lists)[[j]]}")
    )

    p <- ggarrange(p1, p2)

    dir.create(file.path("/home/yi/storage/data/metabolic_reprogramming/",
                         names(gene.lists)[[j]]), showWarnings = FALSE)
    ggsave(p, filename = str_glue(
      "/home/yi/storage/data/metabolic_reprogramming/{names(gene.lists)[[j]]}/{projects[i]}.png"),
           device = "png", height = 5, width = 15, units = "in", dpi = "retina")
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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C              LC_PAPER=en_US.UTF-8       LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
#
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] ggpubr_0.2      magrittr_1.5    forcats_0.3.0   stringr_1.4.0   dplyr_0.7.8     purrr_0.3.0     readr_1.3.1     tidyr_0.8.2     tibble_2.0.1    ggplot2_3.1.0   tidyverse_1.2.1
#
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.0       cellranger_1.1.0 pillar_1.3.1     compiler_3.5.1   plyr_1.8.4       bindr_0.1.1      tools_3.5.1      jsonlite_1.6     lubridate_1.7.4  nlme_3.1-137     gtable_0.2.0     lattice_0.20-38
# [13] pkgconfig_2.0.2  rlang_0.3.1      cli_1.0.1        rstudioapi_0.9.0 yaml_2.2.0       haven_2.0.0      bindrcpp_0.2.2   withr_2.1.2      xml2_1.2.0       httr_1.4.0       generics_0.0.2   hms_0.4.2
# [25] grid_3.5.1       tidyselect_0.2.5 glue_1.3.0       R6_2.3.0         readxl_1.2.0     modelr_0.1.3     backports_1.1.3  scales_1.0.0     rvest_0.3.2      assertthat_0.2.0 colorspace_1.4-0 stringi_1.2.4
# [37] lazyeval_0.2.1   munsell_0.5.0    broom_0.5.1      crayon_1.3.4

