library(KEGGREST)
library(data.table)
library(tidyverse)
library(ggpubr)
library(cowplot)

# Get all reactions associated with ATP (C00002) -----------------------------

ATP_info <- keggGet("cpd:C00002")
rxns <- ATP_info[[1]]$REACTION

# Put reactions into a list of characters
rxns <- unlist(sapply(rxns, function(x) strsplit(x, " ")), use.names = F)
rxns <- paste0("rn:", rxns)

res <- matrix(nrow = length(rxns), ncol = 3)
colnames(res) <- c("reaction", "equation", "enzyme")
# Request 10 reactions at a time due to API limitations
for (i in seq(from=1, to=length(rxns), by=10)) {
  message(stringr::str_glue("Getting reactions {i} to {min(i + 9, length(rxns))}."))
  rxn <- keggGet(rxns[i:min(i + 9, length(rxns))])
  for (j in seq_along(rxn)) {
    res[i+ j - 1, ] <- c(rxn[[j]]$ENTRY,
                         rxn[[j]]$EQUATION,
                         paste(rxn[[j]]$ENZYME, collapse = ","))
  }
}
res <- as.data.table(res)

# 3 columns: reaction, equation, enzyme --------------------------------------

# reactions without EC number
rxn_wo_ec <- res[enzyme == ""]

# reactions with EC numbers like 2.7.2.-
rxn_w_ec_class <- res[grep("-", enzyme)]

# reactions with multiple EC numbers; ones with EC classes are excluded
rxn_w_multiple_ec <- res[grep(",", enzyme)][!reaction %in% rxn_w_ec_class$reaction]

# reactions with only one EC number
rxn_one_ec <- res[grep("^(\\d+\\.){3}\\d+$", enzyme)]

# Reactions with 1 EC number: get matching gene names in human ---------------
# The for loop is really slow!

# Function for getting enzyme - gene information
annotate.ec <- function(df) {
  ECres <- matrix(nrow = nrow(df), ncol = 3)
  colnames(ECres) <- c("EC", "name", "genes")
  df$enzyme <- paste0("ec:", df$enzyme)
  for (i in seq(from=1, to=nrow(df), by=10)) {
    message(stringr::str_glue("Getting enzymes {i} to {min(i + 9, nrow(df))}."))
    ec <- keggGet(df$enzyme[i:min(i + 9, nrow(df))])
    for (j in seq_along(ec)) {
      # Skip entries without matching gene names in human
      ecGene <- ec[[j]]$GENES[grep("^HSA: ", ec[[j]]$GENES)]
      if (length(ecGene) != 0) {
        ECres[i + j - 1, ] <- c(
          ec[[j]]$ENTRY,
          paste(ec[[j]]$NAME, collapse = " "),
          ec[[j]]$GENES[grep("^HSA: ", ec[[j]]$GENES)]
        )
      }
    }
  }
  df <- as.data.table(cbind(df, ECres))
  df <- df[!is.na(genes),.(reaction, equation, EC, name, genes)]
  return(df)
}

hsa_rxn_one_ec <- annotate.ec(rxn_one_ec)

# Reactions with multiple EC numbers: split them and get annotations ---------
#
# for each entry
rxn_w_multiple_ec <- rxn_w_multiple_ec %>%
  mutate(enzyme = strsplit(as.character(enzyme), ",")) %>%
  unnest() %>%
  filter(enzyme != "")


hsa_rxn_w_multiple_ec <- annotate.ec(rxn_w_multiple_ec)
# Save temp results as they take a while to run
save(res, hsa_rxn_one_ec, hsa_rxn_w_multiple_ec,
     file="./storage/R_script/KEGG/KEGG_ATP_reactions.RData")

# Convert entrez IDs to Ensembl gene IDs -------------------------------------

rm(list=ls())
load("/home/yi/storage/R_script/KEGG/KEGG_ATP_reactions.RData")
annot <- data.table::fread("/home/yi/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv")

df <- rbind.data.frame(hsa_rxn_one_ec, hsa_rxn_w_multiple_ec) %>%
  mutate(genes = gsub("^HSA: ", "", genes)) %>%
  mutate(genes = strsplit(as.character(genes), " ")) %>%
  unnest() %>%
  filter(genes != "") %>%
  mutate(entrez = gsub("(^\\d+).*$", "\\1", genes),
         genes = gsub("^.*\\(([^(]+)\\)", "\\1", genes))

annot <- annot %>%
  select(entrez = entrezgene, ensembl = ensembl_gene_id) %>%
  filter(!is.na(entrez)) %>%
  mutate_if(is.numeric, as.character) %>%
  filter(entrez %in% df$entrez)

# Prepare df equations to get number of ATP(s) consumed / produced -----------
parse_equation <- function(df) {
  eq_left <- sapply(df$equation, function(x) strsplit(x, " <=> ")[[1]][1])
  # Only keep the substrate / product C00002 (ATP)
  eq_left <- sapply(strsplit(eq_left, " \\+ "), function(x) ifelse(
    length(x[grep("C00002", x)]) == 0, 0, x[grep("C00002", x)]))
  eq_left[eq_left == "C00002"] <- 1
  eq_left[grep("[^0-9]", eq_left)] <- sapply(eq_left[grep("[^0-9]", eq_left)],
                                             function(x) strsplit(x, " ")[[1]][1])

  eq_right <- sapply(df$equation, function(x) strsplit(x, " <=> ")[[1]][2])
  eq_right <- sapply(strsplit(eq_right, " \\+ "), function(x) ifelse(
    length(x[grep("C00002", x)]) == 0, 0, x[grep("C00002", x)]))
  eq_right[eq_right == "C00002"] <- 1

  res <- data.frame(lhs = eq_left, rhs = eq_right)
  res <- res %>%
    mutate_if(is.factor, as.character) %>%
    mutate(net_ATP = as.numeric(rhs) - as.numeric(lhs))  # NA introduced by n ATP and n-1 ATP
  return(res$net_ATP)
}

df <- df %>%
  inner_join(annot, by = "entrez") %>%
  select(-c(genes, entrez))
df$netATP <- parse_equation(df)

hsa_rxn_w_nATP <- df[is.na(df$netATP), ]
rxn_wo_ec <- res[enzyme == ""]
rxn_w_ec_class <- res[grep("-", enzyme)]
save(rxn_w_ec_class, rxn_wo_ec, hsa_rxn_w_nATP, file="error_rxns.RData")
df <- df[!is.na(df$netATP), ]
data.table::fwrite(df, file="./hsa_ATP_rxns.csv", row.names = F, col.names = T)

# Load corresponding DEA results ---------------------------------------------
rm(list=ls())
df <- data.table::fread("hsa_ATP_rxns.csv")
projects <- list.files("/home/yi/storage/data/TCGA/DEA/csv/")
projects <- gsub("_.*$", "", projects)
projects <- unique(projects)

comparisons <- c("II", "III", "IV")

# Density Plots and scatter plots --------------------------------------------
for (i in seq_along(projects)) {
  # Prepare data frame for plotting ------------------------------------------
  deg <- data.table::fread(
    str_glue("/home/yi/storage/data/TCGA/DEA/csv/{projects[i]}_I_vs_N.csv")
    ) %>%
    rename(ensembl = V1) %>%
    mutate(ensembl = gsub("\\.\\d+$", "", ensembl)) %>%
    inner_join(df, by = "ensembl") %>%
    select(reaction, EC, ensembl, netATP, baseMean, log2FoldChange) %>%
    mutate(
      normalMean = baseMean * netATP,
      tumorMean = baseMean * (2 ^ log2FoldChange) * netATP) %>%
    # Merge genes translated to the same EC
    group_by(reaction, EC) %>%
    summarise(normalMean = mean(normalMean),
              tumorMean = mean(tumorMean)) %>%
    # then merge ECs in the same reaction
    group_by(reaction) %>%
    summarise(Normal = mean(normalMean), I = mean(tumorMean)) %>%
    # long to wide for plotting
    gather("stage", "ATP", -reaction)

  for (j in seq_along(comparisons)) {
    tmp <- data.table::fread(
      str_glue("/home/yi/storage/data/TCGA/DEA/csv/{projects[i]}_{comparisons[j]}_vs_N.csv")
      ) %>%
      rename(ensembl = V1) %>%
      mutate(ensembl = gsub("\\.\\d+$", "", ensembl)) %>%
      inner_join(df, by = "ensembl") %>%
      select(reaction, EC, ensembl, netATP, baseMean, log2FoldChange) %>%
      mutate(tumorMean = baseMean * (2 ^ log2FoldChange) * netATP) %>%
      group_by(reaction, EC) %>%
      summarise(tumorMean = mean(tumorMean)) %>%
      group_by(reaction) %>%
      summarise(Tumor = mean(tumorMean)) %>%
      gather("stage", "ATP", -reaction) %>%
      mutate(stage = comparisons[j])
    deg <- rbind.data.frame(deg, tmp)
  }
  # Density plot with ATPs consumed on the x-axis ----------------------------
  ## replicate the "Normal" group in each facet for plotting
  degT <- deg %>%
    filter(stage != "Normal") %>%
    mutate(facet = stage)

  degN <-deg %>%
    filter(stage == "Normal") %>%
    left_join(data.table(stage = "Normal", facet = unique(degT$facet)),
              by = "stage")
  degDens <- rbind.data.frame(degN, degT)

  ## Density plot
  degDens$ATP <- -degDens$ATP  # larger value means consuming more ATP
  degDens$stage <- ifelse(degDens$stage == "Normal", "Normal", "Tumor")
  p <- ggdensity(degDens, x = "ATP",
                 add = "mean", rug = TRUE,
                 facet.by = "facet",
                 color = "stage", fill = "stage")
  p <- ggpar(p, palette = c("#00AFBB", "#E7B800"),
             title = projects[i], legend = "right")
  ggsave(p, filename = str_glue("/home/yi/storage/data/KEGG/{projects[i]}.png"),
         device = "png", height = 9, width = 21, units = "in", dpi = "retina")

  # Scatter plot with Normal on the x-axis and tumors on the y-axis ----------
  degN <- deg %>%
    filter(stage == "Normal") %>%
    select(-stage) %>%
    rename(Normal = ATP)
  degQQ <- deg %>%
    filter(stage != "Normal") %>%
    left_join(degN, by = "reaction") %>%
    rename(Tumor = ATP) %>%
    mutate(Tumor = -Tumor, Normal = -Normal)

  # Main plot
  pmain <- ggplot(degQQ, aes(x = Normal, y = Tumor, color = stage, shape = stage))+
    geom_point(alpha = 0.7, size = 1.2)+
    geom_abline(slope = 1, intercept = 0, aes(size = 0.7))+
    ggtitle(projects[i])+
    ggpubr::color_palette("jco")+
    theme_minimal()
  # Marginal densities along x axis
  xdens <- cowplot::axis_canvas(pmain, axis = "x")+
    geom_density(data = degQQ, aes(x = Normal, color = stage, fill = stage),
                 alpha = 0.7, size = 0.2)
  # Marginal densities along y axis
  # Need to set coord_flip = TRUE, if you plan to use coord_flip()
  ydens <- cowplot::axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
    geom_density(data = degQQ, aes(x = Tumor, color = stage, fill = stage),
                 alpha = 0.7, size = 0.2)+
    coord_flip()+
    ggpubr::fill_palette("jco")
  p1 <- cowplot::insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
  p2 <- cowplot::insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
  ggsave(p2, filename = str_glue("/home/yi/storage/data/KEGG/{projects[i]}_scatter.png"),
         device = "png", height = 10, width = 10, units = "in", dpi = "retina")
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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C              LC_PAPER=en_US.UTF-8
# [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
#
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] ggExtra_0.8       bindrcpp_0.2.2    forcats_0.3.0     stringr_1.3.1     dplyr_0.7.8       purrr_0.3.0       readr_1.3.1       tidyr_0.8.2       tibble_2.0.1      tidyverse_1.2.1
# [11] data.table_1.12.0 KEGGREST_1.22.0   ggpubr_0.2        magrittr_1.5      ggplot2_3.1.0
#
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.0          lubridate_1.7.4     lattice_0.20-38     png_0.1-7           Biostrings_2.50.2   assertthat_0.2.0    digest_0.6.18       mime_0.6            R6_2.3.0
# [10] cellranger_1.1.0    plyr_1.8.4          backports_1.1.3     stats4_3.5.1        httr_1.4.0          pillar_1.3.1        zlibbioc_1.28.0     rlang_0.3.1         lazyeval_0.2.1
# [19] readxl_1.2.0        rstudioapi_0.9.0    miniUI_0.1.1.1      S4Vectors_0.20.1    labeling_0.3        munsell_0.5.0       shiny_1.2.0         broom_0.5.1         httpuv_1.4.5.1
# [28] compiler_3.5.1      modelr_0.1.2        pkgconfig_2.0.2     BiocGenerics_0.28.0 htmltools_0.3.6     tidyselect_0.2.5    IRanges_2.16.0      later_0.7.5         crayon_1.3.4
# [37] withr_2.1.2         grid_3.5.1          xtable_1.8-3        nlme_3.1-137        jsonlite_1.6        gtable_0.2.0        scales_1.0.0        cli_1.0.1           stringi_1.2.4
# [46] XVector_0.22.0      promises_1.0.1      ggthemes_4.0.1      xml2_1.2.0          generics_0.0.2      cowplot_0.9.4       ggsci_2.9           tools_3.5.1         glue_1.3.0
# [55] hms_0.4.2           parallel_3.5.1      yaml_2.2.0          colorspace_1.4-0    BiocManager_1.30.4  rvest_0.3.2         bindr_0.1.1         haven_2.0.0
