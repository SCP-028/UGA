## Compare the expression of certain genes in different stages of BRCA, and draw boxplots.
# Works on RStudio server
setwd("~/data")
set.seed(1005)
# source("https://raw.githubusercontent.com/ggrothendieck/gsubfn/master/R/list.R")
library(dplyr)
library(biomaRt)
library(edgeR)
library(reshape2)
library(ggplot2)

REGEX <- ""  # default ''
GENELIST <- c("ALCAM", "ADAMTS[14579]3?$", "ACAN", "COMP", "IGF1$", "TGFB1$", "VEGFA",
              "EGFR$", "CYR61", "NOV$", "COL\\d", "P4HB", "DSPP", "ELN", "FBN", "FN1",
              "IBSP", "LAMA", "LAMC", "LTBP", "MATN", "SPP1", "BGN", "DCN", "FMOD",
              "LUM", "TNXB", "TNC", "VCAN$", "VTN")  # default ''

source("./DEG_func.R")

# Load BRCA expression data (counts) and annotation file
dfCount <- data.table::fread("./expression_count/BRCA_tumor.csv",
                              header = T, sep = "\t", data.table = F,
                              stringsAsFactors = F)
annot <- read.csv("./expression_count/annotation/annot.tsv", header=T)
annot <- annot[annot$project == "BRCA" & annot$sample_type == "tumor", ]

# Separate stage i-iii and stage iv data
ensAnnot <- retrieve.ensembAnnot(dfCount)
dfCount <- name.convert(dfCount, ensAnnot, genelist=GENELIST, description=T)
symbolDescription <- cbind.data.frame(rownames(dfCount), dfCount$description)
colnames(symbolDescription) <- c("symbol", "description")
dfCount$description <- NULL
dfCount_iv <- dfCount[ ,colnames(dfCount) %in% annot[grep("iv", annot$tumor_stage), "barcode"]]
dfCount_i_iii <- dfCount[ ,colnames(dfCount) %in% annot[(!grepl("v|x|0|4|not", annot$tumor_stage)) &annot$tumor_stage != '', "barcode"]]

# Differential expression test of stage i-iii and stage iv using `edgeR`
message("Finding differentially expressed genes...")
result <- edgeR.de.test(dfCount_i_iii, dfCount_iv)
result <- inner_join(result, symbolDescription, by = "symbol")

# Load FPKM data
dfFpkm <- data.table::fread("./expression_FPKM/BRCA_tumor.csv",
                             header = T, sep = "\t", data.table = F,
                             stringsAsFactors = F)

# Convert gene name and stage info
dfFpkm <- name.convert(dfFpkm, ensAnnot, genelist=GENELIST)
dfFpkm_iv <- dfFpkm[ ,colnames(dfFpkm) %in% colnames(dfCount_iv)]
dfFpkm_i_iii <- dfFpkm[ ,colnames(dfFpkm) %in% colnames(dfCount_i_iii)]
dfFpkm_iv <- transform.for.plot(dfFpkm_iv, 'stage_iv')
dfFpkm_i_iii <- transform.for.plot(dfFpkm_i_iii, 'stage_i_iii')
dfPlot <- rbind.data.frame(dfFpkm_iv, dfFpkm_i_iii)

# Keep only genes that have a differential expression
dfPlot <- dfPlot[dfPlot$symbol %in% resultFew$symbol, ]
dfPlot$symbol <- factor(dfPlot$symbol, levels = resultFew$symbol)

# Draw boxplot of the expression of the genes
ggplot(dfPlot, aes(symbol, log2_FPKM, fill=stage)) +
    geom_boxplot() +
    theme_grey(base_size = 9) +
    theme(axis.text.x = element_text(size = 8, angle = 330, hjust = 0)) +
    ggtitle("Expression of ECM genes")

ggsave(filename = "./ECM.png", width = 12, height = 9, units = "in")
save(result, resultFew, file = "./ECM.RData")