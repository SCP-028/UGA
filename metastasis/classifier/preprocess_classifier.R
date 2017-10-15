# Use data downloaded by `TCGA_expression_pipeline.py`
# and generate DEG results using edgeR.
# Is required for `*classifier.py`.
library(dplyr)
if (!file.exists("./COAD_iv_i_iii.RData")) {
    # Follow steps in ecm_gene_de.R
    message("Didn't find DEG file, loading counts data...")
    library(edgeR)
    source("./DEG_func.R")
    df <- data.table::fread("./expression_count/COAD_tumor.csv",
                            data.table=F, header=T, sep="\t",
                            stringsAsFactors=F)
    df <- name.convert(df, description=T)
    symbolDescription <- cbind.data.frame(rownames(df), df$description)
    colnames(symbolDescription) <- c("symbol", "description")
    df$description <- NULL
    annot <- read.csv("./expression_count/annotation/annot.tsv", header=T)
    annot <- annot[annot$project == "COAD" & annot$sample_type == "tumor", ]
    df_iv <- df[ ,colnames(df) %in% annot[grep("iv", annot$tumor_stage), "barcode"]]
    df_i_iii <- df[ ,colnames(df) %in% annot[grep("\\si{1,3}[abc]?", annot$tumor_stage), "barcode"]]
    message("Finding differentially expressed genes...")
    result <- edgeR.de.test(df_i_iii, df_iv)
    result <- inner_join(result, symbolDescription, by = "symbol")
    save(result, file="./COAD_iv_i_iii.RData")
}

rm(list=ls())
load("./COAD_iv_i_iii.RData")  # edgeR `result`
source("./DEG_func.R")

df <- data.table::fread("./expression_FPKM/COAD_tumor.csv",
                        data.table=F, header=T, sep="\t")
df <- name.convert(df, description=F)
df <- get.DEG(df, result)

write.csv(df, "./COAD_DEG.csv", row.names=F)
