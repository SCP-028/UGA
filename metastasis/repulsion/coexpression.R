# https://bioconductor.org/packages/devel/bioc/vignettes/CVE/inst/doc/WGCNA_from_TCGA_RNAseq.html
# Co-expression analysis of cadherin genes and sialic acid genes.
try(setwd("C:/Users/jzhou/Desktop/expression_count"), silent=T)
try(setwd("/home/yizhou/data/expression_count"), silent=T)
if(!require(fastcluster)) {
    install.packages(fastcluster)
    library(fastcluster)
}

# Load data
load("BRCA.RData")
stageAnnot <- data.table::fread("./annotation/annot.tsv", data.table=F,
                                stringsAsFactors=F, header=T)
load("../sialic_acid/Sialic_acids_marker.RData")
#cadGenes <- c("P?CDH[AB]?\\d{1,2}", "DSG[1234]")  # cadherin marker genes
cadGenes <- c("CDH1", "CDH2", "CDH3", "CDH4", "CDH5", "CDH6", "CDH7", "CDH8",
              "CDH10", "CDH11", "CDH12", "CDH13", "CDH15", "CDH16", "CDH17",
              "CDH18", "CDH19", "CDH20", "CDH22", "CDH23", "CDH24", "CDH26",
              "DSG1", "DSG1-AS1", "DSG2", "DSG2-AS1", "DSG3", "DSG4",
              "PCDH1", "PCDH7", "PCDH8", "PCDH8P1", "PCDH9", "PCDH9-AS2",
              "PCDH10", "PCDH11X", "PCDH11Y", "PCDH12", "PCDH15", "PCDH17",
              "PCDH18", "PCDH19", "PCDH20", "PCDHA1", "PCDHA2", "PCDHA3",
              "PCDHA4", "PCDHA5", "PCDHA6", "PCDHA7", "PCDHA8", "PCDHA9",
              "PCDHA10", "PCDHA11", "PCDHA12", "PCDHA13",
              "PCDHB1", "PCDHB2", "PCDHB3", "PCDHB4", "PCDHB5", "PCDHB6",
              "PCDHB7", "PCDHB8", "PCDHB9", "PCDHB10", "PCDHB11", "PCDHB12",
              "PCDHB13", "PCDHB14", "PCDHB15", "PCDHB16", "PCDHB17P",
              "PCDHB18P", "PCDHB19P"
              )
genelist <- data.frame(fromNode=c(sia.sia, sia.gal, sia.galnac,  # sialic transferases
                                  "SLC35A1", "SLC17A5",  # more transferases
                                  siaP,  # sialic acid related proteins
                                  siaEnzy, sialidases,  # sialic acids synthesis & sialidases
                                  cadGenes),
                       category=c(rep("transferase_sialic", length(sia.sia)),
                                  rep("transferase_gal", length(sia.gal)),
                                  rep("transferase_galnac", length(sia.galnac)),
                                  rep("co_transferase", 2),
                                  rep("sialic_acid_related", length(siaP)),
                                  rep("sialic_acid_synthesis", length(siaEnzy)),
                                  rep("sialidase", length(sialidases)),
                                  rep("cadherin", length(cadGenes))))
genelist <- genelist[!duplicated(genelist$fromNode), ]

# Functions
source("../DEG_func.R")

edgeR.de.test <- function(df1, df2, group1, group2) {
    #' Use edgeR to perform gene differential expression analysis.
    #'
    #' Require edgeR and dplyr to work.
    #' 
    #' @param df1 First data frame / matrix. Must be counts value!!
    #' @param df2 Second. The result comes as df2:df1.
    #' @param group1 Name of first group.
    #' @param group2 Name of second group.
    library(edgeR)
    df1 <- df1[order(rownames(df1)), ]
    df2 <- df2[order(rownames(df2)), ]
    group <- c(rep(group1, ncol(df1)), rep(group2, ncol(df2)))
    design <- model.matrix(~0+group)
    x <- cbind(df1, df2)
    y <- DGEList(counts=x, group=group)
    y <- y[rowSums(cpm(y) > 1) >= 2, , keep.lib.sizes=F]
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design=design)
    ## logFC logCPM PValue
    et <- exactTest(y)$table
    colnames(et) <- c("log2_fold_change", "log2_CPM", "p_value")
    et$ensembl <- rownames(et)
    rownames(et) <- NULL
    return(et)
}


# Separate before stage iii & stage iii, iv
stage1 <- stageAnnot$barcode[grepl("(^i$)|(\\si[abc]?$)|(1)", stageAnnot$tumor_stage) &
                            stageAnnot$project == "BRCA"]
stage2 <- stageAnnot$barcode[grepl("(^ii$)|(\\si{2}[abc]?$)|(2)", stageAnnot$tumor_stage) &
                            stageAnnot$project == "BRCA"]
stage3 <- stageAnnot$barcode[grepl("(^iii$)|(\\si{3}[abc]?$)|(3)", stageAnnot$tumor_stage) &
                            stageAnnot$project == "BRCA"]
stage4 <- stageAnnot$barcode[grepl("(^iv$)|(\\siv[abc]?$)|(4)", stageAnnot$tumor_stage) &
                            stageAnnot$project == "BRCA"]
datan <- cbind.data.frame(datan, datat[ ,colnames(datat) %in% c(stage1, stage2)])
datat <- datat[ ,colnames(datat) %in% c(stage3, stage4)]
ensemblAnnot <- retrieve.ensembAnnot(datan)
datan <- name.convert(datan, ensemblAnnot)
datat <- name.convert(datat, ensemblAnnot, libSize=T)
libt <- datat$libSize
datat <- datat$exp

# Get differentially expressed genes
deg <- edgeR.de.test(datan, datat, "before_stage_iii", "stage_iii_iv")
colnames(deg)[4] <- "symbol"
deg <- deg[(deg$symbol %in% genelist$fromNode)|
           ((deg$p_value <= 0.01) & (deg$log2_fold_change >= 1)), ]

# Normalize RNA-Seq data
require(limma)
require(dplyr)
df <- datat[rownames(datat) %in% deg$symbol, ]
df <- voom(df, lib.size=libt)$E  # generate weight for each observation

# Take top-k most varient genes using median absolute devision (mad)
df <- t(df[order(apply(df,1,mad), decreasing = T), ])

# Biweight midcorrelation to define an similarity matrix
s <- abs(bicor(df))

# Pick lowest power for which the scale-free topology fit index curve flattens
# out upon reaching a high value (R^2 >= 0.9)
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(df, powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab='Soft Threshold (power)',
     ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence'))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col='red'); abline(h=0.90,col='red')
beta <- 4

# Calculate adjacency matrix
a <- s ^ beta
# Dissimilarity measure
w <- 1 - a

df <- exportNetworkToCytoscape(a, threshold=0.5)$edgeData
df[ ,c("fromAltName", "toAltName")] <- NULL
df <- left_join(df, genelist, by="fromNode")
df$category <- as.character(df$category)
df$category[is.na(df$category)] <- "deg"
write.table(df, file="./sia_cad_edge.txt", row.names=F, quote=F)