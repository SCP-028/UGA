# Find out cadherin-related genes that are differentially expressed in
# different stages of cancer, and calculate the correlation of these
# genes and sialic acid related genes.
try(setwd("C:/Users/jzhou/Desktop/expression_FPKM"), silent=T)
try(setwd("/home/yizhou/data/expression_FPKM"), silent=T)
source("./DEG_func.R")  # all the functions needed
annot <- read.csv("./annotation/annot.tsv")  # stage information
load("../sialic_acid/Sialic_acids_marker.RData")  # sialic acid marker genes
filen <- list.files("./", pattern = "normal.csv$")
filet <- list.files("./", pattern = "tumor.csv$")
project <- intersect(sub("(.*)_.*$", "\\1", filen),
                     sub("(.*)_.*$", "\\1", filet))
df <- data.table::fread(paste0(project[1], "_normal.csv"), header = T,
                               stringsAsFactors = F, sep = "\t", data.table = F)
ensAnnot <- retrieve.ensembAnnot(df)  # ensembl -> gene symbol
rm(df, filen, filet)

for (i in seq_along(project)) {
    datan <- data.table::fread(paste0(project[i], "_normal.csv"), header = T,
                               stringsAsFactors = F, sep = "\t", data.table = F)
    datat <- data.table::fread(paste0(project[i], "_tumor.csv"), header = T,
                               stringsAsFactors = F, sep = "\t", data.table = F)
    if (ncol(datan) <= 5 | ncol(datat) <= 5) {
        message(paste0(project[i], " ignored because it has fewer than 5 cases."))
    }
    else {
        # Boxplot of cadherin-related genes
        dfn <- name.convert(datan, ensAnnot,
                            regex = "(^P?CDH[AB]?\\d{1,2}$)|(^DSG[1234]$)")
        # CDHs, PCDHs, DSG
        dft <- name.convert(datat, ensAnnot,
                            regex = "(^P?CDH[AB]?\\d{1,2}$)|(^DSG[1234]$)")
        cadGenes <- exp.boxplot(dfn, dft, proj=project[i], ensAnnot, annot)
        # Heatmap of cadherin ~ sialic acid genes
        genelist <- c(sia.sia, sia.gal, sia.galnac,  # sialiac transferases
                      siaP,  # sialic acid related proteins
                      siaEnzy, sialidases,  # sialic acids synthesis & sialidases
                      cadGenes)
        genelist <- unique(genelist)
        dfn <- name.convert(datan, ensAnnot, genelist=genelist)
        dft <- name.convert(datat, ensAnnot, genelist=genelist)
        dfn <- prep.cor.heatmap(dfn, genelist, sample="normal", cluster=F)
        dft <- prep.cor.heatmap(dft, genelist, sample="tumor", cluster=F)
        pltn <- nice.heatmap(dfn, paste(project[i], "normal samples", sep=" "))
        pltt <- nice.heatmap(dft, paste(project[i], "tumor samples", sep=" "))
        tiff(filename=paste0("../sialic_acid/", project[i], ".tiff"),
             width=16, height=9, units="in", res=200)
        multiplot(dfn, dft, cols=2)
        dev.off()
        message(paste0(project[i], " finished!"))
    }
}
