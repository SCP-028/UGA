require(tidyr)
require(dplyr)
setwd("~/data/")
load("./methylation/BRCA.RData")
# load and clean probe annotation data
annot <- data.table::fread("./annotation/HumanMethylation450_15017482_v1-2.csv",
                           sep = ",", data.table = F, skip = 7, header = T,
                           stringsAsFactors = F,select = c("IlmnID",   
                                                            "UCSC_RefGene_Name",
                                                            "UCSC_RefGene_Group"))
annot <- annot[grep("(ST3GAL[126])|(ST6GALNAC[24])", annot$UCSC_RefGene_Name), ]
# annot <- data.frame(IlmnID = annot$IlmnID,
#                     gene = paste(annot$UCSC_RefGene_Name,
#                                  annot$UCSC_RefGene_Group, sep = "_"))
# annot <- unique(annot)
# annot <- separate_rows(annot, c("UCSC_RefGene_Name", "UCSC_RefGene_Group"), sep = ";")
# load stage annotation data
stageAnnot <- data.table::fread("./expression_FPKM/annotation/annot.tsv",
                                sep = ",", header = T, data.table = F,
                                stringsAsFactors = F)

methyn <- as.data.frame(methyn[rownames(methyn) %in% annot$IlmnID, ])
methyt <- as.data.frame(methyt[rownames(methyt) %in% annot$IlmnID, ])
methyn <- methyn[rownames(methyn) %in% rownames(methyt), ]
methyt <- methyt[rownames(methyt) %in% rownames(methyn), ]

transform.data <- function(df, stg, project) {
  #' Separate different stages, and melt for ggplot.
  #'
  #' Require reshape2 to work (uses melt function).
  #'  
  #' @param df The data frame whose rownames are to be converted.
  #' @param stg Stage information.
  #' @param project The TCGA project name.
  #'
  #' @return A melted data.frame with different stages.
  library(reshape2)
  df <- as.data.frame(df)
  df$IlmnID <- rownames(df)
  df <- melt(df, id.vars = "IlmnID")
  df$variable <- gsub("_", "-", df$variable)
  df$stage <- 0
  stage1 <- stg$barcode[grepl("(^i$)|(\\si[abc]?$)|(1)", stg$tumor_stage) &
                            stg$project == project]
  stage2 <- stg$barcode[grepl("(^ii$)|(\\si{2}[abc]?$)|(2)", stg$tumor_stage) &
                            stg$project == project]
  stage3 <- stg$barcode[grepl("(^iii$)|(\\si{3}[abc]?$)|(3)", stg$tumor_stage) &
                            stg$project == project]
  stage4 <- stg$barcode[grepl("(^iv$)|(\\siv[abc]?$)|(4)", stg$tumor_stage) &
                            stg$project == project]
  df$stage[df$variable %in% sub("-\\d{2}[ABC]$", "", stage1)] <- "i"
  df$stage[df$variable %in% sub("-\\d{2}[ABC]$", "", stage2)] <- "ii"
  df$stage[df$variable %in% sub("-\\d{2}[ABC]$", "", stage3)] <- "iii"
  df$stage[df$variable %in% sub("-\\d{2}[ABC]$", "", stage4)] <- "iv"
  df <- df[df$stage != 0, ]
  return(df)
}


dfn <- transform.data(methyn, stageAnnot, "BRCA")
dfn <- left_join(dfn, annot, by = "IlmnID")
dfn <- separate_rows(dfn, c("UCSC_RefGene_Name", "UCSC_RefGene_Group"), sep = ";")
dfn <- unique(dfn)
dfn <- dfn %>%
       mutate(gene = paste(UCSC_RefGene_Name, UCSC_RefGene_Group, sep = "_"),
              symbol = UCSC_RefGene_Name, group = UCSC_RefGene_Group) %>%
       select(gene, symbol, group, stage, value, variable)

dft <- transform.data(methyt, stageAnnot, "BRCA")
dft <- left_join(dft, annot, by = "IlmnID")
dft <- separate_rows(dft, c("UCSC_RefGene_Name", "UCSC_RefGene_Group"), sep = ";")
dft <- unique(dft)
dft <- dft %>%
       mutate(gene = paste(UCSC_RefGene_Name, UCSC_RefGene_Group, sep = "_"),
              symbol = UCSC_RefGene_Name, group = UCSC_RefGene_Group) %>%
       select(gene, symbol, group, stage, value, variable)

exp.boxplot <- function(datan, datat) {
    #' Convert ensembl to gene symbol.
    #'
    #' Require dplyr to work (uses inner_join function).
    #' Require ggplot2 & ggpubr to finish the boxplots.
    #'  
    #' @param datan The data frame for normal samples.
    #' @param datat The data frame for tumor samples.
    #'
    #' @return Genes used in the boxplot.

    library(dplyr)
    library(ggplot2)
    library(ggpubr)

    if (nrow(datan) != 0 & nrow(datat) != 0){
        datan$sample <- "control"
        datat$sample <- "tumor"
        df <- rbind.data.frame(datan, datat)
        colnames(df) <- c("gene", "symbol", "group", "stage",
                          "methylation", "case", "sample")
        df$stage[df$sample == "control"] <- "control"
        df$methylation <- log2(df$methylation + 0.1)
        df$gene <- factor(df$gene, levels=(unique(df$gene[order(nchar(df$gene),
                                                                      df$gene)])))
        genelist <- levels(df$gene)
        # Remove genes that ANOVA says there's no difference between groups
        for (j in seq_along(genelist)) {
            temp <- df[df$gene == genelist[j], ]
            temp <- summary(aov(methylation ~ stage, data=temp))
            if (temp[[1]]$`Pr(>F)`[1] > 0.01) {
                df <- df[df$gene != genelist[j], ]
            }
        }
        # symbol variable value stage sample
        # Calculate mean of each group
        df <- df %>% group_by(gene) %>% mutate(meanMethy=mean(methylation))
        # Each group is compared to the mean value
        p <- ggplot(df, aes(stage, methylation))+
                geom_boxplot()+
                facet_wrap(facets = "gene")+
                stat_compare_means(ref.group="control", label="p.signif")+
                geom_hline(aes(yintercept = meanMethy), linetype = 2)+
                ggtitle("Breast cancer sialic acid transferase genes' methylation level")
        return(p)
    }
}

p <- exp.boxplot(dfn, dft)
ggsave(p, filename = "./methylation/sialic_acid_methylation.tiff",
       device = "tiff", width = 16, height = 9, units = "in")