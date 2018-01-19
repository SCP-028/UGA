library(foreach)
library(doSNOW)
cl <- makeCluster(3)
registerDoSNOW(cl)

setwd("C:/Users/yz73026/Desktop/")
genelist <- c(
              "GPI",  # Phosphoglucose isomerase
              "PGM[1235]",  # Phosphoglucomutase
              "G6PD",  # Glucose-6-Phosphate Dehydrogenase
              "PGLS",  # 6-Phosphogluconolactonase
              "MPI",  # Mannose Phosphate Isomerase
              "TKT",  # Transketolase
              "TALDO1",  # Transaldolase 1
              "GFPT[12]",  # Glutamine-Fructose-6-Phosphate Transaminase
              "PFK[LPM]",  # Phosphofructokinase
              "PMM[12]",  # Phosphomannomutase
              "PRPS[12]",  # Phosphoribosyl Pyrophosphate Synthetase
              "UGP2",  # UDP-Glucose Pyrophosphorylase 2
              "GNPNAT1",  # Glucosamine-Phosphate N-Acetyltransferase 1
              "PGM3",  # Phosphoglucomutase 3
              "UAP1",  # UDP-N-Acetylglucosamine Pyrophosphorylase 1
              "NANP",  # N-Acetylneuraminic Acid Phosphatase
              "CMAS",  # CMP-Neu5Ac Synthetase
              "GALT",  # 	Galactose-1-Phosphate Uridylyltransferase
              "GMPP[AB]",  # Mannose-1-Phosphate Guanyltransferase
              "PGD",  # Phosphogluconate Dehydrogenase
              "RPIA",  # Ribose 5-Phosphate Isomerase A
              "GMDS",  # 	GDP-Mannose 4,6-Dehydratase
              "TSTA3",  # GDP-4-Keto-6-Deoxy-D-Mannose-3,5-Epimerase-4-Reductase
              "FPGT",  # Fucose-1-Phosphate Guanylyltransferase
              "GALE",  # 	UDP-Galactose-4-Epimerase / UDP-Glucose 4-Epimerase
              "UGDH",  # UDP-Glucose Dehydrogenase
              "UXS1",  # UDP-Glucuronate Decarboxylase 1
              "GNE",  # UDP-GlcNAc-2-Epimerase / N-Acetylmannosamine Kinase
              "GALK2",  # N-Acetylgalactosamine Kinase
              # CMP-KDN synthase
              "NANS",  # N-Acetylneuraminate-9-Phosphate Synthase
              "FUK",  # Fucokinase
              # N-Acetylglucosamine-1-phosphate uridyltransferase (GlmU), exclusive to prokaryotes
              "NAGK",  # N-Acetylglucosamine Kinase
              # Glucosamine Kinase
              "HK[123]",  # Hexokinase
              "GALK1",  # Galactokinase 1
              "CTPS[12]"  # CTP Synthase
            )
stageAnnot <- data.table::fread("C:/Users/yz73026/Desktop/expression_FPKM/annotation/annot.tsv",
                                data.table=F, stringsAsFactors=F, sep=",")
projects <- c("BRCA", "COAD", "LUAD", "LUSC", "BLCA", "LIHC", "HNSC", "KICH", "THCA", "PRAD")
filenames <- paste0("./expression_FPKM/", projects, ".RData")


name.convert <- function(df, ensemblAnnot, codingOnly=F,
                         regex='', genelist='', description=F,
                         libSize=F) {
    #' Convert ensembl to gene symbol.
    #'
    #' Require dplyr to work (uses inner_join function).
    #'
    #' @param df The data frame whose rownames are to be converted.
    #' @param ensemblAnnot Annotation data.frame from `retrieve.ensembAnnot`.
    #' @param codingOnly Keep only the protein-coding genes and ditch the rest.
    #' @param regex The regular expression for extracting genes of interest.
    #' @param genelist The gene symbols of genes of interest.
    #' @param description Whether to keep description & ensembl of genes.
    #' @param libSize Keep library size for DEG test or not.
    #'
    #' @return A data.frame df with its rownames converted.
    library(dplyr)
    # Check rownames of df
    if ("ensembl" %in% colnames(df)) {
        rownames(df) <- as.character(df$ensembl)
        df$ensembl <- NULL
    }
    ensembl <- sub("(.*)\\..*$", "\\1", rownames(df))
    library_size <- colSums(df)
    if (regex != '') {
        message("Regex found, filtering...")
        ensemblAnnot <- ensemblAnnot[grep(regex, ensemblAnnot$hgnc_symbol), ]
    }
    else if (all(genelist != '')) {
        message("Gene list found, filtering...")
        genelist <- paste0(sub("(.*)", "(^\\1$)", genelist), collapse="|")
        ensemblAnnot <- ensemblAnnot[grep(genelist, ensemblAnnot$hgnc_symbol), ]
    }
    else {
        message("No regex or gene list provided, returning all genes...")
    }
    df$ensembl_gene_id <- ensembl
    df <- inner_join(ensemblAnnot, df, by = "ensembl_gene_id")
    df <- df[order(rowSums(df[ ,grep("TCGA", colnames(df))]), decreasing=T), ]
    df <- df[!duplicated(df$hgnc_symbol), ]
    rownames(df) <- df$hgnc_symbol
    df$hgnc_symbol <- NULL
    if (libSize) {
        if (description) {
            return(list(exp=df, libSize=library_size))
        }
        else {
            return(list(exp=df[ ,grep("TCGA", colnames(df))],
                        libSize=library_size))
        }
    }
    else {
        if (description) {
            return(df)
        }
        else {
            return(df[ ,grep("TCGA", colnames(df))])
        }
    }
}


transform.data <- function(df, annot, project) {
    #' Separate different stages, and melt for ggplot.
    #'
    #' Require reshape2 to work (uses melt function).
    #'
    #' @param df The data frame whose rownames are to be converted.
    #' @param annot Stage information.
    #' @param project The TCGA project name.
    #'
    #' @return A melted data.frame with different stages.
    library(reshape2)
    df <- df[rowMeans(df) >= 1, ]
    df$symbol <- rownames(df)
    df <- melt(df, id.vars = "symbol")
    df$stage <- 0
    stage1 <- annot$barcode[grepl("(^i$)|(\\si[abc]?$)|(1)", annot$tumor_stage) &
                                annot$project == project]
    stage2 <- annot$barcode[grepl("(^ii$)|(\\si{2}[abc]?$)|(2)", annot$tumor_stage) &
                                annot$project == project]
    stage3 <- annot$barcode[grepl("(^iii$)|(\\si{3}[abc]?$)|(3)", annot$tumor_stage) &
                                annot$project == project]
    stage4 <- annot$barcode[grepl("(^iv$)|(\\siv[abc]?$)|(4)", annot$tumor_stage) &
                                annot$project == project]
    df$stage[df$variable %in% stage1] <- "i"
    df$stage[df$variable %in% stage2] <- "ii"
    df$stage[df$variable %in% stage3] <- "iii"
    df$stage[df$variable %in% stage4] <- "iv"
    df <- df[df$stage != 0, ]
    return(df)
}


exp.boxplot <- function(datan, datat, proj, annot) {
    #' Convert ensembl to gene symbol.
    #'
    #' Require dplyr to work (uses inner_join function).
    #' Require ggplot2 & ggpubr to finish the boxplots.
    #'
    #' @param datan The data frame for normal samples.
    #' @param datat The data frame for tumor samples.
    #' @param proj The current cancer project.
    #' @param annot Stage information.
    #'
    #' @return Genes used in the boxplot.

    library(dplyr)
    library(ggplot2)
    library(ggpubr)

    datan <- transform.data(datan, annot, proj)
    datat <- transform.data(datat, annot, proj)
    datan <- datan[datan$symbol %in% datat$symbol, ]
    datat <- datat[datat$symbol %in% datan$symbol, ]

    if (nrow(datan) != 0 & nrow(datat) != 0){
        datan$sample <- "control"
        datat$sample <- "tumor"
        df <- rbind.data.frame(datan, datat)
        colnames(df) <- c("symbol", "barcode", "log2FPKM", "stage", "sample")
        df$stage[df$sample == "control"] <- "control"
        df$log2FPKM <- log2(df$log2FPKM + 0.1)
        df$symbol <- factor(df$symbol, levels=(unique(df$symbol[order(nchar(df$symbol),
                                                                      df$symbol)])))
        genelist <- levels(df$symbol)
        # Remove genes that ANOVA says there's no difference between groups
        for (j in seq_along(genelist)) {
            temp <- df[df$symbol == genelist[j], ]
            temp <- summary(aov(log2FPKM ~ stage, data=temp))
            if (temp[[1]]$`Pr(>F)`[1] > 0.01) {
                df <- df[df$symbol != genelist[j], ]
            }
        }
        # symbol variable value stage sample
        # Calculate mean of each group
        df <- df %>% group_by(symbol) %>% mutate(meanFPKM=mean(log2FPKM))
        # Each group is compared to the mean value
        p <- ggplot(df, aes(stage, log2FPKM))+
                geom_boxplot()+
                facet_wrap(facets = "symbol")+
                stat_compare_means(ref.group="control", label="p.signif")+
                geom_hline(aes(yintercept = meanFPKM), linetype = 2)+
                ggtitle(proj)
        message(paste0(proj, " boxplot finished."))
        return(p)
    }
    else {
        message(paste0(proj, " ignored because none of the genes had a significant difference."))
    }
}

pb <- txtProgressBar(max = length(projects), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

result <- foreach(i = seq_along(filenames),
                  .packages=c("reshape2", "dplyr", "ggplot2", "ggpubr"),
                  .options.snow = opts) %dopar% {
    load(filenames[i])
    if(ncol(datan) >= 10 & ncol(datat) >= 10) {
        datan <- name.convert(datan, annot, genelist=genelist)
        datat <- name.convert(datat, annot, genelist=genelist)
        p <- exp.boxplot(datan, datat, projects[i], stageAnnot)
        ggsave(filename=paste0("./glycosylation/", projects[i], ".tiff"),
               plot=p, device="tiff", width=16, height=9, units="in", dpi=200)
    }
}
close(pb)
stopCluster(cl)
save(result, genelist, file="./glycosylation.RData")