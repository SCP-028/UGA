library(NMF)
library(ggpubr)
library(glue)

setwd("C:/Users/yz73026/Desktop")
stageAnnot <- read.csv("expression_FPKM/annotation/annot.csv")
load("./feb2018_paper_boxplots/genelist.RData")  # a list of genes in 7 pathways


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
    df <- df[order(rowSums(df[ ,grep("(TCGA)|(TARGET)", colnames(df))]), decreasing=T), ]
    df <- df[!duplicated(df$hgnc_symbol), ]
    rownames(df) <- df$hgnc_symbol
    df$hgnc_symbol <- NULL
    if (libSize) {
        if (description) {
            return(list(exp=df, libSize=library_size))
        }
        else {
            return(list(exp=df[ ,grep("(TCGA)|(TARGET)", colnames(df))],
                        libSize=library_size))
        }
    }
    else {
        if (description) {
            return(df)
        }
        else {
            return(df[ ,grep("(TCGA)|(TARGET)", colnames(df))])
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
    df <- as.data.frame(df)
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


paste.vectors <- function(x, y) {
    df <- expand.grid(unique(x), unique(y))
    df <- df[order(df$Var1), ]
    df <- paste(df$Var1, df$Var2, sep="-")
    return(df)
}


prepare.result <- function(expdf, genelist, pathway, project, symbolAnnot, stageAnnot) {
    expdf <- name.convert(expdf, symbolAnnot, genelist=genelist)
    expdf <- t(basis(nmf(t(expdf), rank=1, seed=1005)))
    rownames(expdf) <- project
    expdf <- transform.data(expdf, stageAnnot, project)
    if (nrow(expdf) != 0) {
        expdf$pathway <- pathway
        return (expdf)
    }
    else {
        return (NA)
    }
}


box.plot <- function(df, project) {
  p <- ggboxplot(
    df, x="pathway", y="FPKM", fill="stage",
    rotate = T)
  # Aesthetic settings
  p <- ggpar(
    p, xlab="Pathways", ylab="log2(FPKM + 0.1)", title=project,
    legend="right", legend.title="Stages", palette="npg", ggtheme=theme_light(),
    font.tickslab = "bold")
  return (p)
}


filenames <- list.files("./expression_FPKM/", pattern=".RData$", full.names = T)
projects <- sub(".*/([A-Z-]+)\\.RData", "\\1", filenames)
resultn <- data.frame()
resultt <- data.frame()

for (i in seq_along(projects)) {
    load(filenames[i])
    if (ncol(datan) >= 10 & ncol(datat) >= 10) {
        for (j in seq_along(GENELIST)) {
            resultn <- rbind.data.frame(resultn,
                                        prepare.result(expdf=datan, genelist=GENELIST[[j]],
                                                       pathway=names(GENELIST)[j], project=projects[i],
                                                       symbolAnnot=annot, stageAnnot=stageAnnot)
                                        )
            resultt <- rbind.data.frame(resultn,
                                        prepare.result(expdf=datat, genelist=GENELIST[[j]],
                                                       pathway=names(GENELIST)[j], project=projects[i],
                                                       symbolAnnot=annot, stageAnnot=stageAnnot)
                                        )
        }
    }
}
resultn$stage <- "control"
res <- rbind.data.frame(resultn, resultt)
colnames(res) <- c("project", "sample", "FPKM", "stage", "pathway")
res <- res[rowSums(is.na(res)) == 0, ]
res$FPKM <- log2(res$FPKM + 0.1)
res$stage <- factor(res$stage, levels=c("control", "i", "ii", "iii", "iv"))
res$pathway <- factor(res$pathway, levels=names(GENELIST))
projects <- projects[projects %in% res$project]

for (i in seq_along(projects)) {
    df <- res[res$project == projects[i], ]
    p <- box.plot(df, projects[i])
    ggsave(filename=glue("{projects[i]}.tiff"),
           plot=p, device="tiff", width=6, height=10, units="in", dpi=300)
}
# sessionInfo()
