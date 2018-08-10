retrieve.ensemblAnnot <- function(df, codingOnly=F) {
    #' Prepare for converting ensembl to gene symbol.
    #'
    #' Require biomaRt to work (retrieve annotation).
    #'
    #' @param df The data frame whose rownames are to be converted.
    #' @param codingOnly Keep only the protein-coding genes and ditch the rest.
    #'
    #' @return A data.frame of matching ensembl IDs and symbols.
    require(biomaRt)

    if ("ensembl" %in% colnames(df)) {
        rownames(df) <- as.character(df$ensembl)
    }
    # Retrieve biomaRt annotation data frame
    ensembl <- sub("(.*)\\..*$", "\\1", rownames(df))
    message("Downloading biomaRt manifest...")
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    message("Retrieving annotation data...")
    ensemblAnnot <- getBM(filters= "ensembl_gene_id",
                          attributes= c("ensembl_gene_id", "hgnc_symbol",
                                        "entrezgene", "gene_biotype", "description", "definition_1006"),
                          values=ensembl, mart= mart)
    ensemblAnnot <- ensemblAnnot[ensemblAnnot$hgnc_symbol != '', ]
    if (codingOnly) {
        ensemblAnnot <- ensemblAnnot[ensemblAnnot$gene_biotype == 'protein_coding', ]
        ensemblAnnot$gene_biotype <- NULL
    }
    return(ensemblAnnot)
}


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
    require(dplyr)
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


transform.data <- function(df, stageAnnot, project) {
  #' Separate different stages, and melt for ggplot.
  #'
  #' Require tidyr to work (uses the gather function).
  #'
  #' @param df The data frame whose rownames are to be converted.
  #' @param stageAnnot Stage information.
  #' @param project The TCGA project name.
  #'
  #' @return A melted data.frame with different stages.

  require(tidyr)
  # df <- df[rowMeans(df) >= 1, ]
  df$variable <- rownames(df)
  df <- df %>% gather(key="variable")
  df$stage <- 0
  stage1 <- stageAnnot$barcode[grepl("(^i$)|(\\si[abc]?$)|(1)", stageAnnot$tumor_stage) &
                                 stageAnnot$project == project]
  stage2 <- stageAnnot$barcode[grepl("(^ii$)|(\\si{2}[abc]?$)|(2)", stageAnnot$tumor_stage) &
                                 stageAnnot$project == project]
  stage3 <- stageAnnot$barcode[grepl("(^iii$)|(\\si{3}[abc]?$)|(3)", stageAnnot$tumor_stage) &
                                 stageAnnot$project == project]
  stage4 <- stageAnnot$barcode[grepl("(^iv$)|(\\siv[abc]?$)|(4)", stageAnnot$tumor_stage) &
                                 stageAnnot$project == project]
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

    require(dplyr)
    require(ggplot2)
    require(ggpubr)

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


prep.cor.heatmap <- function(df, genelist, cluster=F,
                             sample=c("normal", "tumor")) {
    #' Calculate spearman correlation coefficients and cluster.
    #'
    #' @param df The input values.
    #' @param genelist Used for arranging the order of the genes.
    #' @param cluster Whether to cluster the genes according to coefficients.
    #' @param sample Normal or tumor samples, used for later plotting.
    #'
    #' @return A melted data.frame of the correlation coefficients.
    require(reshape2)
    df <- df[match(genelist, rownames(df)), ]
    df <- df[apply(df, 1, var) != 0, ]
    df <- cor(t(df), method="spearman")
    if (cluster) {
        dd <- as.dist((1 - df) / 2)
        hc <- hclust(dd)
        df <- df[hc$order, hc$order]
    }
    df[upper.tri(df)] <- NA
    df <- melt(df, na.rm=T)
    df$sample <- sample
    colnames(df) <- c("Var1", "Var2", "value", "sample")
    return(df)
}


nice.heatmap <- function(df, title) {
    #' Plot a nice-looking heatmap.
    #'
    #' @param df The melted data.frame used by ggplot.
    #' @param title The title of the plot.
    #'
    #' @return A ggplot list.
    require(ggplot2)
    ggplt <- ggplot(df, aes(Var1, Var2, fill=value))+
             geom_tile(color="white")+
             ggtitle(title)+
             scale_fill_gradient2(low="#6D9EC1", high="#E46726", mid="white",
                                  midpoint=0, limit=c(-1, 1), space="Lab",
                                  name="Spearman\nCorrelation")+
             coord_fixed()+
             theme_minimal()+
             theme(
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   panel.border=element_blank(),
                   text=element_text(family="HelveticaNeue"),
                   axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
                   legend.justification=c(0, 1),
                   legend.position=c(0, 1),
                   legend.direction="horizontal")+
             guides(fill=guide_colorbar(barwidth=7, barheight=1,
                                        title.position="top",
                                        title.hjust=0.5))
    return(ggplt)
}


# Or just use `ggarrange(p1, p2, col=2)` in `ggpubr`
multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL) {
    #' http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
    #'
    #' @param ... ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
    #' @param cols Number of columns in layout
    #' @param layout A matrix specifying the layout. If present, 'cols' is ignored.
    #'
    #' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
    #' then plot 1 will go in the upper left, 2 will go in the upper right, and
    #' 3 will go all the way across the bottom.
    #'
    #' @examples
    #' multiplot(plot1, plot2, cols=2)
    require(grid)
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
    }
    if (numPlots==1) {
        print(plots[[1]])
    }
    else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout),
                                                   ncol(layout))))
        # Make each plot, in the correct location
        for (i in seq(numPlots)) {
          # Get the matrix positions of the regions that contain this subplot
          matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
          print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                          layout.pos.col = matchidx$col))
        }
    }
}


edgeR.de.test <- function(df1, df2, group1, group2, sepLibSize=F) {
  #' Use edgeR to perform gene differential expression analysis.
  #'
  #' Require edgeR and dplyr to work.
  #'
  #' @param df1 First data frame / matrix. Must be counts value!!
  #' @param df2 Second. The result comes as df2:df1.
  #' @param group1 Name of first group.
  #' @param group2 Name of second group.
  #' @param sepLibSize If `name.convert` has a separate library size,
  #'                   set this to TRUE.
  require(edgeR)
  if(sepLibSize) {
    libSize <- c(df1$libSize, df2$libSize)
    df1 <- df1$exp[order(rownames(df1$exp)), ]
    df2 <- df2$exp[order(rownames(df2$exp)), ]
    group <- c(rep(group1, ncol(df1)), rep(group2, ncol(df2)))
    design <- model.matrix(~0+group)
    x <- cbind(df1, df2)
    y <- DGEList(counts=x, group=group, lib.size=libSize)
  }
  else {
    df1 <- df1[order(rownames(df1)), ]
    df2 <- df2[order(rownames(df2)), ]
    group <- c(rep(group1, ncol(df1)), rep(group2, ncol(df2)))
    design <- model.matrix(~0+group)
    x <- cbind(df1, df2)
    y <- DGEList(counts=x, group=group)
  }
  y <- y[rowSums(cpm(y) >= 1) >= 1, , keep.lib.sizes=F]
  if (nrow(y$counts) > 1) {
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design=design)
    ## logFC logCPM PValue
    et <- exactTest(y)$table
    colnames(et) <- c("log2_fold_change", "log2_CPM", "p_value")
    et$symbol <- rownames(et)
    rownames(et) <- NULL
    return (et)
  }
  else {
    return (NA)
  }
}


transform.for.plot <- function(df, stage='stage_i_iii') {
    #' Melt data frame for ggplot.
    #'
    #' @param df The data frame prepared to melt.
    #' @param stage The tumor stage of `df`.
    library(reshape2)
    df$symbol <- rownames(df)
    df <- melt(df, id.vars = 'symbol', variable.name = 'sample', value.name = 'log2_FPKM')
    df$log2_FPKM <- log2(as.numeric(as.character(df$log2_FPKM)) + 0.1)
    df$stage <- stage
    return(df)
}


get.DEG <- function(df, result, p_value=0.01, logFC=1, logCPM=1,
                    FPKM=0, up=F, strict=F) {
    #' Filter the edgeR results to get only DEGs.
    #'
    #' @param df [data.frame] The FPKM values .
    #' @param result [data.frame] The result given by edgeR.
    #' @param p_value [float] The p-value cutoff.
    #' @param logFC [int] the log2(Fold Change) cutoff.
    #' @param logCPM [int] the log2(Count Per Million) cutoff.
    #' @param FPKM [int] the filtering value for df, useless
    #' if strict=F.
    #' @param up [bool] Whether to keep the down-regulated
    #' genes.
    #' @param strict [bool] Filtering rows in df with too many
    #' zeros.
    #' @return [data.frame] a filtered df.
    if (strict) {
        df <- df[(rowSums(df <= FPKM) <= (ncol(df) / 2)), ]
    }
    result <- result[(abs(result$log2_fold_change) >= logFC) &
                     (result$p_value <= p_value) &
                     (result$log2_CPM >= logCPM), ]
    if (up) {
        result <- result[result$log2_fold_change >= logFC, ]
    }
    df <- df[rownames(df) %in% result$symbol, ]
    df$hgnc_symbol <- rownames(df)
    return(df)
}
