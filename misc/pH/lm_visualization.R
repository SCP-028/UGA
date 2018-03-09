library(reshape2)
library(ggplot2)

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


lm.p <- function (modelobject) {
    #' Return the p-value of a fitted model using F-test.
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}


lm.test <- function(formula, data) {
    #' Return statistics of a linear model.
    model <- lm(formula, data)
    result <- summary(model)
    result$p.value <- lm.p(model)
    return(result)
}


setwd("C:/Users/yz73026/Desktop")
GENELIST <- c(
    "GAPDH", "PYCR1",
    "PYCR2", "G6PD"
)
filenames <- list.files("./expression_FPKM/", pattern=".RData$", full.names = T)
projects <- sub(".*/([A-Z-]+)\\.RData", "\\1", filenames)

model1 <- vector("list", length=length(projects))
model2 <- vector("list", length=length(projects))

for (i in seq_along(projects)) {
    load(filenames[i])
    if (ncol(datat) >= 10) {
        datat <- name.convert(datat, annot, genelist=GENELIST)
        df <- as.data.frame(t(log2(datat + 0.1)))
        model1[[i]] <- lm.test(formula=GAPDH ~ PYCR1, data=df)
        model2[[i]] <- lm.test(formula=G6PD ~ PYCR1 + PYCR2, data=df)
    }
}

names(model1) <- names(model2) <- projects

coef1 <- sapply(model1, function(x) as.data.frame(x$coefficients)$Estimate[-1])
coef2 <- sapply(model2, function(x) as.data.frame(x$coefficients)$Estimate[-1])
rownames(coef2) <- c("PYCR1", "PYCR2")

df1 <- data.frame(symbol="GAPDH~PYCR1", project=names(coef1), value=(coef1))
df2 <- melt(coef2)
colnames(df2) <- c("symbol", "project", "value")
df2$symbol <- paste0("G6PD~",df2$symbol)
df <- rbind.data.frame(df1, df2)

p <- ggplot(df, aes(symbol, project, fill=value))+
        geom_tile(color="white")+
        #  ggtitle("GAPDH ~ PYCR1")+
        scale_fill_gradient2(low="#6D9EC1", high="#E46726", mid="white",
                            midpoint=0, limit=c(-0.5, 1.5), space="Lab",
                            name="Regression\ncoefficients")+
        coord_fixed()+
        theme_minimal()+
        theme(
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            panel.border=element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            text=element_text(family="HelveticaNeue"),
            axis.text.x=element_text(angle=45, hjust=1, vjust=1.0, size=10),
            # legend.justification=c(0, 1),
            # legend.position="right",
            # legend.direction="vertical"
            )+
        guides(fill=guide_colorbar(barwidth=1, barheight=7,
                                title.position="top",
                                title.hjust=0.5))
ggsave(filename="FigureS5.tiff",
           plot=p, device="tiff", width=3, height=10, units="in", dpi=300)

# sessionInfo()
# R version 3.4.3 (2017-11-30)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows >= 8 x64 (build 9200)

# Matrix products: default

# locale:
# [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252
# [3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C
# [5] LC_TIME=English_United States.1252

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base

# other attached packages:
# [1] dplyr_0.7.4          ggplot2_2.2.1        reshape2_1.4.3       RevoUtils_10.0.7
# [5] RevoUtilsMath_10.0.1

# loaded via a namespace (and not attached):
#  [1] Rcpp_0.12.15     assertthat_0.2.0 R6_2.2.2         grid_3.4.3       plyr_1.8.4
#  [6] gtable_0.2.0     magrittr_1.5     scales_0.5.0     pillar_1.2.0     rlang_0.2.0
# [11] stringi_1.1.6    lazyeval_0.2.1   bindrcpp_0.2     tools_3.4.3      stringr_1.3.0
# [16] glue_1.2.0       munsell_0.4.3    compiler_3.4.3   pkgconfig_2.0.1  colorspace_1.3-2
# [21] bindr_0.1        tibble_1.4.2
