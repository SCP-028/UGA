# Calculate fold-change of all genes with symbols.
library(edgeR)
library(dplyr)
library(glue)
library(foreach)
library(doSNOW)
setwd("~/data/expression_count")


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
  library(edgeR)
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
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design=design)
  ## logFC logCPM PValue
  et <- exactTest(y)$table
  colnames(et) <- c("log2_fold_change", "log2_CPM", "p_value")
  et$symbol <- rownames(et)
  rownames(et) <- NULL
  return (et)
}

projects <- c("BLCA", "BRCA", "COAD", "ESCA", "HNSC", "KICH", "KIRC",
              "KIRP", "LIHC", "LUAD", "LUSC", "READ", "STAD", "THCA")
filenames <- glue("{projects}.RData")

pb <- txtProgressBar(max = length(projects), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

cl <- makeCluster(length(projects))
registerDoSNOW(cl)

DEGs <- foreach (
    i = seq_along(projects),
    .combine="rbind.data.frame",
    .packages=c("edgeR", "dplyr", "glue"),
    .options.snow = opts
    ) %dopar% {
    load(filenames[i])
    datan <- name.convert(datan, annot, libSize = T)
    datat <- name.convert(datat, annot, libSize = T)
    tryCatch({
        df <- edgeR.de.test(datan, datat, "normal", "tumor", sepLibSize = T)
        df$project <- projects[i]
        return (df)
    },
    warning=function(w) {
        return (data.frame())
    },
    error=function(e) {
        return (data.frame())
    })
}


close(pb)
stopCluster(cl)
save(DEGs, file="./DEGs.RData")

# sessionInfo()
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux 9 (stretch)

# Matrix products: default
# BLAS/LAPACK: /usr/lib/libopenblasp-r0.2.19.so

# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C              LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base

# other attached packages:
# [1] doSNOW_1.0.16   snow_0.4-2      iterators_1.0.8 foreach_1.4.3   glue_1.2.0      dplyr_0.7.4     edgeR_3.20.1    limma_3.34.1

# loaded via a namespace (and not attached):
#  [1] locfit_1.5-9.1   Rcpp_0.12.13     codetools_0.2-15 lattice_0.20-35  assertthat_0.2.0 grid_3.4.1       R6_2.2.2         magrittr_1.5     rlang_0.1.4      bindrcpp_0.2     tools_3.4.1      compiler_3.4.1   pkgconfig_2.0.1  bindr_0.1
# [15] tibble_1.3.4
