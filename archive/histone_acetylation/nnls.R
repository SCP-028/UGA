## Use the non-negative linear least squares method to calculate
# Initial settings
setwd("C:/Users/jzhou/Desktop/expression_FPKM")
set.seed(1005)
# source(https://raw.githubusercontent.com/ggrothendieck/gsubfn/master/R/list.R)
library(dplyr)

# Prepare to load and save data
load("../ensembl_annot.RData")
file_name <- list.files("./", pattern="csv$")
projectsNormal <- file_name[grep("normal", file_name)]
projectsNormal <- sub("(.*)_normal.*$", "\\1", as.character(projectsNormal))
projectsTumor <- file_name[grep("tumor", file_name)]
projectsTumor <- sub("(.*)_tumor.*$", "\\1", as.character(projectsTumor))
projects <- intersect(projectsNormal, projectsTumor)
projects <- projects[order(projects)]
coefResults <- vector("list", length = length(projects))
names(coefResults) <- projects
results <- vector("list", length = (length(projects) * 2))
names(results) <- as.vector(outer(projects, c("_normal", "_tumor"), paste0))
names(results) <- names(results)[order(names(results))]


filter.df <- function(df, annotation) {
    #' Keep the rows that are wanted, using its ensembl - symbol ids.
    #'
    #' @param df The data frame to filter.
    #' @param annotation The annotation data frame containing only the genes wanted.
    library(dplyr)
    df <- df[which(sub("(.*)\\..*$", "\\1", rownames(df)) %in% annotation$ensembl_gene_id), ]
    df$ensembl_gene_id <- sub("(.*)\\..*$", "\\1", rownames(df))
    df <- inner_join(annotation, df, by="ensembl_gene_id")
    rownames(df) <- df$hgnc_symbol
    df <- df[order(rownames(df)), ]
    return(df[ ,4:ncol(df)])
}


build.nnls <- function(gen.df, con.df) {
    #'  Turn input data frames into suitable ones for package nnls.
    #'
    #' @param gen.df Leaves one gene on the LHS, and the others are taken the opposite number and moved to the RHS.
  library(nnls)
  if (nrow(gen.df >= 1)) {
    y <- as.numeric(as.character(gen.df[1, ]))
    x <- rbind(con.df, -1 * gen.df[-1, , drop=F])
    df <- nnls(t(x), y)
    df$gen.genes <- rownames(gen.df)
    df$con.genes <- rownames(con.df)
    df$geneY <- rownames(gen.df)[1]
    names(df$x) <- rownames(x)
    return(df)
  }
  else {
    return(NA)
  }
}


coef.exp <- function(df, coefDF) {
    #' 
    #'
    #' @param df Expression data frame.
    #' @param coefDF Coefficients data frame.
    df <- df[order(rownames(df)), , drop=F]
    coefDF <- coefDF[rownames(coefDF) %in% rownames(df), , drop=F]
    coefDF <- coefDF[order(rownames(coefDF)), , drop=F]
    coefExp <- matrix(nrow=nrow(df), ncol = ncol(df), dimnames = dimnames(df))
    rownames(coefExp) <- paste0("coef_exp_", rownames(coefExp))
    for (i in seq(nrow(coefExp))) {
        coefExp[i, ] <- as.numeric(as.numeric(coefDF[i, ]) * df[i, ])
    }
    df <- rbind.data.frame(df, coefExp)
    return(df)
}


# Constants
generate_name <- c("ACSS1","PDHA1","PDHB","CPT1A","CPT1B","CPT2")
consume_name <- c("CS","ACOT2","ACOT9","ACOT10","ACOT13")
genelist <- c(generate_name, consume_name)
annotation <- annot[annot$hgnc_symbol %in% genelist, ]

for (i in seq_along(projects)) {
    # load(paste0("./expression/", file_name[i]))
    datan <- data.table::fread(paste0("./", projects[i], "_normal.csv"),
                               data.table=F, header=T, stringsAsFactors=F)
    if (!any(grepl("ENSG", rownames(datan)))) {
        rownames(datan) <- datan$ensembl
        datan$ensembl <- NULL
    }
    
    datat <- data.table::fread(paste0("./", projects[i], "_tumor.csv"),
                               data.table=F, header=T, stringsAsFactors=F)
    if (!any(grepl("ENSG", rownames(datat)))) {
        rownames(datat) <- datat$ensembl
        datat$ensembl <- NULL
    }
    if (ncol(datan) >= 5 & ncol(datat) >= 5) {
        # Get genes that are related to Ac-CoA
        datan <- filter.df(datan, annotation)
        datat <- filter.df(datat, annotation)
        df <- cbind(datan, datat)
        genDF <- df[rownames(df) %in% generate_name, ]
        conDF <- df[rownames(df) %in% consume_name, ]

        # Build non-negative least square linear models
        coefResults[[i]] <- build.nnls(genDF, conDF)
        # Extract coefficients
        coefDF <- cbind(names(coef(coefResults[[i]])), coef(coefResults[[i]]))
        coefDF <- rbind(coefDF, c(coefResults[[i]]$geneY, 1))
        rownames(coefDF) <- coefDF[ ,1]
        coefDF <- coefDF[ ,2, drop=F]
        results[[2*i-1]] <- rbind.data.frame(coef.exp(datan[rownames(datan) %in% generate_name, ], coefDF),
                                             coef.exp(datan[rownames(datan) %in% consume_name, ], coefDF))
        results[[2*i]] <- rbind.data.frame(coef.exp(datat[rownames(datat) %in% generate_name, ], coefDF),
                                             coef.exp(datat[rownames(datat) %in% consume_name, ], coefDF))
    }
}

# Merge coefficients into one data frame
# coef.df <- data.frame(symbol=character())
# for (i in seq_along(results)) {
#     df <- results[[i]]
#     x <- as.data.frame(df$x)
#     colnames(x) <- names(results)[i]
#     x$symbol <- rownames(x)
#     coef.df <- full_join(coef.df, x, by = "symbol")
# }

save(results, coefResults, file="./nnls_mitochondrial.RData")

