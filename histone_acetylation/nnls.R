## Use the non-negative linear least squares method to calculate
# Initial settings
setwd("C:/Users/jzhou/Desktop/")
set.seed(1005)
source("https://raw.githubusercontent.com/ggrothendieck/gsubfn/master/R/list.R")
library(dplyr)

# Prepare to load and save data
file_name <- list.files("./expression/", pattern="RData$")
projects <- sub("(.*)\\.up.*$", "\\1", file_name)
results <- vector("list", length = (length(projects) * 2))
names(results) <- as.vector(outer(projects, c("_normal", "_tumor"), paste0))
names(results) <- names(results)[order(names(results))]


gen.regex <- function(gene_name, ...) {
#   Collapse string array into a regular expression
    gene_name <- c(gene_name, unlist(list(...), recursive=T))
    return(paste0(sub("(.*)", "(^\\1\\\\|)", gene_name), collapse="|"))
}


retrieve.gene <- function(df, regex) {
#   Retrieve genes using a regular expression from large data frames.
    df <- as.data.frame(df[grep(regex, rownames(df)), , drop=F])
    rownames(df) <- sub("(.*)\\|.*$", "\\1", rownames(df))
    return(df)
}


build.nnls <- function(gen.df, con.df) {
#   Turn input data frames into suitable ones for package nnls.
#   gen.df leaves one gene on the LHS, and the others are taken the opposite
#   number and moved to the RHS.
    library(nnls)
    if (nrow(gen.df >= 1)) {
        y <- as.numeric(as.character(gen.df[1, ]))
        x <- rbind(con.df, -1 * gen.df[-1, , drop=F])
        df <- nnls(t(x), y)
        df$gen.genes <- rownames(gen.df)
        df$con.genes <- rownames(con.df)
        names(df$x) <- rownames(x)
        return(df)
    }
    else {
        return(NA)
    }
}


# Constants
generate_name <- gen.regex("ACLY","ACSS2")
consume_name <- gen.regex("ACACA", "ACACB", "HMGCR", "ACOT1", "ACOT7", "ACOT11", "ACOT12", "ACOT13", "CREBBP", "EP300", "KAT2B", "KAT5", "NCOA1", "NCOA3")


for (i in seq_along(projects)) {
    load(paste0("./expression/", file_name[i]))
    # Get genes that are related to Ac-CoA
    gen.normal <- retrieve.gene(datan, generate_name)
    gen.tumor <- retrieve.gene(datat, generate_name)
    con.normal <- retrieve.gene(datan, consume_name)
    con.tumor <- retrieve.gene(datat, consume_name)
    # Build non-negative least square linear models
    results[[2*i-1]] <- build.nnls(gen.normal, con.normal)
    results[[2*i]] <- build.nnls(gen.tumor, con.tumor)
}

# Merge coefficients into one data frame
coef.df <- data.frame(symbol=character())
for (i in seq_along(results)) {
    df <- results[[i]]
    x <- as.data.frame(df$x)
    colnames(x) <- paste0(names(results)[i], "_", df$gen.genes[!df$gen.genes %in% names(df$x)])
    x$symbol <- rownames(x)
    coef.df <- full_join(coef.df, x, by = "symbol")
}

save(results, coef.df, file="./nnls_mitochondrial.RData")

