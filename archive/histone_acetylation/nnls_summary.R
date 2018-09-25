setwd("C:/Users/jzhou/Desktop/")
set.seed(1005)
# source("https://raw.githubusercontent.com/ggrothendieck/gsubfn/master/R/list.R")
library(dplyr)

# Prepare to load and save data
file_name <- list.files("./expression/", pattern="RData$", full.names = T)
projects <- sub("(.*)\\.up.*$", "\\1", file_name)
load("./nnls_mitochondrial.RData")

res.gene <- sub(".*_(.*)$", "\\1", colnames(coef.df)[-1])
colnames(coef.df)[-1] <- sub("(.*)_.*$", "\\1", colnames(coef.df)[-1])

results.mean <- data.frame(symbol=character())
results.median <- data.frame(symbol=character())

for (i in seq_along(projects)) {
    load(file_name[i])
    # Normal #
    rownames(datan) <- sub("(.*)\\|.*$", "\\1", rownames(datan))
    mean.normal <- as.matrix(rowMeans(datan))
    median.normal <- apply(datan, 1, median)

    coef.normal <- coef.df[ ,c(1, (2*i))]  # corresponding coefficients
    coef.normal <- coef.normal[rowSums(is.na(coef.normal)) == 0, ]  # remove NA
    coef.normal <- rbind.data.frame(coef.normal, c(res.gene[2*i-1], 1))  # add LHS gene
    coef.normal[ ,2] <- as.numeric(coef.normal[ ,2])
    coef.normal[coef.normal$symbol %in% results[[2*i-1]]$gen.genes, 2] <- -1 * coef.normal[coef.normal$symbol %in% results[[2*i-1]]$gen.genes, 2]

    mean.normal <- mean.normal[rownames(mean.normal) %in% coef.normal$symbol, , drop=F]
    mean.normal <- mean.normal[order(rownames(mean.normal)), ]
    median.normal <- median.normal[names(median.normal) %in% coef.normal$symbol]
    median.normal <- median.normal[order(names(median.normal))]

    coef.normal <- coef.normal[order(coef.normal$symbol), 2]  # make sure the order is the same
    # mean * coef
    res.normal <- as.data.frame(mean.normal * coef.normal)
    colnames(res.normal) <- names(results)[2*i-1]
    res.normal$symbol <- rownames(res.normal)
    results.mean <- full_join(results.mean, res.normal)
    # median * coef
    res.normal <- as.data.frame(median.normal * coef.normal)
    colnames(res.normal) <- names(results)[2*i-1]
    res.normal$symbol <- rownames(res.normal)
    results.median <- full_join(results.median, res.normal)
    # Tumor #
    rownames(datat) <- sub("(.*)\\|.*$", "\\1", rownames(datat))
    mean.tumor <- as.matrix(rowMeans(datat))
    median.tumor <- apply(datat, 1, median)

    coef.tumor <- coef.df[ ,c(1, (2*i+1))]  # corresponding coefficients
    coef.tumor <- coef.tumor[rowSums(is.na(coef.tumor)) == 0, ]  # remove NA
    coef.tumor <- rbind.data.frame(coef.tumor, c(res.gene[2*i], 1))  # add LHS gene
    coef.tumor[ ,2] <- as.numeric(coef.tumor[ ,2])
    coef.tumor[coef.tumor$symbol %in% results[[2*i]]$gen.genes, 2] <- -1 * coef.tumor[coef.tumor$symbol %in% results[[2*i]]$gen.genes, 2]

    mean.tumor <- mean.tumor[rownames(mean.tumor) %in% coef.tumor$symbol, , drop=F]
    mean.tumor <- mean.tumor[order(rownames(mean.tumor)), ]
    median.tumor <- median.tumor[names(median.tumor) %in% coef.tumor$symbol]
    median.tumor <- median.tumor[order(names(median.tumor))]

    coef.tumor <- coef.tumor[order(coef.tumor$symbol), 2]  # make sure the order is the same
    # mean * coef
    res.tumor <- as.data.frame(mean.tumor * coef.tumor)
    colnames(res.tumor) <- names(results)[2*i]
    res.tumor$symbol <- rownames(res.tumor)
    results.mean <- full_join(results.mean, res.tumor)
    # median * coef
    res.tumor <- as.data.frame(median.tumor * coef.tumor)
    colnames(res.tumor) <- names(results)[2*i]
    res.tumor$symbol <- rownames(res.tumor)
    results.median <- full_join(results.median, res.tumor)
}

save(results.mean, results.median, file="./nnls_summary.RData")