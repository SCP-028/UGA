setwd("C:/Users/jzhou/Desktop/")
set.seed(1005)
# source("https://raw.githubusercontent.com/ggrothendieck/gsubfn/master/R/list.R")
library(dplyr)
library(reshape2)
library(ggplot2)

# Prepare to load and save data
file_name <- list.files("./expression/", pattern="RData$", full.names = T)
projects <- sub("(.*)\\.up.*$", "\\1", file_name)
load("./nnls_mitochondrial.RData")  # coef.df

res.gene <- sub(".*_(.*)$", "\\1", colnames(coef.df)[-1])
colnames(coef.df)[-1] <- sub("(.*)_.*$", "\\1", colnames(coef.df)[-1])

gen.genes <- c()
con.genes <- c()

results.mean <- data.frame(symbol=c(gen.genes, con.genes))

for (i in seq_along(projects)) {
    load(file_name[i])
    # Normal #
    rownames(datan) <- sub("(.*)\\|.*$", "\\1", rownames(datan))
    mean.normal <- as.matrix(rowMeans(datan))

    coef.normal <- coef.df[ ,c(1, (i+1))]  # corresponding coefficients
    coef.normal <- coef.normal[rowSums(is.na(coef.normal)) == 0, ]  # remove NA
    coef.normal <- rbind.data.frame(coef.normal, c(res.gene[i], 1))  # add LHS gene
    coef.normal[ ,2] <- as.numeric(as.character(coef.normal[ ,2]))
    coef.normal[coef.normal$symbol %in% results[[i]]$gen.genes, 2] <- -1 * coef.normal[coef.normal$symbol %in% results[[i]]$gen.genes, 2]

    mean.normal <- mean.normal[rownames(mean.normal) %in% coef.normal$symbol, , drop=F]
    mean.normal <- mean.normal[order(rownames(mean.normal)), ]
    coef.normal <- coef.normal[order(coef.normal$symbol), 2]  # make sure the order is the same
    # mean * coef
    res.normal <- as.data.frame(mean.normal * coef.normal)
    res.normal[res.normal >0] <- res.normal[res.normal >0] / sum(res.normal[res.normal <0])
    res.normal[res.normal <0] <- res.normal[res.normal <0] / sum(res.normal[res.normal <0])  # must do in this order!
    colnames(res.normal) <- paste(names(results)[i], "normal", sep = "_")
    res.normal$symbol <- rownames(res.normal)
    results.mean <- full_join(results.mean, res.normal)
    # Tumor #
    rownames(datat) <- sub("(.*)\\|.*$", "\\1", rownames(datat))
    mean.tumor <- as.matrix(rowMeans(datat))

    coef.tumor <- coef.df[ ,c(1, (i+1))]  # corresponding coefficients
    coef.tumor <- coef.tumor[rowSums(is.na(coef.tumor)) == 0, ]  # remove NA
    coef.tumor <- rbind.data.frame(coef.tumor, c(res.gene[i], 1))  # add LHS gene
    coef.tumor[ ,2] <- as.numeric(as.character(coef.tumor[ ,2]))
    coef.tumor[coef.tumor$symbol %in% results[[i]]$gen.genes, 2] <- -1 * coef.tumor[coef.tumor$symbol %in% results[[i]]$gen.genes, 2]

    mean.tumor <- mean.tumor[rownames(mean.tumor) %in% coef.tumor$symbol, , drop=F]
    mean.tumor <- mean.tumor[order(rownames(mean.tumor)), ]
    coef.tumor <- coef.tumor[order(coef.tumor$symbol), 2]  # make sure the order is the same
    # mean * coef
    res.tumor <- as.data.frame(mean.tumor * coef.tumor)
    res.tumor[res.tumor > 0] <- res.tumor[res.tumor > 0] / sum(res.tumor[res.tumor > 0])
    res.tumor[res.tumor < 0] <- res.tumor[res.tumor < 0] / sum(res.tumor[res.tumor < 0])
    colnames(res.tumor) <- paste(names(results)[i], "tumor", sep = "_")
    res.tumor$symbol <- rownames(res.tumor)
    results.mean <- full_join(results.mean, res.tumor)
}

save(results.mean, file="./nnls_summary.RData")

df <- melt.data.frame(results.mean, id.vars = "symbol", variable_name = "project")
colnames(df) <- c("symbol", "project", "value")
df$symbol <- factor(df.m$symbol, levels = c(gen.genes, con.genes))

base_size <- 9
ggplot(df, aes(symbol, project, fill=value)) +
    geom_tile() +
    # scale_fill_brewer(type = "div", palette = 5) +  # Red & Blue
    scale_fill_gradient2(low="#1fddff", high="#ff4b1f", na.value = "gray95") +
    theme_grey(base_size = base_size) +
    labs(x = "", y = "") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(axis.text.x = element_text(size = base_size * 0.8, angle = 330, hjust = 0)) +
    ggtitle("ax/sum(ax) result for mitochondrial histone acetylation genes")