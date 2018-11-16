## Calculates the expression of UCPs in normal vs. cancer tissues, and draw density plots.
setwd("C:/Users/jzhou/Desktop/")
source("https://raw.githubusercontent.com/ggrothendieck/gsubfn/master/R/list.R")
library(ggplot2)
library(reshape)

file_name <- list.files("./SLC25A14/", pattern="RData$")
projects <- sub("(.*)\\.up.*$", "\\1", file_name)

df <- matrix(nrow = length(projects), ncol = 5,
                dimnames = list(projects, c("exp.normal", "exp.tumor", "fold.change", "p.value", "sig")))

result <- vector("list", 5)
gene_name <- names(result) <- c("UCP1", "UCP2", "UCP3", "SLC25A27", "SLC25A14")
result <- lapply(result, function(x) return(df))
result <- result[order(names(result))]


prepare.plot <- function(df, tissue="normal") {
    x <- as.data.frame(rowSums(df) / ncol(df))
    x$tissue <- tissue
    x$gene <- rownames(x)
    df <- as.data.frame(df)
    df$gene <- rownames(df)
    df$tissue <- tissue
    df <- melt(df)
    return(list(df, x))
}


for (i in seq_along(projects)) {
    load(paste0("./SLC25A14/", file_name[i]))
    ucp.normal <- datan[grep("(^UCP1\\|)|(^UCP2\\|)|(^UCP3\\|)|(^SLC25A27\\|)|(^SLC25A14\\|)", rownames(datan)), , drop=F]
    ucp.tumor <- datat[grep("(^UCP1\\|)|(^UCP2\\|)|(^UCP3\\|)|(^SLC25A27\\|)|(^SLC25A14\\|)", rownames(datat)), , drop=F]
    rownames(ucp.normal) <- sub("(.*)\\|.*$", "\\1", rownames(ucp.normal))
    rownames(ucp.tumor) <- sub("(.*)\\|.*$", "\\1", rownames(ucp.tumor))
    ucp.normal <- ucp.normal[order(rownames(ucp.normal)), , drop=F]
    ucp.tumor <- ucp.tumor[order(rownames(ucp.tumor)), , drop=F]
    for (j in seq(nrow(ucp.normal))) {  # t-test of normal and tumor tissues
        x <- t.test(ucp.normal[j, ], ucp.tumor[j, ])
        result[[rownames(ucp.normal)[j]]][i, 1:3] <- c(x$estimate[1], x$estimate[2], (x$estimate[2] / x$estimate[1]))
        if(x$estimate[2] > x$estimate[1]) {
            result[[rownames(ucp.normal)[j]]][i,4] <- x$p.value
            result[[rownames(ucp.normal)[j]]][i,5] <- ifelse(x$p.value <= 0.05, "Y", "N")
        }
    }
    list[ucp.normal, mean.normal] <- prepare.plot(ucp.normal, "normal")
    list[ucp.tumor, mean.tumor] <- prepare.plot(ucp.tumor, "tumor")
    df <- rbind(ucp.normal, ucp.tumor)
    mean.df <- rbind(mean.normal, mean.tumor)
    colnames(mean.df) <- c("mean.value", "tissue", "gene")
    for (j in seq_along(gene_name)) {  # Plot for each gene in certain cancer types
        ggplot(df[df$gene == gene_name[j], ], aes(x=value, fill=tissue))+
        geom_density(alpha=.3)+
        geom_vline(data=mean.df[mean.df$gene == gene_name[j], ], aes(xintercept=mean.value, color=tissue), linetype="dashed", size=1)+
        ggtitle(paste0(projects[i], "_", gene_name[j]))
        ggsave(filename=paste0("./UCP/",projects[i], "_", gene_name[j], ".png"), width = 12, height = 9, units = "in")
    }
}
