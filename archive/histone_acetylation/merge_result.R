gene_name <- character()
for (i in seq_along(result.normal)) {
    gene_name <- union(gene_name, rownames(result.normal[[i]]))
    gene_name <- union(gene_name, rownames(result.tumor[[i]]))
}

# df <- matrix(nrow=length(gene_name), ncol=(length(result.normal) * 2 + 1), 
             # dimnames=list(gene_name, as.vector(outer(projects, c("_normal", "_tumor"), paste0))))
df <- data.frame(row.names = gene_name)
df$gene_name <- gene_name

library(dplyr)
for (i in seq_along(result.normal)) {
    x <- as.data.frame(result.normal[[i]])
    x$gene_name <- rownames(x)
    y <- as.data.frame(result.tumor[[i]])
    y$gene_name <- rownames(y)

    df <- left_join(df, x, by = "gene_name")
    df <- left_join(df, y, by = "gene_name")
}

z <- colnames(df)[-1]
temp <- as.vector(outer(projects, c("_normal", "_tumor"), paste0))
colnames(df) <- c("gene_name", paste0(temp, "_", sub("(.*)\\|.*$", "\\1", z)))