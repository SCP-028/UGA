setwd("/lustre1/yz73026/array/GSEA_TF/all")
filenames <- list.files("./", pattern = "\\.txt")
filenames <- sub("\\.txt", "", filenames)
filepaths <- list.files("./", pattern = "\\.txt", full.names = T)
all_GSEA <- vector("list", length(filenames))
for (i in seq_along(filepaths)){
    df <- read.table(filepaths[i])
    df <- df[order(df$ES, decreasing = T), ]
    x <- df$ES
    names(x) <- rownames(df)
    all_GSEA[[i]] <- x
}
names(all_GSEA) <- filenames

setwd("/lustre1/yz73026/array/GSEA_TF/TF_only")
filenames <- list.files("./", pattern = "\\.txt")
filenames <- sub("\\.txt", "", filenames)
filepaths <- list.files("./", pattern = "\\.txt", full.names = T)
tf_GSEA <- vector("list", length(filenames))
for (i in seq_along(filepaths)){
    df <- read.table(filepaths[i])
    df <- df[order(df$ES, decreasing = T), ]
    x <- df$ES
    names(x) <- rownames(df)
    all_GSEA[[i]] <- x
}
names(tf_GSEA) <- filenames

save(all_GSEA, tf_GSEA, file = "../results.RData")