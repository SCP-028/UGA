library(foreach)
library(doSNOW)
cl <- makeCluster(3)
registerDoSNOW(cl)
setwd("C:/Users/jzhou/Desktop/")
source('C:/Users/jzhou/Git/UGA/metastasis/DEG_func.R')

annotation <- data.table::fread('./annotation/ensembl_symbol.tsv',
                                data.table=F, stringsAsFactors=F,
                                sep='\t', header=T)
dehydrogenase <- list(annotation[grep("dehydrogenase", annotation$description), "hgnc_symbol"])
filenames <- list.files("./expression_count", pattern="[A-Z\\-]\\.RData$", full.names=T)
projects <- sub(".*/([A-Z\\-]+)\\.RData$", "\\1", filenames)
# Find DEGs in all cancer types
pb <- txtProgressBar(max = length(projects), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
result <- foreach(i = seq_along(filenames), .combine="rbind.data.frame",
                  .packages=c("edgeR", "dplyr"), .options.snow = opts) %dopar% {
    load(filenames[i])
    if(ncol(datan) >= 10 & ncol(datat) >= 10) {
        dfn <- name.convert(datan, annotation, libSize=T)
        dft <- name.convert(datat, annotation, libSize=T)
        x <- edgeR.de.test(dfn, dft, "normal", "tumor", sepLibSize=T)
        x$project <- projects[i]
        x[x$p_value <= 0.05 & abs(x$log2_fold_change) >= 1, ]
    }
}
close(pb)
stopCluster(cl)
save(result, file="./differentially_expressed_genes.RData")
# GSEA
library(fgsea)
projects <- unique(result$project)
fgseaRes <- data.frame()
pb <- txtProgressBar(max = length(projects), style = 3)
for(i in seq_along(projects)) {
    df <- result[result$project == projects[i], ]
    stats <- df$log2_fold_change
    names(stats) <- df$symbol
    x <- fgsea(pathways=dehydrogenase, stats=stats, nperm=10000)
    fgseaRes <- rbind.data.frame(fgseaRes, as.data.frame(x))
    setTxtProgressBar(pb, i)
}
rownames(fgseaRes) <- projects
save(fgseaRes, file="./fgsea_result.RData")