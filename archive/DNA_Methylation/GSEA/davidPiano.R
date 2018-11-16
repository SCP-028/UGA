calc.david <- function(genes,universe){
#        load("/Users/shacao/pathways.list.GOMaco.CORUM.RCT.HGNC.Msigdb.gsc.RData")
        library(piano)
#        gsc = loadGSC("/Users/shacao/c2.c5.c6.gmt.txt", type = "gmt")
        member <- genes
        p.hyper <- runGSAhyper(genes = member, rep(0, length(member)), 0.05,universe = universe , gsc, gsSizeLim = c(10, Inf), adjMethod = "fdr")
        return(cbind(p.hyper$pvalues, p.hyper$p.adj))
}

grep.df <- function(keyword, df, genes, gsc){
    temp <- df[grep(keyword, as.character(rownames(df)), ignore.case = T), ]
#   temp <- temp[as.numeric(temp[ ,1]) <= 0.05, ]
    x <- gsc[[1]]
    x <- x[names(x) %in% as.character(rownames(temp))]
    y <- list()
    for (i in seq_along(x)){
        z <- unlist(x[i])
        y[[i]] <- z[z %in% genes]
    }
#   Return Geneset and genesInGeneset
    return(list(x, y))
}
