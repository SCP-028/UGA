rlist <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
   args <- as.list(match.call())
   args <- args[-c(1:2,length(args))]
   length(value) <- length(args)
   for(i in seq(along=args)) {
     a <- args[[i]]
     if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
   }
   x
}

GSEA.EnrichmentScore2 <- function(gene.list, gene.set,
                                  t.stat, weighted.score.type) {
#gene.list is the location 

    N <- length(gene.list)
    Nh <- length(gene.set)
    Nm <-  N - Nh

    loc.vector <- vector(length=N, mode="numeric")
    peak.res.vector <- vector(length=Nh, mode="numeric")
    valley.res.vector <- vector(length=Nh, mode="numeric")
    tag.correl.vector <- vector(length=Nh, mode="numeric")
    tag.diff.vector <- vector(length=Nh, mode="numeric")
    tag.loc.vector <- vector(length=Nh, mode="numeric")

    loc.vector[gene.list] <- seq(1, N)
    tag.loc.vector <- loc.vector[gene.set]

    tag.loc.vector <- sort(tag.loc.vector, decreasing = F)

    correl.vector=t.stat
    # tag.correl.vector is the weight of each gene in
    # the gene set
    if (weighted.score.type == 0) {
        tag.correl.vector <- rep(1, Nh)
    } else if (weighted.score.type == 1) {
        tag.correl.vector <- correl.vector[tag.loc.vector]
        tag.correl.vector <- abs(tag.correl.vector)
    }

    norm.tag <- 1.0/sum(tag.correl.vector)
    tag.correl.vector <- tag.correl.vector * norm.tag
    norm.no.tag <- 1.0/Nm
    tag.diff.vector[1] <- (tag.loc.vector[1] - 1)
    tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] -
                             tag.loc.vector[1:(Nh - 1)] - 1
    tag.diff.vector <- tag.diff.vector * norm.no.tag
    peak.res.vector <- cumsum(tag.correl.vector -
                              tag.diff.vector)
    valley.res.vector <- peak.res.vector - tag.correl.vector
    max.ES <- max(peak.res.vector)
    min.ES <- min(valley.res.vector)
    ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES),
                 digits=5)
    return(list(ES = ES))
}


#load("/Users/shacao/errands/test_fenton_ribosome/pathways.list.GOMaco.CORUM.RCT.HGNC.Msigdb.RData")

calc.gsea.1<-function(t.stat, weighted.score.type){
    gene.labels=names(t.stat)
    Ng <- length(pathways)
    temp.size.G=sapply(pathways, function(x) length(x))
    max.size.G <- max(temp.size.G)
    gs <- matrix(rep("null", Ng*max.size.G),
                     nrow=Ng, ncol= max.size.G)
    for (i in 1:Ng) {
        gene.set.size <- length(pathways[[i]])
        gene.set.name <- names(pathways)[i]
        gene.set.tags=pathways[[i]]
        existing.set <- is.element(gene.set.tags, gene.labels)
        set.size <- length(existing.set[existing.set == T])
        gs[i,] <- c(gene.set.tags[existing.set],
                    rep("null", max.size.G - set.size))
    }
    rownames(gs) = names(pathways)
    gs.new = gs[rowSums(gs!="null")>5 &
                rowSums(gs!="null") <500, ]
    Ng = nrow(gs.new)

    rows=length(t.stat)
    nperm=1000
    scores=vector("numeric",Ng)
    pvals=vector("numeric",Ng)
    pval.indic=0  # Permutation test

    obs.gene.list2 <- order(t.stat,decreasing=TRUE)
    phi <- matrix(nrow = Ng, ncol = nperm)
    for (i in 1:Ng) {
        gene.set <- gs.new[i,gs.new[i,] != "null"]
        gene.set2 <- match(gene.set, gene.labels)
        if(pval.indic==1) {
            for (r in 1:nperm) {
                reshuffled.gene.labels <- sample(1:rows)
                GSEA.results <- GSEA.EnrichmentScore2(
                                gene.list=reshuffled.gene.labels,
                                gene.set=gene.set2,
                                t.stat=t.stat[reshuffled.gene.labels],
                                weighted.score.type=weighted.score.type)
                phi[i, r] <- GSEA.results$ES
            }
        }
        GSEA.results <- GSEA.EnrichmentScore2(
                        gene.list=obs.gene.list2,
                        gene.set=gene.set2,
                        t.stat=t.stat,
                        weighted.score.type=weighted.score.type)
        scores[i] <- GSEA.results$ES
        if(pval.indic==1) {
            pvals[i] = ifelse(GSEA.results$ES > 0,
                        sum(phi[i,] > GSEA.results$ES)/nperm,
                        -sum(phi[i,] < GSEA.results$ES)/nperm)
        }
        gc()
    }
    data.ES=cbind(scores, pvals)
    colnames(data.ES)=c("ES","pval")
    rownames(data.ES)=rownames(gs.new)
    #names(pvals)=rownames(gs.new)
    return(data.ES)
}

merge.result <- function(pathways, directory = "./"){
    library(dplyr)
    resultFiles <- list.files(directory, pattern = "\\.txt")
    fileNames <- sub("(.*)\\.txt", "\\1", resultFiles)
    resultFiles <- paste0("./", resultFiles)
    df <- data.frame(annotation = names(pathways))
    df$annotation <- as.character(df$annotation)
    df <- df[!duplicated(df$annotation), , drop = F]
    pb <- txtProgressBar(min = 0, max = length(resultFiles), style = 3, char = "#")
    for (i in seq_along(resultFiles)){
        temp <- data.table::fread(resultFiles[i], sep = " ",
                                header = F, data.table = F)
        temp$V3 <- NULL
        colnames(temp) <- c("annotation", "ES")
        temp$annotation <- as.character(temp$annotation)
        temp <- temp[order(temp$annotation, decreasing = T), ]
        temp <- temp[!duplicated(temp$annotation), ]
        df <- left_join(df, temp, by = "annotation")
        colnames(df)[ncol(df)] <- fileNames[i]
        setTxtProgressBar(pb, i)
        gc()
    }
    return(df)
}

extractR2 <- function(coef2, coef3, coef4, coef5, results,
                      coef.num, valley.r2.value, only.tf = T){
#   Extract %dev value from results dataframe, and corresponding gene symbols.
    library(dplyr)
    if (only.tf){
        coefList <- list(coef2=coef2, coef3=coef3, coef4=coef4, coef5=coef5)
        df <- coefList[[paste0("coef", coef.num)]]
        coef.num <- paste0("coef", coef.num, "_%dev")
        df <- df[ ,c("gene_symbol", coef.num, "TFmethy")]
        colnames(df) <- c("gene_symbol", "r_squared", "TFmethy")
        df <- df %>% 
                filter(!is.na(r_squared)) %>%
                filter(r_squared >= valley.r2.value) %>%
                filter(TFmethy >= 2) %>%
                select(gene_symbol, r_squared) %>%
                arrange(-r_squared)
    }
    else {
        coef.num <- paste0("coef", coef.num, "_%dev")
            df <- results[ ,c("gene_symbol", coef.num)]
            colnames(df) <- c("gene_symbol", "r_squared")
            df <- df %>% 
                    filter(!is.na(r_squared)) %>%
                    filter(r_squared >= valley.r2.value) %>%
                    arrange(-r_squared)
    }

    return(df)
}

getGenelist <- function(dirpath){
#   Get gene list from gene expression data frames.
    filename <- list.files(dirpath, pattern = "\\.RData", full.names = T)
    gene.list <- character()
    for (i in seq_along(filename)){
        load(filename[i])
        gene.list <- as.character(c(gene.list, 
                                    rownames(datan), rownames(datat)))
    }
    gene.list <- gene.list[!duplicated(gene.list)]
    return(gene.list)
}