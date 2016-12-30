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

calc.gsea.1<-function(t.stat,weighted.score.type){
    t.stat=t.stat[!is.na(t.stat)]
    gene.labels=names(t.stat)

    #pathways=sapply(pathways, function(x)x[-c(1:2)])
    #pathways=pathways.list.GOMaco.CORUM.RCT.HGNC.Msigdb
    Ng <- length(pathways)
    temp.size.G=sapply(pathways, function(x)length(x))
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
    nperm=5000
    scores=vector("numeric",Ng)
    pvals=vector("numeric",Ng)
    pval.indic=0

    obs.gene.list2 <- order(t.stat,decreasing=TRUE)
    phi <- matrix(nrow = Ng, ncol = nperm)
    for (i in 1:Ng) {
        if(i %%100==1) {
            print(i)
        }
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
