## Correlation analysis of UCPs and mitochondrial outer membrane ion transporters.
setwd("C:/Users/jzhou/Desktop/expression")
library(ggplot2)
library(reshape)
library(plyr)

file_name <- list.files("./", pattern="RData$")
projects <- sub("(.*)\\.up.*$", "\\1", file_name)
results <- vector("list", length(projects) * 2)
names(results) <- as.vector(outer(projects, c("_normal", "_tumor"), paste0))
names(results) <- names(results)[order(names(results))]


cor.mtest <- function(df1, df2, conf.level = 0.95){
    df1 <- as.matrix(df1)
    df2 <- as.matrix(df2)
    p.mat <- matrix(NA, nrow(df1), nrow(df2))
    dimnames(p.mat) <- list(rownames(df1), rownames(df2))
    # diag(p.mat) <- 0
    # diag(lowCI.mat) <- diag(uppCI.mat) <- 1
    for(i in seq(nrow(df1))){
        for(j in seq(nrow(df2))){
            p.mat[i,j] <- cor.test(df1[i, ], df2[j, ], conf.level = conf.level)$p.value
            # lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
            # uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
        }
    }
    return(p.mat)
}


exp.genes <- function(df) {
    x <- rowSums(df)
    x <- x[x >= summary(x)["3rd Qu."]]
    df <- df[rownames(df) %in% names(x), ]
    return(df)
}


ucp5.t.test <- function(datan, datat) {
    x <- datan[grep("SLC25A14", rownames(datan)), ]
    y <- datat[grep("SLC25A14", rownames(datat)), ]
    if (class(x) == "numeric") {  # Got hits, otherwise would be empty matrix
        z <- t.test(x, y, na.rm=T)
        return(z$p.value <= 0.01 & z$estimate[1] < z$estimate[2])
    }
    else {
        return(FALSE)
    }
}


cor.heatmap <- function(cor.df, p.df, cancer_type, tissue="normal", filepath) {
    library(ggplot2)
    library(reshape)
    library(plyr)
    cor.df.m <- melt(cor.df)
    p.df.m <- melt(p.df)
    df <- cbind(cor.df.m, p.df.m$value)
    colnames(df) <- c("ptdase", "cardio", "cor.value", "p.value")
    df$stars <- cut(df$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
    # cor.df.s <- ddply(cor.df.m, .(variable), transform, rescale = scale(value))

    base_size <- 9
    # color_palette <- c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac')
    ggplot(df, aes(ptdase, cardio, fill=cor.value)) +
    geom_tile() +
    scale_fill_gradient2(low="#1fddff", high="#ff4b1f") +
    # scale_fill_brewer(type="div", palette = 1) +
    # scale_fill_manual(values=color_palette) +
    geom_text(aes(label=stars), color="black", size=5) +
    theme_grey(base_size = base_size) +
    labs(x = "", y = "") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(axis.text.x = element_text(size = base_size * 0.8, angle = 330, hjust = 0)) +
    ggtitle(paste0(cancer_type, "_", tissue))

    ggsave(filename = filepath, width = 12, height = 9, units = "in")
    return(df)
}



# df <- matrix(nrow=length(file_name), ncol=3, dimnames=list(projects, c("UCP5", "1st", "3rd")))
# library(corrplot)


for (i in seq_along(projects)) {
    load(file_name[i])
    # ucp5_normal <- as.numeric(datan[grep("SLC25A14", rownames(datan)), , drop=F])
    # ucp5_tumor <- as.numeric(datat[grep("SLC25A14", rownames(datat)), , drop=F])
    # df[i,1] <- sum(datan[grep("SLC25A14", rownames(datan)), ]) / ncol(datan)
    # df[i,2:3] <- summary(rowSums(datan))[c("1st Qu.", "3rd Qu.")] / ncol(datan)
    # datan <- exp.genes(datan)
    # datat <- exp.genes(datat)
    if (ucp5.t.test(datan, datat)) {
        cardio_normal <- datan[grep("(^AGPAT)|(^\\|51099)|(^PGS1)|(^CRLS)", rownames(datan)), , drop = F]
        
        ptdase_normal <- datan[grep("(^CLPP\\|)|(^LONP1\\|)|(^THOP1\\|)|(^HTRA2\\|)|(^PMPCA\\|)|(^PMPCB\\|)|(^SPG7\\|)|CLPX\\|)|(^AFG3L2\\|)", rownames(datan)), , drop = F]
        cardio_tumor <- datat[grep("(^AGPAT)|(^\\|51099)|(^PGS1)|(^CRLS)", rownames(datat)), , drop = F]
        
        ptdase_tumor <- datat[grep("(^CLPP\\|)|(^LONP1\\|)|(^THOP1\\|)|(^HTRA2\\|)|(^PMPCA\\|)|(^PMPCB\\|)|(^SPG7\\|)|CLPX\\|)|(^AFG3L2\\|)", rownames(datat)), , drop = F]
        # ucp_normal <- datan[grep("(^SLC25A14)|(^SLC25A27)|(^UCP3)|(^UCP2)", rownames(datan)), , drop = F]
        # transporter_tumor <- datat[grep("(^SLC25A14)|(^SLC25A27)|(^UCP3)|(^UCP2)|(^AGPAT)|(^\\|51099)|(^PGS1)|(^CRLS)", rownames(datat)), , drop = F]
        rownames(cardio_normal) <- sub("(.*)\\|.*$", "\\1", rownames(cardio_normal)) 
        rownames(ptdase_normal) <- sub("(.*)\\|.*$", "\\1", rownames(ptdase_normal)) 
        rownames(cardio_tumor) <- sub("(.*)\\|.*$", "\\1", rownames(cardio_tumor)) 
        rownames(ptdase_tumor) <- sub("(.*)\\|.*$", "\\1", rownames(ptdase_tumor))

        dfn <- cor(t(ptdase_normal), t(cardio_normal))
        dft <- cor(t(ptdase_tumor), t(cardio_tumor))
        p.n <- cor.mtest(ptdase_normal, cardio_normal, 0.95)
        p.t <- cor.mtest(ptdase_tumor, cardio_tumor, 0.95)

        results[[2*i-1]] <- cor.heatmap(cor.df=dfn, p.df=p.n, cancer_type=projects[i],
                                        tissue="normal", filepath=paste0("../ptdase_cardio/", projects[i], "_", "normal.png"))
        results[[2*i]] <- cor.heatmap(cor.df=dft, p.df=p.t, cancer_type=projects[i],
                                      tissue="tumor", filepath=paste0("../ptdase_cardio/", projects[i], "_", "tumor.png"))

        # png(filename=paste0("../cardiolipin/", names(results)[2*i-1], ".png"), width = 1024, height = 1024, units = "px")
        # corrplot.mixed(dfn, lower="number", upper="square", order="alphabet", p.mat=resn[[1]], sig.level=0.05, tl.pos="lt")
        # dev.off()
        # png(filename=paste0("../cardiolipin/", names(results)[2*i], ".png"), width = 1024, height = 1024, units = "px")
        # corrplot.mixed(dft, lower="number", upper="square", order="alphabet", p.mat=rest[[1]], sig.level=0.05, tl.pos="lt")
        # dev.off()
    }
}
save(results, file="../ptdase_cardio/results.RData")