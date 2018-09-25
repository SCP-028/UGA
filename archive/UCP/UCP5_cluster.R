## UCP5 (SLC25A14) clustering with all other genes in all cancer types
setwd("C:/Users/jzhou/Desktop/expression")
library(ggplot2)
library(reshape)
library(dplyr)
library(plyr)

# Read file names
file_name <- list.files("./", pattern="RData$")
projects <- sub("(.*)\\.up.*$", "\\1", file_name)

# Initiate empty data frames used for full join later
df.normal <- df.tumor <- data.frame(entrez=character())

# Constants
base_size <- 9  # for ggplot


ucp5.t.test <- function(datan, datat) {
#   Only keep cancer types that actually had hits of SLC25A14
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


cor.mtest <- function(df1, df2, conf.level = 0.95){
#   Get p-value of correlation tests
    df1 <- as.matrix(df1)
    df2 <- as.matrix(df2)
    p.mat <- matrix(NA, nrow(df1), nrow(df2))
    dimnames(p.mat) <- list(rownames(df1), rownames(df2))
    for(i in seq(nrow(df1))){
        for(j in seq(nrow(df2))){
            p.mat[i,j] <- cor.test(df1[i, ], df2[j, ], conf.level = conf.level)$p.value
        }
    }
    return(p.mat)
}


cluster.exp <- function(df.all, df.tissue, num=150) {
#   Cluster and plot `num` genes.
    df.tissue <- df.tissue[order(abs(df.tissue[ ,1]), decreasing=T), ]
    df.tissue <- df.tissue[1:min(num, nrow(df.tissue)), ]
    # Put SLC25A14 in
    gene_list <- c(df.tissue$entrez, "9016")
    rownames(df.all) <- sub(".*\\|(.*)$", "\\1", rownames(df.all))
    df.all <- as.data.frame(df.all[rownames(df.all) %in% gene_list, , drop=F])
    rownames(df.all)[rownames(df.all) == "9016"] <- "SLC25A14"

    # Cluster genes against different cases
    ord <- hclust(dist(df.all, method="euclidean"), method="ward.D2")$order
    df.all$entrez <- rownames(df.all)
    df.m <- melt(df.all, id.vars="entrez", variable_name="case")
    colnames(df.m)[3] <- "expression"
    df.m$entrez <- factor(df.m$entrez, levels=df.all$entrez[ord])
    return(df.m)
}


# Calculate correlation value and p-value for SLC25A14 and all other genes
for (i in seq_along(projects)) {
    load(file_name[i])
    if (ucp5.t.test(datan, datat)) {
        ucp5_normal <- datan[grep("(^SLC25A14\\|)", rownames(datan)), , drop = F]
        other_normal <- datan[rownames(datan) != rownames(ucp5_normal), , drop = F]
        ucp5_tumor <- datat[grep("(^SLC25A14\\|)", rownames(datat)), , drop = F]
        other_tumor <- datat[rownames(datat) != rownames(ucp5_normal), , drop = F]

        rownames(ucp5_normal) <- sub(".*\\|(.*)$", "\\1", rownames(ucp5_normal)) 
        rownames(other_normal) <- sub(".*\\|(.*)$", "\\1", rownames(other_normal)) 
        rownames(ucp5_tumor) <- sub(".*\\|(.*)$", "\\1", rownames(ucp5_tumor)) 
        rownames(other_tumor) <- sub(".*\\|(.*)$", "\\1", rownames(other_tumor))

        dfn <- cor(t(ucp5_normal), t(other_normal))
        dft <- cor(t(ucp5_tumor), t(other_tumor))
        p.n <- cor.mtest(ucp5_normal, other_normal, 0.95)
        p.t <- cor.mtest(ucp5_tumor, other_tumor, 0.95)
        
        # Keep only genes that are significantly correlated (p <= 0.05 & cor >= 0.5) 
        dfn <- as.data.frame(t(dfn[ ,which(p.n <= 0.05 & abs(dfn) >= 0.5), drop=F]))
        dft <- as.data.frame(t(dft[ ,which(p.t <= 0.05 & abs(dft) >= 0.5), drop=F]))
        colnames(dfn) <- colnames(dft) <- projects[i]
        dfn$entrez <- rownames(dfn)
        dft$entrez <- rownames(dft)
        # p.n <- p.n[ ,colnames(p.n) %in% colnames(dfn), drop=F]
        # p.t <- p.t[ ,colnames(p.t) %in% colnames(dft), drop=F]

        # Merge with the big data frame, keeping all rows
        df.normal <- full_join(df.normal, dfn, by="entrez")
        df.tumor <- full_join(df.tumor, dft, by="entrez")

        # Cluster
        dfn.m <- cluster.exp(datan, dfn, num=150)
        dft.m <- cluster.exp(datat, dft, num=150)

        # Heatmap of these genes' expression
        ggplot(dfn.m, aes(entrez, case, fill=expression)) +
        geom_tile() +
        scale_fill_gradient2(low="#1fddff", high="#ff4b1f",
                             na.value="gray95") +
        theme_grey(base_size=base_size) +
        labs(x="", y="") +
        scale_x_discrete(expand=c(0, 0)) +
        scale_y_discrete(expand=c(0, 0)) +
        theme(axis.text.x=element_text(size=base_size * 0.8,
                                       angle=270, hjust=0),
               axis.text.y=element_blank(),
               axis.ticks.y=element_blank()) +
        ggtitle(paste0(projects[i], " normal"))
        ggsave(filename=paste0("../SLC25A14/clustering/",
                               projects[i], "_normal.png"),
               width=12, height=9, units="in")

        ggplot(dft.m, aes(entrez, case, fill=expression)) +
        geom_tile() +
        scale_fill_gradient2(low="#1fddff", high="#ff4b1f",
                             na.value="gray95") +
        theme_grey(base_size=base_size) +
        labs(x="", y="") +
        scale_x_discrete(expand=c(0, 0)) +
        scale_y_discrete(expand=c(0, 0)) +
        theme(axis.text.x=element_text(size=base_size * 0.8,
                                       angle=270, hjust=0),
               axis.text.y=element_blank(),
               axis.ticks.y=element_blank()) +
        ggtitle(paste0(projects[i], " tumor"))
        ggsave(filename=paste0("../SLC25A14/clustering/",
                               projects[i], "_tumor.png"),
               width=12, height=9, units="in")

        print(paste0("Finished: ", projects[i]))
    }
}

# Remove entrez column
rownames(df.normal) <- df.normal$entrez
rownames(df.tumor) <- df.tumor$entrez
df.normal$entrez <- NULL
df.tumor$entrez <- NULL

save(df.normal, df.tumor, file="../SLC25A14/clustering/results.RData")