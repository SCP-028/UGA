## Calculates the expression of UCPs in normal vs. cancer tissues, and draw density plots.
setwd("C:/Users/jzhou/Desktop/")
set.seed(1005)
source("https://raw.githubusercontent.com/ggrothendieck/gsubfn/master/R/list.R")
library(ggplot2)

file_name <- list.files("./expression/", pattern="RData$")
projects <- sub("(.*)\\.up.*$", "\\1", file_name)

result <- data.frame()
p.value.cutoff <- 0.05

gene_name <- c("MCT1","MCT4","SMCT1","SMCT4","HDAC\\d","HDAC10","HDAC11","SIRT\\d","CREBBP","EP300","KAT2B","KAT5","NCOA1","NCOA3","ACSS1", "ACSS2","ACSS1","HMGR","PDHA1","PDHB","CPT1A","CPT1B","CPT2","CS","ACOT12","ACOT2","ACOT9","ACOT10","ACOT13","ACOT1","ACOT4","ACOT7", "ACOT8","ACOT9" ,"ACOT11","ACOT13")
regex_gene_name <- paste0(sub("(.*)", "(^\\1\\\\|)", gene_name), collapse="|")


for (i in seq_along(projects)) {
    load(paste0("./expression/", file_name[i]))
    ucp.normal <- datan[grep(regex_gene_name, rownames(datan)), , drop=F]
    ucp.tumor <- datat[grep(regex_gene_name, rownames(datat)), , drop=F]
    rownames(ucp.normal) <- sub("(.*)\\|.*$", "\\1", rownames(ucp.normal))
    rownames(ucp.tumor) <- sub("(.*)\\|.*$", "\\1", rownames(ucp.tumor))
    temp <- matrix(nrow = nrow(ucp.normal), ncol = 6)
    colnames(temp) <- c("symbol", "project", "normal_exp", "tumor_exp", "p_value", "regulation")
    for (j in seq(nrow(ucp.normal))) {  # t-test of normal and tumor tissues
        x <- t.test(ucp.normal[j, ], ucp.tumor[j, ])
        temp[j,1:5] <- c(rownames(ucp.normal)[j], projects[i], x$estimate[1], x$estimate[2], x$p.value)
        temp[j,6] <- ifelse(x$p.value <= p.value.cutoff, (x$estimate[2] / x$estimate[1]), 1)
    }
    result <- rbind(result, as.data.frame(temp))
}

save(result, file="./t_test.RData")
df <- result[ ,c("symbol", "project", "regulation")]
df <- reshape(df, idvar = "symbol", timevar = "project", direction = "wide")
colnames(df) <- sub("regulation.", "", colnames(df))
ord <- hclust(dist(df[ ,-1], method = "euclidean"), method = "ward.D2")$order
df.m <- melt(df, id.vars = "symbol", variable_name = "project")
colnames(df.m) <- c("symbol", "project", "log2_fold_change")
df.m$symbol <- factor(df.m$symbol, levels = df$symbol[ord])
df.m$log2_fold_change <- log2(as.numeric(as.character(df.m$log2_fold_change)))


base_size <- 9
ggplot(df.m, aes(symbol, project, fill=log2_fold_change)) +
    geom_tile() +
    # scale_fill_brewer(type = "div", palette = 5) +  # Red & Blue
    scale_fill_gradient2(low="#1fddff", high="#ff4b1f", na.value = "gray95") +
    theme_grey(base_size = base_size) +
    labs(x = "", y = "") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(axis.text.x = element_text(size = base_size * 0.8, angle = 330, hjust = 0)) +
    ggtitle("T-test for histone acetylation genes")

ggsave(filename = "./t_test_acetylation.png", width = 12, height = 9, units = "in")
save(result, file="./t_test.RData")