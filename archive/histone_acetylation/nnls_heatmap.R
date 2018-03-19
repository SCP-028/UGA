generate_name <- c("ACSS1","PDHA1","PDHB","CPT1A","CPT1B","CPT2")
consume_name <- c("CS","ACOT2","ACOT9","ACOT10","ACOT13")

res.mean.gen <- results.mean[rownames(results.mean) %in% gen.gene, ]
res.mean.con <- results.mean[rownames(results.mean) %in% con.gene, ]


prepare.plot <- function(df) {
    library(reshape)
    df <- apply(df, 2, function(x) abs(x / sum(x, na.rm=T)))
    df <- df[ ,grep("tumor", colnames(df))] - df[ ,grep("normal", colnames(df))]
    colnames(df) <- sub("_tumor", "", colnames(df))
    
    df <- melt(df)
    colnames(df) <- c("symbol", "project", "value")
    return(df)
}


res.mean.gen <- prepare.plot(res.mean.gen)
res.mean.con <- prepare.plot(res.mean.con)

library(ggplot2)
base_size <- 9
ggplot(res.mean.con, aes(symbol, project, fill=value)) +
    geom_tile() +
    scale_fill_gradient2(low="#1fddff", high="#ff4b1f", na.value = "gray95") +
    theme_grey(base_size = base_size) +
    labs(x = "", y = "") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(axis.text.x = element_text(size = base_size * 0.8, angle = 330, hjust = 0)) +
    ggtitle("T-test for histone acetylation genes")

ggsave(filename = "./t_test_acetylation.png", width = 12, height = 9, units = "in")
