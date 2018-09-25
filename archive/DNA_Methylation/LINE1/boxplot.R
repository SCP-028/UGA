gene <- read.table("./methyLevel.txt", sep = "\t", header = T)
drawBox <- function(df, group, filename){
    df <- df[df$Group == group, ]
    df$Group <- NULL
    df$mCG <- NULL
    df$Type <- paste(sub("TCGA-(.*)", "\\1", df$Type.1), 
                            df$Type, sep = ".")
    df$Type <- sub("Normal", "normal", df$Type)
    df$Type <- sub("Cancer", "cancer", df$Type)
    df$Type.1 <- NULL

    png(filename = filename, width = 1600, height = 900)
    par(mar=c(8,8,4,2) + 0.1, cex.lab = 1.3, cex.axis = 1.3)
    boxplot(df$Methy_Percent ~ df$Type, las = 2)
    mtext(side = 2, line = 5, cex = 1.3,
            text = "Methylation percentage")
    mtext(side = 3, line = 1, cex = 2.0,
            text = paste("Methylation Percentage of Gene", group))
    dev.off()
}
