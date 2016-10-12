setwd("/home/yi/data/DNA_Methylation/LINE1/Data/GSE54210")
df <- read.table("GSE54210_methylated_unmethylated_signals.txt.txt", header = T, sep = "\t")
annot <- read.table("annot", header = T, sep = "\t")
methy_cutoff <- 0.7
unmethy_cutoff <- 0.3

beta.value <- function(methy, unmethy)
{
    temp <- methy / (methy + unmethy + 100)
    return(temp)
}

methy.percent <- function(df)
{
    class(df) <- "numeric"
    temp <- ifelse (df > methy_cutoff, 1, df)
    temp <- ifelse(temp < unmethy_cutoff, 0, temp)
    temp <- ifelse(temp == 1 | temp == 0, temp, 0.5)
    percentage <- colMeans(temp)
    return(percentage)
}

colnames(df) <- sub("(\\w{2})\\.(\\d*.*)", "\\1-\\2", colnames(df))
result <- matrix(nrow = nrow(df), ncol = (ncol(df) + 2) / 3)
colnames(result) <- c("ID", as.character(annot$ID))
result[ ,1] <- as.character(df[ ,1])
for (i in 2:(ncol(result)))
{
    temp <- df[ ,grep(colnames(result)[i], colnames(df))]
    x <- temp[ ,grep("Methylated", colnames(temp))]
    y <- temp[ ,grep("Unmethylated", colnames(temp))]
    result[ ,i] <- beta.value(x,y)
}
percentage <- methy.percent(result)

write.table(result, "./result.txt")
write.table(percentage, "./methy_percent.txt")